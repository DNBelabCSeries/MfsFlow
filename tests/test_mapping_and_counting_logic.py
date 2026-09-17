import os
import gzip
import subprocess
import tempfile
import unittest
import json
import collections
from types import SimpleNamespace
from unittest import mock

from mfsflow.scripts.mapping_analysis import (
    build_mapping_producer,
    build_star_command,
    build_star_misc_base,
    run_star_pipe,
    setup_gtf,
)
from mfsflow.scripts.dge_utils import balance_reference_chunks, summarize_exon_intron_counts
from mfsflow.scripts.read_utils import is_pair_representative
from mfsflow.scripts.run_featurecounts import (
    _config_bool,
    build_featurecounts_cmd,
    normalize_read_category,
    resolve_counting_strand_modes,
    run_parallel_fc_postprocessing,
    should_count_read,
    update_read_stats,
)
from mfsflow.path_layout import stage_state_dir
from mfsflow.stages.mapping import run_mapping_stage


class MappingAndCountingLogicTests(unittest.TestCase):
    def test_counting_parallel_option_parses_yaml_style_booleans(self):
        self.assertTrue(_config_bool("yes"))
        self.assertTrue(_config_bool(None, default=True))
        self.assertFalse(_config_bool("no", default=True))

    def test_parallel_fc_postprocessing_preserves_source_order_and_budget(self):
        created_processes = []

        class FakeProcess:
            def __init__(self, target, args, name):
                self.target = target
                self.args = args
                self.name = name
                self.exitcode = None
                self.pid = None
                created_processes.append(self)

            def start(self):
                self.pid = len(created_processes)
                source = self.args[6]
                with open(self.args[1], "wb") as handle:
                    handle.write(b"BAM")
                with open(self.args[2], "w", encoding="utf-8") as handle:
                    json.dump(
                        {
                            "read_stats": {source: {"Exon": 1}},
                            "coverage": [1] * 100,
                            "coverage_reads": 1,
                        },
                        handle,
                    )
                self.exitcode = 0

            def join(self):
                return None

            def is_alive(self):
                return False

            def terminate(self):
                self.exitcode = -15

        fake_context = SimpleNamespace(Process=FakeProcess)
        with tempfile.TemporaryDirectory() as tmpdir:
            internal = os.path.join(tmpdir, "internal.bam")
            umi = os.path.join(tmpdir, "umi.bam")
            for path in (internal, umi):
                with open(path, "wb") as handle:
                    handle.write(b"FC")

            with mock.patch(
                "mfsflow.scripts.run_featurecounts._parallel_postprocess_supported",
                return_value=True,
            ), mock.patch(
                "mfsflow.scripts.run_featurecounts.multiprocessing.get_context",
                return_value=fake_context,
            ):
                results = run_parallel_fc_postprocessing(
                    [("Internal", internal, 0), ("UMI", umi, 1)],
                    "samtools",
                    20,
                    {},
                    {},
                    True,
                    {},
                    {},
                    {"Internal": 0.5, "UMI": 0.25},
                    42,
                )

            self.assertEqual([result[0] for result in results], ["Internal", "UMI"])
            self.assertEqual([process.args[4] for process in created_processes], [4, 4])
            self.assertEqual(results[0][3]["Internal"]["Exon"], 1)
            self.assertFalse(os.path.exists(created_processes[0].args[2]))

    def test_parallel_fc_postprocessing_cleans_outputs_after_worker_failure(self):
        created_paths = []

        class FailingProcess:
            def __init__(self, target, args, name):
                self.args = args
                self.name = name
                self.exitcode = None
                self.pid = None

            def start(self):
                self.pid = 1
                for path in self.args[1:3]:
                    with open(path, "wb") as handle:
                        handle.write(b"partial")
                    created_paths.append(path)
                self.exitcode = 1 if self.name.endswith("umi") else 0

            def join(self):
                return None

            def is_alive(self):
                return False

            def terminate(self):
                self.exitcode = -15

        fake_context = SimpleNamespace(Process=FailingProcess)
        with tempfile.TemporaryDirectory() as tmpdir:
            internal = os.path.join(tmpdir, "internal.bam")
            umi = os.path.join(tmpdir, "umi.bam")
            for path in (internal, umi):
                with open(path, "wb") as handle:
                    handle.write(b"FC")

            with mock.patch(
                "mfsflow.scripts.run_featurecounts._parallel_postprocess_supported",
                return_value=True,
            ), mock.patch(
                "mfsflow.scripts.run_featurecounts.multiprocessing.get_context",
                return_value=fake_context,
            ):
                with self.assertRaisesRegex(RuntimeError, "UMI rc=1"):
                    run_parallel_fc_postprocessing(
                        [("Internal", internal, 0), ("UMI", umi, 1)],
                        "samtools",
                        8,
                        {},
                        {},
                        False,
                        {},
                        {},
                        {"Internal": 1.0, "UMI": 1.0},
                        42,
                    )

            self.assertTrue(created_paths)
            self.assertTrue(all(not os.path.exists(path) for path in created_paths))

    def test_build_star_misc_base_uses_target_specific_overhang(self):
        misc = build_star_misc_base(
            "/ref/star",
            8,
            "PE",
            False,
            "/tmp/genes.gtf",
            101,
            "samtools",
        )
        overhang_index = misc.index("--sjdbOverhang")
        self.assertEqual(misc[overhang_index + 1], "100")
        # stream_corrector outputs BAM; STAR must decode via samtools view
        command_index = misc.index("--readFilesCommand")
        self.assertEqual(misc[command_index + 1], "samtools view -@ 1")

    def test_pre_corrected_chunks_bypass_python_corrector(self):
        command, label = build_mapping_producer(
            ["chunk1.bam", "chunk2.bam"],
            False,
            "umi",
            "python3",
            "stream_corrector.py",
            "binning.tsv",
            "ids.tsv",
            "samtools",
        )
        self.assertEqual(command, ["samtools", "cat", "-o", "-", "chunk1.bam", "chunk2.bam"])
        self.assertEqual(label, "samtools cat")

    def test_raw_chunks_still_use_streaming_corrector(self):
        command, label = build_mapping_producer(
            ["chunk.raw.tagged.bam"],
            True,
            "internal",
            "python3",
            "stream_corrector.py",
            "binning.tsv",
            "ids.tsv",
            "samtools",
        )
        self.assertEqual(command[0:2], ["python3", "stream_corrector.py"])
        self.assertIn("internal", command)
        self.assertEqual(label, "stream_corrector")

    def test_build_star_misc_base_skips_overhang_when_index_has_sjdb(self):
        misc = build_star_misc_base(
            "/ref/star",
            8,
            "PE",
            True,
            "/tmp/genes.gtf",
            0,
            "samtools",
        )
        self.assertNotIn("--sjdbOverhang", misc)
        self.assertNotIn("--sjdbGTFfile", misc)

    def test_star_command_keeps_paths_with_spaces_as_single_arguments(self):
        command = build_star_command(
            "/tools with spaces/STAR",
            ["--genomeDir", "/reference with spaces/star"],
            "--clip3pAdapterSeq ACGT",
            True,
            "/output with spaces/sample.",
        )
        self.assertEqual(command[0], "/tools with spaces/STAR")
        self.assertEqual(command[2], "/reference with spaces/star")
        self.assertEqual(command[-1], "/output with spaces/sample.")
        self.assertIn("--twopassMode", command)

    def test_stream_corrector_failure_fails_mapping_even_when_star_succeeds(self):
        producer = mock.Mock(returncode=7)
        producer.stdout = mock.Mock()
        consumer = mock.Mock(returncode=0)
        with mock.patch(
            "mfsflow.scripts.mapping_analysis.subprocess.Popen",
            side_effect=[producer, consumer],
        ) as popen:
            with self.assertRaisesRegex(RuntimeError, "stream_corrector exited with code 7"):
                run_star_pipe(["python3", "stream_corrector.py"], ["STAR", "--readFilesIn", "/dev/stdin"])

        self.assertNotIn("shell", popen.call_args_list[1].kwargs)

    def test_mapping_timeout_kills_star_and_stream_corrector(self):
        producer = mock.Mock(returncode=None)
        producer.stdout = mock.Mock()
        producer.poll.return_value = None
        producer.wait.side_effect = [None]
        consumer = mock.Mock(returncode=None)
        consumer.poll.return_value = None
        consumer.wait.side_effect = [
            subprocess.TimeoutExpired(["STAR"], 0.1),
            None,
        ]
        with mock.patch(
            "mfsflow.scripts.mapping_analysis.subprocess.Popen",
            side_effect=[producer, consumer],
        ):
            with self.assertRaisesRegex(RuntimeError, "timed out"):
                run_star_pipe(
                    ["python3", "stream_corrector.py"],
                    ["STAR", "--readFilesIn", "/dev/stdin"],
                    timeout=0.1,
                )

        producer.kill.assert_called_once()
        consumer.kill.assert_called_once()

    def test_mapping_resume_uses_filtering_manifest_after_tmp_root_change(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = os.path.join(tmpdir, "XPRESS_PROCESSING")
            old_tmp = os.path.join(tmpdir, "old-tmp")
            new_tmp = os.path.join(tmpdir, "new-tmp")
            os.makedirs(old_tmp)
            os.makedirs(new_tmp)
            raw_bam = os.path.join(old_tmp, "sample.group_000.raw.tagged.bam")
            with open(raw_bam, "wb") as handle:
                handle.write(b"BAM")
            state_dir = stage_state_dir(outdir)
            os.makedirs(state_dir)
            with open(os.path.join(state_dir, "Filtering.manifest.json"), "w", encoding="utf-8") as handle:
                json.dump({"artifacts": [{"path": raw_bam}]}, handle)

            runtime = SimpleNamespace(
                project="sample",
                analysis_dir=outdir,
                tmp_merge_path=new_tmp,
                out_dir=outdir,
                python_exec="python3",
                yaml_file=os.path.join(outdir, "config", "run_config.yaml"),
                resolve_script=lambda name: name,
                which_stage="Mapping",
            )
            commands = []

            umi_chunks, int_chunks = run_mapping_stage(
                runtime,
                lambda command, _name: commands.append(command),
            )

            self.assertEqual(umi_chunks, [raw_bam])
            self.assertEqual(int_chunks, [raw_bam])
            self.assertIn(raw_bam, commands[0][commands[0].index("--umi_bam") + 1])

    def test_mapping_rerun_removes_stale_star_outputs_but_keeps_legacy_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            processing = os.path.join(tmpdir, "XPRESS_PROCESSING")
            os.makedirs(os.path.join(processing, "config"))
            stale = os.path.join(processing, "sample.filtered.tagged.internal.Aligned.out.bam")
            legacy = os.path.join(processing, "sample.filtered.tagged.internal.unmapped.bam")
            umi_chunk = os.path.join(tmpdir, "sample.raw.tagged.bam")
            for path in (stale, legacy, umi_chunk):
                with open(path, "wb") as handle:
                    handle.write(b"BAM")
            runtime = SimpleNamespace(
                project="sample",
                analysis_dir=processing,
                tmp_merge_path=os.path.join(processing, "intermediate", "tmp_merge"),
                out_dir=processing,
                python_exec="python3",
                yaml_file=os.path.join(processing, "config", "run_config.yaml"),
                resolve_script=lambda name: name,
                which_stage="Filtering",
            )

            run_mapping_stage(runtime, lambda _command, _name: None, umi_chunks=[umi_chunk])

            self.assertFalse(os.path.exists(stale))
            self.assertTrue(os.path.exists(legacy))

    def test_dge_reference_chunks_are_balanced_by_mapped_reads(self):
        chunks = balance_reference_chunks(
            ["chr1", "chr2", "chr3", "chr4"],
            {"chr1": 1000, "chr2": 600, "chr3": 400, "chr4": 1},
            2,
        )
        self.assertEqual(len(chunks), 2)
        loads = [sum({"chr1": 1000, "chr2": 600, "chr3": 400, "chr4": 1}[c] for c in chunk) for chunk in chunks]
        self.assertEqual(sorted(loads), [1000, 1001])

    def test_dge_cell_summary_preserves_existing_combined_stats_semantics(self):
        summary = summarize_exon_intron_counts(
            {"BC1": {"G1": 3, "G2": 2}},
            {"BC1": {"G1": 1, "G3": 4}},
        )
        self.assertEqual(summary["exon"]["BC1"], {"umis": 5, "genes": 2})
        self.assertEqual(summary["intron"]["BC1"], {"umis": 5, "genes": 2})
        self.assertEqual(summary["inex"]["BC1"], {"umis": 10, "genes": 3})

    def test_setup_gtf_decompresses_gzipped_gtf(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gz_gtf = os.path.join(tmpdir, "genes.gtf.gz")
            with gzip.open(gz_gtf, "wt") as handle:
                handle.write('chr1\tT\texon\t1\t10\t.\t+\t.\tgene_id "g1";\n')

            final_gtf, extra = setup_gtf(
                {"reference": {"GTF_file": gz_gtf, "STAR_index": tmpdir}},
                "Sample01",
                tmpdir,
                "samtools",
            )

            self.assertEqual(extra, "")
            self.assertTrue(final_gtf.endswith("Sample01.final_annot.gtf"))
            with open(final_gtf) as handle:
                self.assertIn('gene_id "g1"', handle.read())

    def test_featurecounts_pe_layout_enables_pair_options(self):
        cmd = build_featurecounts_cmd(
            "featureCounts",
            "in.bam",
            "features.saf",
            "out.txt",
            4,
            1,
            "PE",
        )
        self.assertIn("-p", cmd)
        self.assertIn("--countReadPairs", cmd)
        self.assertIn("-C", cmd)

    def test_featurecounts_se_layout_does_not_enable_pair_options(self):
        cmd = build_featurecounts_cmd(
            "featureCounts",
            "in.bam",
            "features.saf",
            "out.txt",
            4,
            1,
            "SE",
        )
        self.assertNotIn("-p", cmd)
        self.assertNotIn("--countReadPairs", cmd)
        self.assertNotIn("-C", cmd)

    def test_resolve_counting_strand_modes_preserves_historical_internal_default(self):
        umi_mode, internal_mode = resolve_counting_strand_modes({"strand": 2})
        self.assertEqual((umi_mode, internal_mode), (2, 0))

    def test_resolve_counting_strand_modes_supports_explicit_internal_override(self):
        umi_mode, internal_mode = resolve_counting_strand_modes({"strand": 1, "internal_strand": 0})
        self.assertEqual((umi_mode, internal_mode), (1, 0))

    def test_should_count_read_filters_secondary_and_supplementary(self):
        self.assertTrue(should_count_read("q1\t99\tchr1\t1\t255\t50M\t=\t1\t0\tACGT\tFFFF"))
        self.assertFalse(should_count_read("q2\t355\tchr1\t1\t255\t50M\t=\t1\t0\tACGT\tFFFF"))
        self.assertFalse(should_count_read("q3\t2147\tchr1\t1\t255\t50M\t=\t1\t0\tACGT\tFFFF"))

    def test_pair_representative_counts_r1_once(self):
        r1 = "q1\t99\tchr1\t1\t255\t50M\t=\t1\t0\tACGT\tFFFF"
        r2 = "q1\t147\tchr1\t1\t255\t50M\t=\t1\t0\tACGT\tFFFF"
        self.assertTrue(is_pair_representative(r1))
        self.assertFalse(is_pair_representative(r2))

    def test_pair_representative_keeps_unmapped_r1_like_zumis(self):
        r1 = "q1\t69\t*\t0\t0\t*\t=\t1\t0\tACGT\tFFFF"
        r2 = "q1\t137\tchr1\t1\t255\t50M\t=\t0\t0\tACGT\tFFFF"
        self.assertTrue(is_pair_representative(r1))
        self.assertFalse(is_pair_representative(r2))

    def test_pair_representative_rejects_secondary_and_supplementary(self):
        secondary_r1 = "q1\t355\tchr1\t1\t255\t50M\t=\t1\t0\tACGT\tFFFF"
        supplementary_r1 = "q1\t2147\tchr1\t1\t255\t50M\t=\t1\t0\tACGT\tFFFF"
        self.assertFalse(is_pair_representative(secondary_r1))
        self.assertFalse(is_pair_representative(supplementary_r1))

    def test_single_end_record_is_always_representative(self):
        self.assertTrue(is_pair_representative("q1\t0\tchr1\t1\t255\t50M\t*\t0\t0\tACGT\tFFFF"))

    def test_normalize_read_category_collapses_unassigned_reasons(self):
        self.assertEqual(normalize_read_category("Exon"), "Exon")
        self.assertEqual(normalize_read_category("MappingQuality"), "Other_Unassigned")
        self.assertEqual(normalize_read_category("FragmentLength"), "Other_Unassigned")

    def test_sam_fallback_counts_missing_cb_as_unused(self):
        stats = collections.defaultdict(lambda: collections.defaultdict(int))
        line = "q1\t64\tchr1\t1\t255\t50M\t*\t0\t0\tACGT\tFFFF"
        update_read_stats(stats, line, "Exon", "UMI")
        self.assertEqual(1, stats["__NO_CB__"]["Unused BC"])


if __name__ == "__main__":
    unittest.main()

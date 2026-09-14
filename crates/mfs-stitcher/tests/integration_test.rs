use rust_htslib::bam::header::Header;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{self, Read};
use std::fs::File;
use std::io::Write;

use mfs_stitcher::pipeline::{run_pipeline, PipelineConfig};

#[test]
fn test_end_to_end_stitching() {
    let tmp = std::env::temp_dir().join(format!("mfs_stitcher_test_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).expect("create tempdir");
    let in_bam_path = tmp.join("input.bam");
    let out_bam_path = tmp.join("stitched.bam");
    let gtf_path = tmp.join("annotation.gtf");

    // 1. Create minimal BAM header
    let mut header = Header::new();
    let mut header_rec = rust_htslib::bam::header::HeaderRecord::new(b"SQ");
    header_rec.push_tag(b"SN", "chr1");
    header_rec.push_tag(b"LN", "10000");
    header.push_record(&header_rec);

    // 2. Write input BAM with 2 reads for UMI1 (overlapping), 1 read for
    // UMI2, and a long-spliced read whose stitched position precedes the
    // GTF interval of its later gene.  The latter catches per-cluster-only
    // sorting: GENE3 is processed after GENE2 but its molecule starts at
    // position 2500, before GENE2's molecule at position 3000.
    {
        let mut writer = bam::Writer::from_path(&in_bam_path, &header, bam::Format::Bam)
            .expect("open bam writer");

        // Read 1 for CellA, UMI1: pos 1000, 50M
        let mut r1 = bam::Record::new();
        r1.set(
            b"read1",
            Some(&CigarString(vec![Cigar::Match(50)])),
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC",
            &[30u8; 50],
        );
        r1.set_tid(0);
        r1.set_pos(1000);
        r1.set_flags(0); // single-end / mapped
        r1.push_aux(b"CB", Aux::String("CellA")).unwrap();
        r1.push_aux(b"UB", Aux::String("UMI1")).unwrap();
        r1.push_aux(b"GE", Aux::String("GENE1")).unwrap();
        writer.write(&r1).unwrap();

        // Read 2 for CellA, UMI1: pos 1020, 50M (overlaps 1020-1049, extends to 1069)
        let mut r2 = bam::Record::new();
        r2.set(
            b"read2",
            Some(&CigarString(vec![Cigar::Match(50)])),
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC",
            &[35u8; 50],
        );
        r2.set_tid(0);
        r2.set_pos(1020);
        r2.set_flags(0);
        r2.push_aux(b"CB", Aux::String("CellA")).unwrap();
        r2.push_aux(b"UB", Aux::String("UMI1")).unwrap();
        r2.push_aux(b"GE", Aux::String("GENE1")).unwrap();
        writer.write(&r2).unwrap();

        // Read 3 for CellB, UMI2: pos 1500, 30M
        let mut r3 = bam::Record::new();
        r3.set(
            b"read3",
            Some(&CigarString(vec![Cigar::Match(30)])),
            b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            &[40u8; 30],
        );
        r3.set_tid(0);
        r3.set_pos(1500);
        r3.set_flags(0);
        r3.push_aux(b"CB", Aux::String("CellB")).unwrap();
        r3.push_aux(b"UB", Aux::String("UMI2")).unwrap();
        r3.push_aux(b"GE", Aux::String("GENE1")).unwrap();
        writer.write(&r3).unwrap();

        // Read 4 for a later gene with a long reference skip. It overlaps
        // GENE3's 5000-6000 fetch interval but starts at 2500.
        let mut r4 = bam::Record::new();
        r4.set(
            b"read4",
            Some(&CigarString(vec![
                Cigar::Match(30),
                Cigar::RefSkip(2500),
                Cigar::Match(30),
            ])),
            &[b'G'; 60],
            &[40u8; 60],
        );
        r4.set_tid(0);
        r4.set_pos(2500);
        r4.set_flags(0);
        r4.push_aux(b"CB", Aux::String("CellD")).unwrap();
        r4.push_aux(b"UB", Aux::String("UMI4")).unwrap();
        r4.push_aux(b"GE", Aux::String("GENE3")).unwrap();
        writer.write(&r4).unwrap();

        // Read 5 for a second, non-overlapping gene.
        let mut r5 = bam::Record::new();
        r5.set(
            b"read5",
            Some(&CigarString(vec![Cigar::Match(30)])),
            b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            &[40u8; 30],
        );
        r5.set_tid(0);
        r5.set_pos(3000);
        r5.set_flags(0);
        r5.push_aux(b"CB", Aux::String("CellC")).unwrap();
        r5.push_aux(b"UB", Aux::String("UMI3")).unwrap();
        r5.push_aux(b"GE", Aux::String("GENE2")).unwrap();
        writer.write(&r5).unwrap();
    }

    // Build BAM index
    bam::index::build(&in_bam_path, None, bam::index::Type::Bai, 1).expect("build bai index");

    // 3. Write minimal GTF with gene, transcript, and exons
    {
        let mut gtf_file = File::create(&gtf_path).unwrap();
        writeln!(
            gtf_file,
            "chr1\tensembl\tgene\t3000\t4000\t.\t+\t.\tgene_id \"GENE2\"; gene_name \"GENE2\";\n\
             chr1\tensembl\ttranscript\t3000\t4000\t.\t+\t.\tgene_id \"GENE2\"; transcript_id \"TRANSCRIPT2\";\n\
             chr1\tensembl\texon\t3000\t4000\t.\t+\t.\tgene_id \"GENE2\"; transcript_id \"TRANSCRIPT2\";\n\
             chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n\
             chr1\tensembl\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TRANSCRIPT1\";\n\
             chr1\tensembl\texon\t1000\t1200\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TRANSCRIPT1\";\n\
             chr1\tensembl\texon\t1400\t2000\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TRANSCRIPT1\";\n\
             chr1\tensembl\tgene\t5000\t6000\t.\t+\t.\tgene_id \"GENE3\"; gene_name \"GENE3\";\n\
             chr1\tensembl\ttranscript\t5000\t6000\t.\t+\t.\tgene_id \"GENE3\"; transcript_id \"TRANSCRIPT3\";\n\
             chr1\tensembl\texon\t5000\t6000\t.\t+\t.\tgene_id \"GENE3\"; transcript_id \"TRANSCRIPT3\";"
        )
        .unwrap();
    }

    // 4. Run pipeline with zero-friction on-the-fly indexing, index dumping, and matrix export
    let dump_prefix = tmp.join("test_index");
    let mex_dir = tmp.join("mex_out");
    let counts_tsv = tmp.join("counts.tsv.gz");
    let molecules_tsv = tmp.join("molecules.tsv.gz");

    let config = PipelineConfig {
        input_bam: in_bam_path,
        output_bam: out_bam_path.clone(),
        gtf_file: gtf_path,
        isoform_file: None,
        junction_file: None,
        dump_index: Some(dump_prefix.clone()),
        matrix_out: Some(mex_dir.clone()),
        counts_tsv: Some(counts_tsv.clone()),
        molecules_tsv: Some(molecules_tsv.clone()),
        threads: 2,
        single_end: true,
        skip_iso: false, // tests zero-friction on-the-fly index calculation
        umi_tag: "UB".to_string(),
        cell_tag: "CB".to_string(),
        cells_file: None,
        genes_file: None,
        contig: None,
        gene_identifier: "gene_id".to_string(),
    };

    run_pipeline(config).expect("pipeline run successful");

    // 5. Inspect stitched output BAM
    let mut reader = bam::Reader::from_path(&out_bam_path).expect("open stitched bam");
    let mut stitched_records = Vec::new();
    let mut rec = bam::Record::new();
    while let Some(Ok(())) = reader.read(&mut rec) {
        stitched_records.push(rec.clone());
    }

    assert_eq!(stitched_records.len(), 4, "Expected 4 stitched molecules");

    for pair in stitched_records.windows(2) {
        assert!(
            pair[0].tid() < pair[1].tid()
                || (pair[0].tid() == pair[1].tid() && pair[0].pos() <= pair[1].pos()),
            "Output BAM is not coordinate sorted: ({}, {}) followed by ({}, {})",
            pair[0].tid(),
            pair[0].pos(),
            pair[1].tid(),
            pair[1].pos()
        );
    }

    let qnames: Vec<String> = stitched_records
        .iter()
        .map(|r| String::from_utf8_lossy(r.qname()).to_string())
        .collect();
    assert!(qnames.contains(&"CellA:GENE1:UMI1".to_string()));
    assert!(qnames.contains(&"CellB:GENE1:UMI2".to_string()));
    assert!(qnames.contains(&"CellC:GENE2:UMI3".to_string()));
    assert!(qnames.contains(&"CellD:GENE3:UMI4".to_string()));

    // Find CellA molecule: it merged 1000..1049 and 1020..1069 -> 1000..1069 (70 bp total)
    let cell_a_rec = stitched_records
        .iter()
        .find(|r| r.qname() == b"CellA:GENE1:UMI1")
        .unwrap();
    assert_eq!(cell_a_rec.pos(), 1000);
    assert_eq!(cell_a_rec.seq_len(), 70);
    assert_eq!(cell_a_rec.cigar().to_string(), "70M");

    // Verify tags
    assert_eq!(cell_a_rec.aux(b"NR").unwrap(), Aux::I32(2)); // 2 reads
    assert_eq!(cell_a_rec.aux(b"ER").unwrap(), Aux::I32(2)); // 2 exonic
    assert_eq!(cell_a_rec.aux(b"CB").unwrap(), Aux::String("CellA"));
    assert_eq!(cell_a_rec.aux(b"UB").unwrap(), Aux::String("UMI1"));
    assert_eq!(cell_a_rec.aux(b"GX").unwrap(), Aux::String("GENE1"));

    // Verify CT tag generated by zero-friction on-the-fly indexing
    assert_eq!(cell_a_rec.aux(b"CT").unwrap(), Aux::String("TRANSCRIPT1"));

    // Verify index dumping produced expected files
    let iso_file = format!("{}.intervals.json.gz", dump_prefix.display());
    let jun_file = format!("{}.refskip.json.gz", dump_prefix.display());
    assert!(
        std::path::Path::new(&iso_file).exists(),
        "Intervals JSON file should exist"
    );
    assert!(
        std::path::Path::new(&jun_file).exists(),
        "Refskip JSON file should exist"
    );

    // Verify MEX export files
    assert!(
        mex_dir.join("matrix.mtx.gz").exists(),
        "matrix.mtx.gz should exist"
    );
    assert!(
        mex_dir.join("barcodes.tsv.gz").exists(),
        "barcodes.tsv.gz should exist"
    );
    assert!(
        mex_dir.join("features.tsv.gz").exists(),
        "features.tsv.gz should exist"
    );

    // Verify TSV export files
    assert!(counts_tsv.exists(), "counts.tsv.gz should exist");
    assert!(molecules_tsv.exists(), "molecules.tsv.gz should exist");

    // Verify output BAM index was created automatically and is readable
    let bai_path = format!("{}.bai", out_bam_path.display());
    assert!(
        std::path::Path::new(&bai_path).exists(),
        "Output BAM index (.bai) should exist"
    );
    let mut indexed_reader =
        bam::IndexedReader::from_path(&out_bam_path).expect("open indexed output bam");
    indexed_reader
        .fetch((0, 1000, 2000))
        .expect("fetch region from output bam");
}

#[test]
fn test_edge_cases_and_fixes() {
    let tmp = std::env::temp_dir().join(format!("mfs_stitcher_edge_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).expect("create tempdir");
    let in_bam_path = tmp.join("input.bam");
    let out_bam_path = tmp.join("stitched.bam");
    let gtf_path = tmp.join("annotation.gtf");

    // 1. Header with 2 chromosomes: chr1 and chr2
    let mut header = Header::new();
    let mut sq1 = rust_htslib::bam::header::HeaderRecord::new(b"SQ");
    sq1.push_tag(b"SN", "chr1");
    sq1.push_tag(b"LN", "10000");
    header.push_record(&sq1);

    let mut sq2 = rust_htslib::bam::header::HeaderRecord::new(b"SQ");
    sq2.push_tag(b"SN", "chr2");
    sq2.push_tag(b"LN", "10000");
    header.push_record(&sq2);

    {
        let mut writer = bam::Writer::from_path(&in_bam_path, &header, bam::Format::Bam)
            .expect("open bam writer");

        // --- Edge Case 1: Soft-clipped read (10S30M) on chr1 at pos 500 ---
        // Query seq: 10 'T's (soft-clipped) + 30 'A's (matching reference 500..529)
        let mut r_soft = bam::Record::new();
        let mut soft_seq = Vec::new();
        soft_seq.extend(vec![b'T'; 10]);
        soft_seq.extend(vec![b'A'; 30]);
        r_soft.set(
            b"read_soft",
            Some(&CigarString(vec![Cigar::SoftClip(10), Cigar::Match(30)])),
            &soft_seq,
            &[30u8; 40],
        );
        r_soft.set_tid(0); // chr1
        r_soft.set_pos(500);
        r_soft.set_flags(0); // mapped
        r_soft.push_aux(b"CB", Aux::String("Cell1")).unwrap();
        r_soft.push_aux(b"UB", Aux::String("UMI_SOFT")).unwrap();
        r_soft.push_aux(b"GE", Aux::String("GENE1")).unwrap();
        writer.write(&r_soft).unwrap();

        // --- Edge Case 2: Insertion read (10M2I20M) on chr1 at pos 600 ---
        // Query seq: 10 'A's (ref 600..609) + 2 'C's (insertion) + 20 'G's (ref 610..629)
        let mut r_ins = bam::Record::new();
        let mut ins_seq = Vec::new();
        ins_seq.extend(vec![b'A'; 10]);
        ins_seq.extend(vec![b'C'; 2]);
        ins_seq.extend(vec![b'G'; 20]);
        r_ins.set(
            b"read_ins",
            Some(&CigarString(vec![
                Cigar::Match(10),
                Cigar::Ins(2),
                Cigar::Match(20),
            ])),
            &ins_seq,
            &[30u8; 32],
        );
        r_ins.set_tid(0); // chr1
        r_ins.set_pos(600);
        r_ins.set_flags(0);
        r_ins.push_aux(b"CB", Aux::String("Cell1")).unwrap();
        r_ins.push_aux(b"UB", Aux::String("UMI_INS")).unwrap();
        r_ins.push_aux(b"GE", Aux::String("GENE1")).unwrap();
        writer.write(&r_ins).unwrap();

        // --- Edge Case 3: Secondary (0x100) and Supplementary (0x800) filtering ---
        // Read with UMI_FILTER: 1 primary read + 1 secondary read
        let mut r_primary = bam::Record::new();
        r_primary.set(
            b"read_prim",
            Some(&CigarString(vec![Cigar::Match(30)])),
            &[b'C'; 30],
            &[30u8; 30],
        );
        r_primary.set_tid(0);
        r_primary.set_pos(650);
        r_primary.set_flags(0); // primary
        r_primary.push_aux(b"CB", Aux::String("Cell1")).unwrap();
        r_primary
            .push_aux(b"UB", Aux::String("UMI_FILTER"))
            .unwrap();
        r_primary.push_aux(b"GE", Aux::String("GENE1")).unwrap();
        writer.write(&r_primary).unwrap();

        let mut r_sec = bam::Record::new();
        r_sec.set(
            b"read_sec",
            Some(&CigarString(vec![Cigar::Match(30)])),
            &[b'T'; 30],
            &[30u8; 30],
        );
        r_sec.set_tid(0);
        r_sec.set_pos(650);
        r_sec.set_flags(0x100); // secondary alignment flag!
        r_sec.push_aux(b"CB", Aux::String("Cell1")).unwrap();
        r_sec.push_aux(b"UB", Aux::String("UMI_FILTER")).unwrap();
        r_sec.push_aux(b"GE", Aux::String("GENE1")).unwrap();
        writer.write(&r_sec).unwrap();

        // --- Edge Case 4: Read on chr2 at pos 100 for GENE2 ---
        let mut r_chr2 = bam::Record::new();
        r_chr2.set(
            b"read_chr2",
            Some(&CigarString(vec![Cigar::Match(40)])),
            &[b'G'; 40],
            &[35u8; 40],
        );
        r_chr2.set_tid(1); // chr2
        r_chr2.set_pos(100);
        r_chr2.set_flags(0);
        r_chr2.push_aux(b"CB", Aux::String("Cell2")).unwrap();
        r_chr2.push_aux(b"UB", Aux::String("UMI_CHR2")).unwrap();
        r_chr2.push_aux(b"GE", Aux::String("GENE2")).unwrap();
        writer.write(&r_chr2).unwrap();
    }

    bam::index::build(&in_bam_path, None, bam::index::Type::Bai, 1).expect("build bai index");

    // Write GTF for GENE1 (chr1) and GENE2 (chr2)
    {
        let mut gtf_file = File::create(&gtf_path).unwrap();
        writeln!(
            gtf_file,
            "chr1\tensembl\tgene\t400\t2000\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n\
             chr1\tensembl\ttranscript\t400\t2000\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\";\n\
             chr1\tensembl\texon\t400\t2000\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\";\n\
             chr2\tensembl\tgene\t50\t1000\t.\t+\t.\tgene_id \"GENE2\"; gene_name \"GENE2\";\n\
             chr2\tensembl\ttranscript\t50\t1000\t.\t+\t.\tgene_id \"GENE2\"; transcript_id \"TX2\";\n\
             chr2\tensembl\texon\t50\t1000\t.\t+\t.\tgene_id \"GENE2\"; transcript_id \"TX2\";"
        )
        .unwrap();
    }

    let config = PipelineConfig {
        input_bam: in_bam_path,
        output_bam: out_bam_path.clone(),
        gtf_file: gtf_path,
        isoform_file: None,
        junction_file: None,
        dump_index: None,
        matrix_out: None,
        counts_tsv: None,
        molecules_tsv: None,
        threads: 2,
        single_end: true,
        skip_iso: false,
        umi_tag: "UB".to_string(),
        cell_tag: "CB".to_string(),
        cells_file: None,
        genes_file: None,
        contig: None,
        gene_identifier: "gene_id".to_string(),
    };

    run_pipeline(config).expect("run pipeline for edge cases");

    // Open output BAM and check records
    let mut reader = bam::Reader::from_path(&out_bam_path).expect("open output bam");
    let mut recs = Vec::new();
    let mut rec = bam::Record::new();
    while let Some(Ok(())) = reader.read(&mut rec) {
        recs.push(rec.clone());
    }

    assert_eq!(
        recs.len(),
        4,
        "Expected 4 molecules (UMI_SOFT, UMI_INS, UMI_FILTER, UMI_CHR2)"
    );

    // Check that records are strictly sorted by tid and pos
    for i in 0..recs.len() - 1 {
        let tid_a = recs[i].tid();
        let tid_b = recs[i + 1].tid();
        let pos_a = recs[i].pos();
        let pos_b = recs[i + 1].pos();
        assert!(
            tid_a < tid_b || (tid_a == tid_b && pos_a <= pos_b),
            "Output BAM not sorted! rec[{}] = ({}, {}), rec[{}] = ({}, {})",
            i,
            tid_a,
            pos_a,
            i + 1,
            tid_b,
            pos_b
        );
    }

    // 1. Verify Soft-clipped molecule: pos 500, seq should be 30 'A's (NOT 10 'T's!)
    let soft_mol = recs
        .iter()
        .find(|r| r.qname() == b"Cell1:GENE1:UMI_SOFT")
        .unwrap();
    assert_eq!(soft_mol.pos(), 500);
    assert_eq!(soft_mol.seq_len(), 30);
    let expected_soft_seq = "A".repeat(30);
    assert_eq!(
        String::from_utf8_lossy(&soft_mol.seq().as_bytes()),
        expected_soft_seq,
        "Soft-clipped read bases must not shift consensus sequence"
    );

    // 2. Verify Insertion molecule: pos 600, seq should be 10 'A's + 20 'G's (30 bp)
    let ins_mol = recs
        .iter()
        .find(|r| r.qname() == b"Cell1:GENE1:UMI_INS")
        .unwrap();
    assert_eq!(ins_mol.pos(), 600);
    assert_eq!(ins_mol.seq_len(), 30);
    let expected_ins_seq = format!("{}{}", "A".repeat(10), "G".repeat(20));
    assert_eq!(
        String::from_utf8_lossy(&ins_mol.seq().as_bytes()),
        expected_ins_seq,
        "Insertion must not displace downstream reference coordinate consensus"
    );

    // 3. Verify Secondary alignment filtering: UMI_FILTER should have NR == 1 (secondary ignored)
    let filter_mol = recs
        .iter()
        .find(|r| r.qname() == b"Cell1:GENE1:UMI_FILTER")
        .unwrap();
    assert_eq!(filter_mol.aux(b"NR").unwrap(), Aux::I32(1));
    assert_eq!(
        String::from_utf8_lossy(&filter_mol.seq().as_bytes()),
        "C".repeat(30)
    );

    // 4. Verify chr2 molecule exists and has tid == 1
    let chr2_mol = recs
        .iter()
        .find(|r| r.qname() == b"Cell2:GENE2:UMI_CHR2")
        .unwrap();
    assert_eq!(chr2_mol.tid(), 1);
    assert_eq!(chr2_mol.pos(), 100);

    // 5. Verify IndexedReader works across both chromosomes on the output BAM
    let mut indexed =
        bam::IndexedReader::from_path(&out_bam_path).expect("open indexed stitched bam");
    indexed.fetch((0, 450, 700)).expect("fetch chr1 region");
    let mut chr1_fetched = 0;
    while let Some(Ok(_)) = indexed.read(&mut rec) {
        chr1_fetched += 1;
    }
    assert_eq!(chr1_fetched, 3); // UMI_SOFT, UMI_INS, UMI_FILTER

    indexed.fetch((1, 50, 200)).expect("fetch chr2 region");
    let mut chr2_fetched = 0;
    while let Some(Ok(_)) = indexed.read(&mut rec) {
        chr2_fetched += 1;
    }
    assert_eq!(chr2_fetched, 1); // UMI_CHR2
}

#[test]
fn test_paired_end_dedup_and_filtering() {
    let tmp = std::env::temp_dir().join(format!("mfs_stitcher_pe_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).expect("create tempdir");
    let in_bam_path = tmp.join("input_pe.bam");
    let out_bam_path = tmp.join("stitched_pe.bam");
    let gtf_path = tmp.join("annotation_pe.gtf");

    let mut header = Header::new();
    let mut sq1 = rust_htslib::bam::header::HeaderRecord::new(b"SQ");
    sq1.push_tag(b"SN", "chr1");
    sq1.push_tag(b"LN", "10000");
    header.push_record(&sq1);

    {
        let mut writer = bam::Writer::from_path(&in_bam_path, &header, bam::Format::Bam)
            .expect("open bam writer");

        // 1. Valid PE pair: QNAME "frag1", R1 and R2
        let mut r1 = bam::Record::new();
        r1.set(
            b"frag1",
            Some(&CigarString(vec![Cigar::Match(50)])),
            &[b'A'; 50],
            &[30u8; 50],
        );
        r1.set_tid(0);
        r1.set_pos(1000);
        r1.set_flags(0x43); // Paired (0x1), Proper pair (0x2), First in template (0x40)
        r1.push_aux(b"CB", Aux::String("CellPE")).unwrap();
        r1.push_aux(b"UB", Aux::String("UMI_VALID")).unwrap();
        r1.push_aux(b"GE", Aux::String("GENE_PE")).unwrap();
        writer.write(&r1).unwrap();

        let mut r2 = bam::Record::new();
        r2.set(
            b"frag1",
            Some(&CigarString(vec![Cigar::Match(50)])),
            &[b'A'; 50],
            &[30u8; 50],
        );
        r2.set_tid(0);
        r2.set_pos(1020);
        r2.set_flags(0x83); // Paired (0x1), Proper pair (0x2), Second in template (0x80)
        r2.push_aux(b"CB", Aux::String("CellPE")).unwrap();
        r2.push_aux(b"UB", Aux::String("UMI_VALID")).unwrap();
        r2.push_aux(b"GE", Aux::String("GENE_PE")).unwrap();
        writer.write(&r2).unwrap();

        // 2. Mate unmapped read: should be filtered in PE mode
        let mut r_mate_unmapped = bam::Record::new();
        r_mate_unmapped.set(
            b"frag_bad1",
            Some(&CigarString(vec![Cigar::Match(50)])),
            &[b'C'; 50],
            &[30u8; 50],
        );
        r_mate_unmapped.set_tid(0);
        r_mate_unmapped.set_pos(1100);
        r_mate_unmapped.set_flags(0x49); // Paired (0x1), Mate unmapped (0x8), First (0x40)
        r_mate_unmapped
            .push_aux(b"CB", Aux::String("CellPE"))
            .unwrap();
        r_mate_unmapped
            .push_aux(b"UB", Aux::String("UMI_UNMAPPED_MATE"))
            .unwrap();
        r_mate_unmapped
            .push_aux(b"GE", Aux::String("GENE_PE"))
            .unwrap();
        writer.write(&r_mate_unmapped).unwrap();

        // 3. Improper pair: should be filtered in PE mode
        let mut r_improper = bam::Record::new();
        r_improper.set(
            b"frag_bad2",
            Some(&CigarString(vec![Cigar::Match(50)])),
            &[b'G'; 50],
            &[30u8; 50],
        );
        r_improper.set_tid(0);
        r_improper.set_pos(1200);
        r_improper.set_flags(0x41); // Paired (0x1), First (0x40), NOT proper pair
        r_improper.push_aux(b"CB", Aux::String("CellPE")).unwrap();
        r_improper
            .push_aux(b"UB", Aux::String("UMI_IMPROPER"))
            .unwrap();
        r_improper.push_aux(b"GE", Aux::String("GENE_PE")).unwrap();
        writer.write(&r_improper).unwrap();
    }

    bam::index::build(&in_bam_path, None, bam::index::Type::Bai, 1).expect("build bai index");

    {
        let mut gtf_file = File::create(&gtf_path).unwrap();
        writeln!(
            gtf_file,
            "chr1\tensembl\tgene\t900\t2000\t.\t+\t.\tgene_id \"GENE_PE\"; gene_name \"GENE_PE\";\n\
             chr1\tensembl\ttranscript\t900\t2000\t.\t+\t.\tgene_id \"GENE_PE\"; transcript_id \"TX_PE\";\n\
             chr1\tensembl\texon\t900\t2000\t.\t+\t.\tgene_id \"GENE_PE\"; transcript_id \"TX_PE\";"
        )
        .unwrap();
    }

    let config = PipelineConfig {
        input_bam: in_bam_path,
        output_bam: out_bam_path.clone(),
        gtf_file: gtf_path,
        isoform_file: None,
        junction_file: None,
        dump_index: None,
        matrix_out: None,
        counts_tsv: None,
        molecules_tsv: None,
        threads: 2,
        single_end: false, // PE mode!
        skip_iso: false,
        umi_tag: "UB".to_string(),
        cell_tag: "CB".to_string(),
        cells_file: None,
        genes_file: None,
        contig: None,
        gene_identifier: "gene_id".to_string(),
    };

    run_pipeline(config).expect("run PE pipeline");

    let mut reader = bam::Reader::from_path(&out_bam_path).expect("open stitched PE bam");
    let mut recs = Vec::new();
    let mut rec = bam::Record::new();
    while let Some(Ok(())) = reader.read(&mut rec) {
        recs.push(rec.clone());
    }

    // Only the valid pair (UMI_VALID) should be stitched. UMI_UNMAPPED_MATE and UMI_IMPROPER must be discarded.
    assert_eq!(
        recs.len(),
        1,
        "Only 1 valid molecule should be produced in PE mode"
    );

    let valid_rec = &recs[0];
    assert_eq!(valid_rec.qname(), b"CellPE:GENE_PE:UMI_VALID");

    // Check that PE reads sharing the same QNAME are deduplicated so NR == 1 and ER == 1!
    assert_eq!(
        valid_rec.aux(b"NR").unwrap(),
        Aux::I32(1),
        "NR must be 1 (deduplicated paired-end fragment/template count), NOT 2!"
    );
    assert_eq!(
        valid_rec.aux(b"ER").unwrap(),
        Aux::I32(1),
        "ER must be 1 (deduplicated paired-end fragment count), NOT 2!"
    );
}

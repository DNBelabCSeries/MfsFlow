"""
HTML report generation: assembles multi-omics analysis metrics and renders interactive Plotly-based reports.

This module reads pipeline output files, computes summary statistics, and
generates interactive HTML reports with embedded Plotly charts for
visualization of sequencing quality, barcode statistics, and expression data.
"""

import json
import csv
import statistics
import gzip
import hashlib
import html
import logging
import shutil
import base64
import os
from datetime import datetime
from pathlib import Path
from string import Template
import re

from mfsflow import __version__
from mfsflow.path_layout import config_dir, expression_dir, outputs_dir, stats_dir

logger = logging.getLogger(__name__)

_JS_TEMPLATE_PLACEHOLDERS = {
    "id",
    "plotId",
    "clusterNum",
    "libraryId",
    "coef",
    "exp",
    "cluster",
    "emptyMessage",
    "plate",
    "titleMetric",
}

_JSON_ARRAY_PLACEHOLDERS = {
    "bead_count_x",
    "bead_count_y",
    "rna_umap_x",
    "rna_umap_y",
    "rna_leiden",
    "rna_numi",
    "rna_clusters",
    "rna_gene_body_percentile",
    "rna_gene_body_umi_maxnorm",
    "rna_gene_body_internal_maxnorm",
    "rna_gene_body_all_maxnorm",
    "rna_saturation_fraction",
    "rna_saturation_lib_pct",
    "rna_saturation_gene_pct",
    "rna_saturation_median_genes_umi",
    "rna_saturation_median_genes_read",
    "rna_saturation_sampling_fraction",
    "rna_saturation_seq_pct",
    "rna_saturation_median_genes",
    "rna_rankplot_data",
    "rna_stats_table_data",
    "rna_marker_data",
    "barcode_mode_cards_data",
    "manual_barcode_table_data",
    "auto_plate_summary_data",
    "sequencing_quality_summary_data",
}

_JSON_OBJECT_PLACEHOLDERS = {
    "rna_read_distribution_bar_data",
    "rna_read_distribution_box_data",
    "barcode_report_summary_data",
    "well_qc_status_data",
}

_CANONICAL_STATS_KEYS = {
    "exon_reads": "Exon_reads",
    "intron_reads": "intron_reads",
    "intergenic_reads": "Intergenic_reads",
    "ambiguity_reads": "Ambiguity_reads",
    "unmapped_reads": "Unmapped_reads",
    "antisense_reads": "Antisense_reads",
    "exon_genes": "Exon_genes",
    "intron_genes": "Intron_genes",
    "intron_exon_genes": "Intron_Exon_genes",
    "exon_read_genes": "Exon_read_genes",
    "intron_read_genes": "Intron_read_genes",
    "intron_exon_read_genes": "Intron_Exon_read_genes",
    "exon_umis": "Exon_umis",
    "intron_umis": "Intron_umis",
    "intron_exon_umis": "Intron_Exon_umis",
    "umi_reads": "umi_reads",
    "internal_reads": "internal_reads",
    "wellid": "wellID",
}


def _default_placeholder_value(name):
    if name in _JSON_ARRAY_PLACEHOLDERS:
        return "[]"
    if name in _JSON_OBJECT_PLACEHOLDERS:
        return "{}"
    if name.endswith("_enabled") or name.endswith("_available"):
        return "false"
    if name.endswith("_html"):
        return ""
    if name.endswith("_data"):
        return "[]"
    if name.endswith("_pct") or name.endswith("_percent") or name.endswith("_percentage"):
        return ""
    if name == "rna_species":
        return "Unknown"
    if name == "rna_fraction_transcriptome_in_cells":
        return ""
    if name == "rna_saturation":
        return ""
    return ""


def _load_logo_data_uri(template_dir):
    for name in ("logo.png", "logo.svg", "logo.jpg", "logo.jpeg", "logo.ico"):
        path = template_dir / name
        if not path.exists():
            continue
        ext = path.suffix.lower().lstrip(".")
        if ext == "svg" or ext == "ico":
            try:
                raw = path.read_bytes()
                mime = f"image/{ext}"
                return f"data:{mime};base64,{base64.b64encode(raw).decode('ascii')}"
            except Exception:
                pass
            continue
        compressed = _compress_logo_image(path)
        if compressed:
            return compressed
        try:
            raw = path.read_bytes()
            mime = f"image/{'jpeg' if ext == 'jpg' else ext}"
            return f"data:{mime};base64,{base64.b64encode(raw).decode('ascii')}"
        except Exception:
            pass
    return ""


def _compress_logo_image(path):
    try:
        from PIL import Image
        import io
    except ImportError:
        return ""
    try:
        img = Image.open(path)
        if img.mode in ("RGBA", "P"):
            img = img.convert("RGBA")
        else:
            img = img.convert("RGB")
        w, h = img.size
        if h > 0 and w > 0:
            aspect = w / h
            if 0.7 <= aspect <= 1.4:
                crop_h = h // 2
                top = (h - crop_h) // 2
                img = img.crop((0, top, w, top + crop_h))
                w, h = img.size
        max_height = 120
        if h > max_height:
            ratio = max_height / h
            img = img.resize((int(w * ratio), max_height), Image.LANCZOS)
        buf = io.BytesIO()
        img.save(buf, format="PNG", optimize=True)
        raw = buf.getvalue()
        if len(raw) > 8192:
            buf2 = io.BytesIO()
            img = img.convert("RGB")
            img.save(buf2, format="JPEG", quality=70, optimize=True)
            raw = buf2.getvalue()
        if len(raw) > 12000:
            img2 = img.convert("RGB")
            max_h2 = 80
            w2, h2 = img2.size
            if h2 > max_h2:
                ratio2 = max_h2 / h2
                img2 = img2.resize((int(w2 * ratio2), max_h2), Image.LANCZOS)
            buf3 = io.BytesIO()
            img2.save(buf3, format="JPEG", quality=50, optimize=True)
            raw = buf3.getvalue()
            mime = "image/jpeg"
        else:
            mime = "image/png"
        return f"data:{mime};base64,{base64.b64encode(raw).decode('ascii')}"
    except Exception:
        return ""


def _normalize_stats_record(record):
    if not isinstance(record, dict):
        return {}
    out = {}
    for k, v in record.items():
        if k is None:
            continue
        key = str(k).strip()
        out[key] = v
        lk = key.lower()
        canon = _CANONICAL_STATS_KEYS.get(lk)
        if canon:
            if canon not in out or out.get(canon) in (None, ""):
                out[canon] = v
    return out


def _to_float(value):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    s = str(value).strip()
    if not s:
        return None
    try:
        return float(s.replace(",", ""))
    except Exception:
        return None


_DEFAULT_WELL_QC_THRESHOLDS = {
    "min_reads": 1000.0,
    "min_mapping_ratio": 0.30,
    "min_genes": 100.0,
    "min_umis": 100.0,
}


def _resolve_well_qc_thresholds(config, sample_type=""):
    """Return validated report-level thresholds for Active wells.

    A configured null/empty value disables that criterion. Invalid values use
    the default so a malformed optional report setting cannot break report
    generation.

    Manual and custom runs list their wells explicitly, so those wells are the
    user's own selection and are not re-filtered unless the user configured
    ``well_qc`` values. Plate/auto runs use the defaults because they contain
    empty wells that must be separated from real samples.
    """
    mode = str(
        sample_type
        or (config.get("sample") or {}).get("sample_type")
        or ""
    ).strip().lower()
    apply_defaults = mode in ("", "auto", "plate", "discover")

    configured = config.get("well_qc", {}) if isinstance(config, dict) else {}
    if not isinstance(configured, dict):
        configured = {}

    thresholds = {}
    for key, default in _DEFAULT_WELL_QC_THRESHOLDS.items():
        if key in configured:
            raw = configured[key]
        else:
            raw = default if apply_defaults else None
        if raw is None or str(raw).strip().lower() in {"", "none", "null", "~"}:
            thresholds[key] = None
            continue

        value = _to_float(raw)
        valid = value is not None and value >= 0
        if key == "min_mapping_ratio":
            valid = valid and value <= 1
        if not valid:
            logger.warning(
                "Invalid well_qc.%s=%r; using default %s",
                key,
                raw,
                default,
            )
            value = default
        thresholds[key] = value

    # Intron_Exon_genes and Intron_Exon_umis are derived from the same UMI
    # matrix. If the gene threshold is higher, it is the stricter criterion
    # and the UMI threshold becomes redundant; preserve the user's settings.
    if (
        thresholds["min_genes"] is not None
        and thresholds["min_umis"] is not None
        and thresholds["min_genes"] > thresholds["min_umis"]
    ):
        logger.warning(
            "well_qc.min_genes=%s exceeds min_umis=%s; the UMI criterion is redundant",
            thresholds["min_genes"],
            thresholds["min_umis"],
        )
    return thresholds


def _evaluate_well_qc(reads, mapping, genes, umis, thresholds):
    """Evaluate one well and return ``(passed, failed_criteria)``."""
    values = {
        "reads": (reads, "min_reads"),
        "mapping": (mapping, "min_mapping_ratio"),
        "genes": (genes, "min_genes"),
        "umis": (umis, "min_umis"),
    }
    failed = []
    for label, (value, key) in values.items():
        threshold = thresholds.get(key)
        if threshold is None:
            continue
        if value is None or value < threshold:
            failed.append(label)
    return not failed, failed


def _format_well_qc_threshold(key, value):
    if value is None:
        return "disabled"
    if key == "min_mapping_ratio":
        return f"{value * 100:.1f}%"
    if float(value).is_integer():
        return f"{int(value):,}"
    return f"{value:g}"


def _well_qc_help_text(thresholds):
    """Build the explanation shown by the Active-wells help button."""
    return (
        "A well is counted as Active only when it passes every enabled criterion.\n"
        f"Reads: >= {_format_well_qc_threshold('min_reads', thresholds.get('min_reads'))}\n"
        f"Mapping ratio: >= {_format_well_qc_threshold('min_mapping_ratio', thresholds.get('min_mapping_ratio'))}\n"
        f"Genes: >= {_format_well_qc_threshold('min_genes', thresholds.get('min_genes'))}\n"
        f"UMIs: >= {_format_well_qc_threshold('min_umis', thresholds.get('min_umis'))}\n"
        "PE data are counted as read pairs. A missing metric does not pass an enabled criterion."
    )


def _fmt_int(value):
    n = _to_float(value)
    if n is None:
        return ""
    return f"{int(round(n)):,}"


def _fmt_pct(value):
    n = _to_float(value)
    if n is None:
        return ""
    return f"{n * 100.0:.1f}%"


def _format_quality_value(value, value_type):
    if value is None:
        return "NA"
    if value_type == "pct":
        return _fmt_pct(value) or "NA"
    if value_type == "int":
        return _fmt_int(value) or "NA"
    return str(value)


def _max_end_from_base_definition(definition):
    if isinstance(definition, str):
        parts = definition.split(";")
    elif isinstance(definition, list):
        parts = definition
    else:
        parts = []
    max_end = None
    for part in parts:
        for _name, body in re.findall(r"(\w+)\(([^)]*)\)", str(part)):
            for rng in str(body).split(","):
                rng = rng.strip()
                if "-" not in rng:
                    continue
                try:
                    end = int(rng.split("-", 1)[1])
                except Exception:
                    continue
                max_end = end if max_end is None else max(max_end, end)
    return max_end


def _configured_read_length(config, file_key):
    seq_files = (config or {}).get("sequence_files", {}) or {}
    entry = seq_files.get(file_key, {}) if isinstance(seq_files, dict) else {}
    return _max_end_from_base_definition(entry.get("base_definition", []))


def _resolve_project_artifact(sample_outdir, project, relative_dir, suffix):
    """Resolve a project artifact exactly, with an unambiguous legacy fallback."""
    root = Path(sample_outdir)
    base = root / relative_dir if relative_dir else root
    project = str(project or "").strip()
    if project:
        exact = base / f"{project}{suffix}"
        if exact.is_file():
            return exact

        # Older layouts occasionally placed the same project-named artifact in
        # another subdirectory. Never fall back to a different project's file.
        candidates = sorted(path for path in root.rglob(f"{project}{suffix}") if path.is_file())
        candidates = list(dict.fromkeys(candidates))
        if len(candidates) == 1:
            logger.warning("Using legacy artifact location for project %s: %s", project, candidates[0])
            return candidates[0]
        if len(candidates) > 1:
            logger.warning(
                "Refusing ambiguous artifact selection for project %s and suffix %s: %s",
                project,
                suffix,
                ", ".join(str(path) for path in candidates),
            )
        return None

    candidates = sorted(path for path in base.glob(f"*{suffix}") if path.is_file()) if base.exists() else []
    if not candidates:
        candidates = sorted(path for path in root.rglob(f"*{suffix}") if path.is_file())
    candidates = list(dict.fromkeys(candidates))
    if len(candidates) == 1:
        logger.warning("Using legacy unscoped artifact for project %s: %s", project or "<unknown>", candidates[0])
        return candidates[0]
    if len(candidates) > 1:
        logger.warning(
            "Refusing ambiguous artifact selection for project %s and suffix %s: %s",
            project or "<unknown>",
            suffix,
            ", ".join(str(path) for path in candidates),
        )
    return None


def _load_q30_rows(sample_outdir, project=""):
    path = _resolve_project_artifact(sample_outdir, project, "stats", ".q30_stats.tsv")
    if path is None:
        return {}
    rows = {}
    with open(path, "r", encoding="utf-8", errors="ignore", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metric = str(row.get("metric") or "").strip()
            if not metric:
                continue
            rows[metric] = {
                "total_bases": _to_float(row.get("total_bases")),
                "q30_bases": _to_float(row.get("q30_bases")),
                "q30_rate": _to_float(row.get("q30_rate")),
            }
    return rows


def _sum_bcstats_reads(sample_outdir, project=""):
    path = _resolve_project_artifact(sample_outdir, project, "", ".BCstats.txt")
    if path is None:
        return None
    total = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            try:
                total += int(parts[1])
            except Exception:
                continue
    return total if total > 0 else None


def _load_read_stats_json(sample_outdir, project=""):
    path = _resolve_project_artifact(sample_outdir, project, "stats", ".read_stats.json")
    if path is None:
        return {}
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        data = json.load(handle)
    return data.get("read_stats", {}) or {}


def _read_count_from_read_stats_bucket(counts):
    if not isinstance(counts, dict):
        return 0
    umi_reads = _to_float(counts.get("UMI_Reads"))
    internal_reads = _to_float(counts.get("Internal_Reads"))
    if umi_reads is not None or internal_reads is not None:
        return int((umi_reads or 0) + (internal_reads or 0))

    total = 0
    for key, value in counts.items():
        if key in {"UMI_Reads", "Internal_Reads"}:
            continue
        numeric = _to_float(value)
        if numeric is not None:
            total += int(numeric)
    return total


def _barcode_read_counts_from_read_stats(read_stats):
    if not isinstance(read_stats, dict) or not read_stats:
        return None, None, None, None
    valid_reads = 0
    for barcode, counts in read_stats.items():
        if not isinstance(barcode, str) or barcode.startswith("__"):
            continue
        valid_reads += _read_count_from_read_stats_bucket(counts)

    unused_reads = int(_to_float((read_stats.get("__NO_CB__", {}) or {}).get("Unused BC")) or 0)
    total_reads = valid_reads + unused_reads
    valid_rate = (valid_reads / total_reads) if total_reads else None
    return total_reads, valid_reads, unused_reads, valid_rate


def _infer_transcriptome_label(reference_config):
    if not isinstance(reference_config, dict):
        return "Unknown"
    for key in ("transcriptome_name", "transcriptome", "reference_name"):
        value = str(reference_config.get(key) or "").strip()
        if value:
            return value

    star_index = str(reference_config.get("STAR_index") or "").strip()
    if not star_index:
        return "Unknown"
    path = Path(star_index.rstrip("/"))
    name = path.name
    if name.lower() in {"star", "star_index", "starindex", "index", "genome"} and path.parent.name:
        return path.parent.name
    return name or "Unknown"


def _median(values):
    vals = [float(v) for v in values if v is not None]
    if not vals:
        return None
    return statistics.median(vals)


def _load_expected_barcode_rows(sample_outdir):
    path = Path(config_dir(str(sample_outdir))) / "expect_id_barcode.tsv"
    if not path.exists():
        return []
    rows = []
    try:
        with open(path, "r", encoding="utf-8", errors="ignore", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if not row:
                    continue
                rows.append({
                    "wellID": str(row.get("wellID", "")).strip(),
                    "umi_barcodes": str(row.get("umi_barcodes", "")).strip(),
                    "internal_barcodes": str(row.get("internal_barcodes", "")).strip(),
                })
    except Exception:
        return []
    return rows


def _well_plate_id(well_id):
    m = re.match(r"^P(\d+)[A-P]\d+$", str(well_id or "").strip())
    return m.group(1) if m else ""


def _process_barcode_report_data(sample_outdir, combined_context, config):
    records = []
    try:
        records = json.loads(combined_context.get("rna_stats_table_data", "[]") or "[]")
    except Exception:
        records = []
    if not isinstance(records, list):
        records = []

    configured_sample_type = str((config.get("sample") or {}).get("sample_type") or "").strip().lower()
    sample_type = str(combined_context.get("sample_type") or configured_sample_type).strip().lower()
    barcode_source = str(config.get("barcode_source") or "").strip()
    expected_rows = _load_expected_barcode_rows(sample_outdir)
    # Normalise case so "p1a1" from expect_id_barcode.tsv still matches the
    # "P1A1" well IDs used in stats.tsv; a case mismatch would otherwise mark
    # every well as unexpected and silently zero the Active-well count.
    expected_by_well = {
        str(row["wellID"]).strip().upper(): row
        for row in expected_rows
        if row.get("wellID")
    }
    # Use the configured input mode for QC defaults. A discover run may select
    # the manual report template after barcode discovery, but it still needs
    # the plate/discover defaults for Active-well filtering.
    well_qc_mode = configured_sample_type or sample_type
    well_qc_thresholds = _resolve_well_qc_thresholds(config, well_qc_mode)
    active_help = _well_qc_help_text(well_qc_thresholds)

    total_reads = []
    umi_reads = []
    active_reads = []
    active_genes = []
    active_read_genes = []
    active_umis = []
    active_mapping_ratios = []
    # Fallback pools: metrics for every well regardless of QC outcome. Used only
    # when no well passes, so the report shows data instead of blank medians.
    all_genes = []
    all_read_genes = []
    all_umis = []
    all_mapping_ratios = []
    active = 0
    plate_ids = set()
    manual_rows = []
    plate_summary = {}
    well_qc_status = {}

    for raw in records:
        row = _normalize_stats_record(raw)
        well = str(row.get("wellID") or "").strip()
        if not well:
            continue
        well_key = well.upper()
        is_expected = not expected_rows or well_key in expected_by_well

        internal = _to_float(row.get("internal_reads")) or 0.0
        umi = _to_float(row.get("umi_reads")) or 0.0
        all_reads = _to_float(row.get("all_reads"))
        if all_reads is None:
            all_reads = internal + umi
        total_reads.append(all_reads)
        umi_reads.append(umi)

        gene_val = _to_float(row.get("Intron_Exon_genes"))
        if gene_val is None:
            gene_val = _to_float(row.get("Exon_genes"))
        read_gene_val = _to_float(row.get("Intron_Exon_read_genes"))
        if read_gene_val is None:
            read_gene_val = _to_float(row.get("Exon_read_genes"))
        umi_val = _to_float(row.get("Intron_Exon_umis"))
        if umi_val is None:
            umi_val = _to_float(row.get("Exon_umis"))

        mapping = _to_float(row.get("MappingRatio"))
        if mapping is None:
            exon_r = _to_float(row.get("Exon_reads")) or 0.0
            intron_r = _to_float(row.get("intron_reads")) or _to_float(row.get("Intron_reads")) or 0.0
            intergenic_r = _to_float(row.get("Intergenic_reads")) or 0.0
            ambiguity_r = _to_float(row.get("Ambiguity_reads")) or 0.0
            unmapped_r = _to_float(row.get("Unmapped_reads")) or 0.0
            other_unassigned_r = _to_float(row.get("Other_Unassigned_reads")) or 0.0
            denom = exon_r + intron_r + intergenic_r + ambiguity_r + other_unassigned_r + unmapped_r
            mapping = (
                (exon_r + intron_r + intergenic_r + ambiguity_r + other_unassigned_r) / denom
                if denom > 0
                else None
            )
        is_active, active_failures = _evaluate_well_qc(
            all_reads,
            mapping,
            gene_val,
            umi_val,
            well_qc_thresholds,
        )
        reported_active = is_active and is_expected
        if reported_active:
            active += 1
            active_reads.append(all_reads)
            if gene_val is not None:
                active_genes.append(gene_val)
            if read_gene_val is not None:
                active_read_genes.append(read_gene_val)
            if umi_val is not None:
                active_umis.append(umi_val)
            if mapping is not None:
                active_mapping_ratios.append(mapping)
        # The fallback is still restricted to expected wells. Unexpected rows
        # may be stale or diagnostic artifacts and must not affect report
        # medians even when no expected well passes QC.
        if is_expected:
            if gene_val is not None:
                all_genes.append(gene_val)
            if read_gene_val is not None:
                all_read_genes.append(read_gene_val)
            if umi_val is not None:
                all_umis.append(umi_val)
            if mapping is not None:
                all_mapping_ratios.append(mapping)
        if not is_expected and expected_rows:
            active_failures = list(active_failures) + ["unexpected"]
        well_qc_status[well] = {
            "active": reported_active,
            "expected": is_expected,
            "failed": active_failures,
        }

        plate = _well_plate_id(well)
        if plate:
            plate_ids.add(plate)
            bucket = plate_summary.setdefault(plate, {"plate": plate, "expected_wells": 0, "active_wells": 0, "reads": [], "genes": [], "read_genes": [], "umis": []})
            bucket["active_wells"] += 1 if reported_active else 0
            if reported_active:
                bucket["reads"].append(all_reads)
                bucket["genes"].append(gene_val)
                bucket["read_genes"].append(read_gene_val)
                bucket["umis"].append(umi_val)

        expected = expected_by_well.get(well_key, {})
        # Calculate the genic fraction used by the report. ExonIntronRatio is
        # retained separately as the legacy exon/intron quotient.
        umi_fraction = (umi / all_reads) if all_reads and all_reads > 0 else None
        exon_r = _to_float(row.get("Exon_reads")) or 0.0
        intron_r = _to_float(row.get("intron_reads")) or _to_float(row.get("Intron_reads")) or 0.0
        genic_ratio = _to_float(row.get("GenicRatio"))
        if genic_ratio is None:
            genic_ratio = ((exon_r + intron_r) / all_reads) if all_reads and all_reads > 0 else None
        exon_intron_ratio = _to_float(row.get("ExonIntronRatio"))
        if exon_intron_ratio is None:
            exon_intron_ratio = (exon_r / intron_r) if intron_r > 0 else None
        manual_rows.append({
            "wellID": well,
            "internal_barcodes": expected.get("internal_barcodes") or row.get("internal_barcodes") or "",
            "umi_barcodes": expected.get("umi_barcodes") or row.get("umi_barcodes") or "",
            "all_reads": int(all_reads),
            "umi_reads": int(umi),
            "internal_reads": int(internal),
            "genes": int(gene_val) if gene_val is not None else "",
            "read_genes": int(read_gene_val) if read_gene_val is not None else "",
            "umis": int(umi_val) if umi_val is not None else "",
            "mapping_ratio": mapping,
            "UMIfrac": umi_fraction,
            "umi_fraction": umi_fraction,
            "Exon_reads": exon_r,
            "intron_reads": intron_r,
            "Intron_reads": intron_r,
            "ExonIntronRatio": exon_intron_ratio,
            "exon_intron_ratio": exon_intron_ratio,
            "GenicRatio": genic_ratio,
            "genic_ratio": genic_ratio,
            "active": reported_active,
            "active_status": "Active" if reported_active else "Not active",
            "active_failures": ", ".join(active_failures),
        })

    if expected_rows:
        for row in expected_rows:
            plate = _well_plate_id(row.get("wellID"))
            if plate:
                plate_ids.add(plate)
                bucket = plate_summary.setdefault(plate, {"plate": plate, "expected_wells": 0, "active_wells": 0, "reads": [], "genes": [], "read_genes": [], "umis": []})
                bucket["expected_wells"] += 1
    else:
        for plate in plate_summary:
            plate_summary[plate]["expected_wells"] = len([r for r in records if _well_plate_id((r or {}).get("wellID")) == plate])

    total_read_sum = sum(total_reads)
    umi_read_sum = sum(umi_reads)
    expected_count = len(expected_rows) if expected_rows else len(records)
    summary = {
        "sample_type": sample_type or "unknown",
        "configured_sample_type": configured_sample_type or "unknown",
        "sample_type_label": {"auto": "Auto 384-well plate", "manual": "Manual samples", "custom": "Custom barcode set"}.get(sample_type, sample_type or "Unknown"),
        "barcode_source": barcode_source,
        "expected_wells": expected_count,
        "active_wells": active,
        "active_fraction": (active / expected_count) if expected_count else None,
        "active_well_thresholds": well_qc_thresholds,
        "plate_count": len(plate_ids),
        "total_reads": total_read_sum,
        "median_reads": _median(active_reads),
        "median_genes": _median(active_genes),
        "median_read_genes": _median(active_read_genes),
        "median_umis": _median(active_umis),
        "median_mapping_ratio": _median(active_mapping_ratios),
        "umi_fraction": (umi_read_sum / total_read_sum) if total_read_sum > 0 else None,
    }

    # No well passed QC: median cards based on Active wells would all be blank,
    # which reads like "the run produced nothing" even though data exist. Fall
    # back to every well and flag it so the report can say so explicitly.
    qc_all_failed = bool(expected_count) and active == 0 and bool(records)
    summary["qc_all_failed"] = qc_all_failed
    # Mirror the flag into well_qc_status so the JS summary tables can use
    # all-wells medians too instead of rendering blanks for the Active subset.
    well_qc_status["__all_failed__"] = qc_all_failed
    if qc_all_failed:
        summary["median_reads"] = _median(total_reads)
        summary["median_genes"] = _median(all_genes)
        summary["median_read_genes"] = _median(all_read_genes)
        summary["median_umis"] = _median(all_umis)
        summary["median_mapping_ratio"] = _median(all_mapping_ratios)
        logger.warning(
            "No well passed the Active-well QC thresholds (%s); median metrics "
            "fall back to all %d well(s).",
            well_qc_thresholds,
            expected_count,
        )

    # When nothing passed QC the medians come from every well, so the card help
    # must say that instead of claiming they are Active-well medians.
    median_scope = (
        "all wells (no well passed QC)" if qc_all_failed else "Active wells"
    )
    qc_warning = None
    if qc_all_failed:
        qc_warning = (
            f"No well passed the Active-well QC thresholds, so median metrics "
            f"below are computed across all {expected_count} well(s) instead of "
            f"Active wells. Review the well_qc thresholds in the configuration."
        )

    cards = [
        {
            "label": "Barcode mode",
            "value": summary["sample_type_label"],
        },
        {
            "label": "Expected wells",
            "value": _fmt_int(summary["expected_wells"]),
            "help_title": "Expected wells",
            "help": "Number of wells listed in expect_id_barcode.tsv. Wells with zero reads remain in this denominator.",
        },
        {
            "label": "Active wells",
            "value": _fmt_int(summary["active_wells"]),
            "help_title": "Active wells",
            "help": active_help,
            "warning": qc_warning,
        },
        {
            "label": "Median reads",
            "value": _fmt_int(summary["median_reads"]),
            "help_title": "Median reads",
            "help": f"Median all_reads across {median_scope}. In PE mode, one read pair is one unit.",
        },
        {
            "label": "Median genes (UMIs)",
            "value": _fmt_int(summary["median_genes"]),
            "help_title": "Median genes (UMIs)",
            "help": f"Median detected genes across {median_scope}, using the exon plus intron UMI matrix when available.",
        },
        {
            "label": "Median UMIs",
            "value": _fmt_int(summary["median_umis"]),
            "help_title": "Median UMIs",
            "help": f"Median unique molecular identifiers across {median_scope}, using the exon plus intron UMI matrix when available.",
        },
        {
            "label": "Median mapping",
            "value": _fmt_pct(summary["median_mapping_ratio"]),
            "help_title": "Median mapping",
            "help": f"Median MappingRatio across {median_scope}. MappingRatio is the fraction of counted reads that are not unmapped.",
        },
    ]

    plate_rows = []
    for plate in sorted(plate_summary, key=lambda x: int(x) if str(x).isdigit() else str(x)):
        item = plate_summary[plate]
        expected = item["expected_wells"] or len(item["reads"])
        plate_rows.append({
            "plate": plate,
            "expected_wells": expected,
            "active_wells": item["active_wells"],
            "active_fraction": (item["active_wells"] / expected) if expected else None,
            "median_reads": _median(item["reads"]),
            "median_genes": _median(item["genes"]),
            "median_read_genes": _median(item.get("read_genes", [])),
            "median_umis": _median(item["umis"]),
        })

    manual_rows.sort(key=lambda x: str(x.get("wellID", "")))
    combined_context["barcode_report_summary_data"] = json.dumps(summary)
    combined_context["barcode_mode_cards_data"] = json.dumps(cards)
    combined_context["manual_barcode_table_data"] = json.dumps(manual_rows)
    combined_context["auto_plate_summary_data"] = json.dumps(plate_rows)
    combined_context["well_qc_status_data"] = json.dumps(well_qc_status)


def _read_discovered_sample_type(sample_outdir, project=""):
    path = _resolve_project_artifact(sample_outdir, project, "barcodes", ".barcode_discovery.tsv")
    candidates = [path] if path is not None else []

    for path in candidates:
        in_summary = False
        try:
            with open(path, "r", encoding="utf-8") as handle:
                reader = None
                for line in handle:
                    line = line.strip()
                    if not line:
                        continue
                    if line == "[summary]":
                        in_summary = True
                        reader = None
                        continue
                    if line.startswith("[") and line != "[summary]":
                        in_summary = False
                        continue
                    if not in_summary:
                        continue
                    if reader is None:
                        reader = line.split("\t")
                        continue
                    values = line.split("\t")
                    row = dict(zip(reader, values))
                    sample_type = str(row.get("candidate_type") or "").strip().lower()
                    if sample_type in ("auto", "manual"):
                        return sample_type
        except Exception:
            continue
    return ""


def _select_report_template(sample_type, sample_outdir, template_dir, config=None):
    sample_type = str(sample_type or "").strip().lower()
    template_auto = template_dir / "template_auto.html"
    template_manual = template_dir / "template_manual.html"

    sample_cfg = (config or {}).get("sample") or {}
    discovered_type = str(sample_cfg.get("discovered_sample_type") or "").strip().lower()
    if discovered_type in ("auto", "manual"):
        report_mode = discovered_type
    elif sample_type in ("auto", "manual"):
        report_mode = sample_type
    elif sample_type == "discover":
        report_mode = _read_discovered_sample_type(sample_outdir, (config or {}).get("project")) or "auto"
    else:
        report_mode = "manual"

    template_path = template_auto if report_mode == "auto" else template_manual
    return template_path, report_mode


def _count_lines_in_tsv_gz(path):
    try:
        p = Path(path)
        if not p.exists():
            return None
        with gzip.open(p, "rt", encoding="utf-8", errors="ignore") as f:
            n = 0
            for _ in f:
                n += 1
            return n
    except Exception:
        return None


def _infer_total_genes_from_expression(sample_outdir, project):
    if not project:
        return None
    base = Path(expression_dir(sample_outdir))
    candidates = [
        base / f"{project}.inex.umi" / "genes.tsv.gz",
        base / f"{project}.inex.umi" / "features.tsv.gz",
        base / f"{project}.exon.umi" / "genes.tsv.gz",
        base / f"{project}.exon.umi" / "features.tsv.gz",
    ]
    for p in candidates:
        n = _count_lines_in_tsv_gz(p)
        if n is not None and n > 0:
            return n
    return None


def _process_rna_cluster_assignment(_sample_outdir, combined_context):
    combined_context["rna_umap_x"] = "[]"
    combined_context["rna_umap_y"] = "[]"
    combined_context["rna_leiden"] = "[]"
    combined_context["rna_numi"] = "[]"
    combined_context["rna_clusters"] = "[]"
    combined_context["cell_annotation_available"] = "false"
    return


def _process_rna_bead_count_data(_sample_outdir, combined_context):
    combined_context["bead_count_x"] = "[]"
    combined_context["bead_count_y"] = "[]"
    return


def _process_rna_rankplot_data(_sample_outdir, combined_context):
    combined_context["rna_rankplot_data"] = "[]"
    return


def infer_enabled_sections(config):
    """
    Infer which sections (rna, vdj, etc.) are enabled based on config
    or directory structure.
    """
    enabled = []
    # For now, we mainly support RNA in this pipeline
    if config.get('sequence_files'):
        enabled.append('rna')
    
    # Check if VDJ/ATAC outputs exist (future proofing)
    if 'out_dir' in config:
        out_dir = Path(config['out_dir'])
        if (out_dir / 'VDJ-T_ANALYSIS_WORKFLOW_PROCESSING').exists():
            enabled.append('vdj-t')
        if (out_dir / 'VDJ-B_ANALYSIS_WORKFLOW_PROCESSING').exists():
            enabled.append('vdj-b')
        if (out_dir / 'ATAC_ANALYSIS_WORKFLOW_PROCESSING').exists():
            enabled.append('atac')
            
    return enabled

def create_report_directories(sample_outdir, _config=None):
    """
    Create necessary REPORT directory structure for multi-omics analysis
    based on the configuration
    """
    sample_path = Path(sample_outdir)

    # For MfsFlow, we might not need strict directory structure for REPORT 
    # if we are just generating one HTML, but we'll keep it for consistency.
    # Locate pipeline output tables from the current directory layout.
    
    # Create final output directory 'outs'
    # In the pipeline, out_dir is XPRESS_PROCESSING. 
    # sample_outdir passed here is likely XPRESS_PROCESSING.
    # The final report should go to 'outs' in the project root.
    
    outs_dir = Path(outputs_dir(str(sample_path)))
    outs_dir.mkdir(parents=True, exist_ok=True)
    
    return outs_dir


def _copy_or_link(src, dst):
    """Hard-link a file when possible, falling back to a metadata-preserving copy."""
    try:
        os.link(str(src), str(dst))
    except OSError:
        shutil.copy2(str(src), str(dst))
    return str(dst)


def _export_file_if_exists(src, dst):
    src = Path(src)
    dst = Path(dst)
    if not src.exists() or not src.is_file():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    temporary = dst.with_name(f".{dst.name}.{os.getpid()}.tmp")
    try:
        if temporary.exists():
            temporary.unlink()
        _copy_or_link(src, temporary)
        os.replace(temporary, dst)
    finally:
        if temporary.exists():
            temporary.unlink()
    return True


def _export_tree_if_exists(src, dst):
    src = Path(src)
    dst = Path(dst)
    if not src.exists() or not src.is_dir():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    temporary = dst.with_name(f".{dst.name}.{os.getpid()}.tmp")
    if temporary.exists():
        shutil.rmtree(temporary)
    try:
        shutil.copytree(src, temporary, copy_function=_copy_or_link)
        if dst.exists():
            shutil.rmtree(dst)
        os.replace(temporary, dst)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return True


def export_deliverables_to_outs(sample_outdir, outs_dir, project):
    sample_outdir = Path(sample_outdir)
    outs_dir = Path(outs_dir)
    project = str(project or "").strip()
    exported_paths = []

    expression_out = Path(expression_dir(str(sample_outdir)))
    if expression_out.exists():
        for h5ad in sorted(expression_out.glob(f"{project}.h5ad")):
            dst = outs_dir / "expression" / h5ad.name
            if _export_file_if_exists(h5ad, dst):
                exported_paths.append(str(dst))
        for mex_dir in sorted(expression_out.glob(f"{project}.*")):
            if not mex_dir.is_dir():
                continue
            if (
                (mex_dir / "matrix.mtx.gz").exists()
                and (mex_dir / "features.tsv.gz").exists()
                and (mex_dir / "barcodes.tsv.gz").exists()
            ):
                dst = outs_dir / "expression" / mex_dir.name
                if _export_tree_if_exists(mex_dir, dst):
                    exported_paths.append(str(dst))

    stats_out = Path(stats_dir(str(sample_outdir)))
    if stats_out.exists():
        for src in sorted(stats_out.glob(f"{project}.*")):
            if not src.is_file():
                continue
            dst = outs_dir / "stats" / src.name
            if _export_file_if_exists(src, dst):
                exported_paths.append(str(dst))

    bam_patterns = [
        f"{project}.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam*",
        f"{project}.filtered.Aligned.GeneTagged.UBcorrected.bam*",
    ]
    seen = set()
    for pattern in bam_patterns:
        for src in sorted(sample_outdir.glob(pattern)):
            if src in seen or not src.is_file():
                continue
            seen.add(src)
            dst = outs_dir / "bam" / src.name
            if _export_file_if_exists(src, dst):
                exported_paths.append(str(dst))

    cfg_out = Path(config_dir(str(sample_outdir)))
    if cfg_out.exists():
        for pattern in ("run_config.yaml", "expect_id_barcode.tsv"):
            for src in sorted(cfg_out.glob(pattern)):
                dst = outs_dir / "config" / src.name
                if _export_file_if_exists(src, dst):
                    exported_paths.append(str(dst))

    exported = {os.path.abspath(path) for path in exported_paths}
    stale_candidates = []
    for category in ("expression", "stats", "bam"):
        category_dir = outs_dir / category
        if category_dir.exists():
            stale_candidates.extend(category_dir.glob(f"{project}.*"))
    for stale in stale_candidates:
        if os.path.abspath(str(stale)) in exported:
            continue
        if stale.is_dir():
            shutil.rmtree(stale)
        elif stale.is_file():
            stale.unlink()

    logger.info(f"Exported {len(exported_paths)} deliverable file(s)/dir(s) to: {outs_dir}")
    return exported_paths


def _copy_snapshot_file(src, dst):
    src = Path(src)
    dst = Path(dst)
    if not src.is_file():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    temporary = dst.with_name(f".{dst.name}.{os.getpid()}.tmp")
    try:
        shutil.copy2(src, temporary)
        os.replace(temporary, dst)
    finally:
        if temporary.exists():
            temporary.unlink()
    return True


def export_diagnostics_to_outs(sample_outdir, outs_dir, project, include_config=False):
    """Copy a stable diagnostic snapshot without modifying resume inputs."""
    sample_outdir = Path(sample_outdir)
    outs_dir = Path(outs_dir)
    diagnostic_root = outs_dir / "diagnostics"
    copied = []

    fixed_files = [
        (sample_outdir / "logs" / "pipeline.log", diagnostic_root / "logs" / "pipeline.log"),
        (sample_outdir / "logs" / "pipeline_timing.tsv", diagnostic_root / "logs" / "pipeline_timing.tsv"),
        (sample_outdir / f"{project}.BCstats.txt", diagnostic_root / "barcodes" / f"{project}.BCstats.txt"),
    ]
    for src, dst in fixed_files:
        if _copy_snapshot_file(src, dst):
            copied.append(str(dst))

    stage_dir = sample_outdir / "logs" / "stages"
    if stage_dir.exists():
        for src in sorted(stage_dir.iterdir()):
            dst = diagnostic_root / "stages" / src.name
            if _copy_snapshot_file(src, dst):
                copied.append(str(dst))

    barcode_output = sample_outdir / "barcodes"
    if barcode_output.exists():
        for src in sorted(barcode_output.glob(f"{project}.*")):
            dst = diagnostic_root / "barcodes" / src.name
            if _copy_snapshot_file(src, dst):
                copied.append(str(dst))

    if include_config:
        for name in ("run_config.yaml", "expect_id_barcode.tsv"):
            src = sample_outdir / "config" / name
            dst = diagnostic_root / "config" / name
            if _copy_snapshot_file(src, dst):
                copied.append(str(dst))
    logger.info("Exported %d diagnostic file(s) to: %s", len(copied), diagnostic_root)
    return copied


def _file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_run_manifest(outs_dir, project, config, status, report_path):
    """Write a machine-readable inventory and provenance record for exported outputs."""
    outs_dir = Path(outs_dir)
    manifest_name = "run_manifest.json" if status == "success" else "partial_run_manifest.json"
    manifest_path = outs_dir / manifest_name
    artifacts = []
    checksum_limit = 64 * 1024 * 1024
    for path in sorted(item for item in outs_dir.rglob("*") if item.is_file()):
        if path == manifest_path or path.name.endswith(".tmp"):
            continue
        size = path.stat().st_size
        artifacts.append({
            "path": str(path.relative_to(outs_dir)),
            "size_bytes": size,
            "sha256": _file_sha256(path) if size <= checksum_limit else None,
            "checksum_note": None if size <= checksum_limit else "omitted for files larger than 64 MiB",
        })
    payload = {
        "schema_version": 1,
        "project": project,
        "status": status,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "mfsflow_version": __version__,
        "report": Path(report_path).name,
        "analysis_error": config.get("_analysis_error") if status != "success" else None,
        "tool_versions": config.get("tool_versions", {}),
        "reference": config.get("reference", {}),
        "artifacts": artifacts,
    }
    temporary = manifest_path.with_name(f".{manifest_path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
        os.replace(temporary, manifest_path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return manifest_path

def calculate_summary_metrics(sample_outdir, project=""):
    """
    Calculate summary metrics from stats.tsv
    """
    stats_path = _resolve_project_artifact(sample_outdir, project, "stats", ".stats.tsv")
        
    metrics = {
        "rna_estm_Num_cell": "",
        "rna_median_umis_percell": "",
        "rna_mean_umis_percell": "",
        "rna_median_genes_percell": "",
        "rna_mean_genes_percell": "",
        "rna_mean_reads_percell": "",
        "rna_total_gene": "",
        "rna_fraction_transcriptome_in_cells": "",
        "rna_species": "Unknown",
    }
    
    if stats_path is None:
        logger.warning("No stats.tsv found for summary metrics.")
        return metrics
        
    def to_float(v):
        if v is None:
            return None
        s = str(v).strip()
        if not s:
            return None
        try:
            return float(s.replace(",", ""))
        except Exception:
            return None

    def fmt_int(n):
        try:
            return f"{int(n):,}"
        except Exception:
            return ""

    try:
        with open(stats_path, "r", encoding="utf-8", errors="ignore", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = [_normalize_stats_record(r) for r in reader if r]
        if not rows:
            return metrics

        metrics["rna_estm_Num_cell"] = fmt_int(len(rows))

        umis = []
        genes = []
        reads = []
        exon_reads = []
        intron_reads = []
        intergenic_reads = []
        ambiguity_reads = []
        unmapped_reads = []

        for r in rows:
            u = to_float(r.get("Intron_Exon_umis")) if r.get("Intron_Exon_umis") is not None else None
            if u is None:
                u = to_float(r.get("Exon_umis")) if r.get("Exon_umis") is not None else None
            if u is not None:
                umis.append(u)

            g = to_float(r.get("Intron_Exon_genes")) if r.get("Intron_Exon_genes") is not None else None
            if g is None:
                g = to_float(r.get("Exon_genes")) if r.get("Exon_genes") is not None else None
            if g is not None:
                genes.append(g)

            ru = to_float(r.get("umi_reads")) if r.get("umi_reads") is not None else 0.0
            ri = to_float(r.get("internal_reads")) if r.get("internal_reads") is not None else 0.0
            if ru is not None or ri is not None:
                reads.append(float(ru or 0.0) + float(ri or 0.0))

            v = to_float(r.get("Exon_reads"))
            if v is not None:
                exon_reads.append(v)
            v = to_float(r.get("intron_reads"))
            if v is not None:
                intron_reads.append(v)
            v = to_float(r.get("Intergenic_reads"))
            if v is not None:
                intergenic_reads.append(v)
            v = to_float(r.get("Ambiguity_reads"))
            if v is not None:
                ambiguity_reads.append(v)
            v = to_float(r.get("Unmapped_reads"))
            if v is not None:
                unmapped_reads.append(v)

        if umis:
            metrics["rna_median_umis_percell"] = fmt_int(statistics.median(umis))
            metrics["rna_mean_umis_percell"] = fmt_int(sum(umis) / len(umis))
        if genes:
            metrics["rna_median_genes_percell"] = fmt_int(statistics.median(genes))
            metrics["rna_mean_genes_percell"] = fmt_int(sum(genes) / len(genes))
        if reads:
            metrics["rna_mean_reads_percell"] = fmt_int(sum(reads) / len(reads))

        if exon_reads or intron_reads or intergenic_reads or ambiguity_reads or unmapped_reads:
            exon_sum = sum(exon_reads) if exon_reads else 0.0
            intron_sum = sum(intron_reads) if intron_reads else 0.0
            intergenic_sum = sum(intergenic_reads) if intergenic_reads else 0.0
            ambig_sum = sum(ambiguity_reads) if ambiguity_reads else 0.0
            unmapped_sum = sum(unmapped_reads) if unmapped_reads else 0.0
            denom = exon_sum + intron_sum + intergenic_sum + ambig_sum + unmapped_sum
            transcriptome = exon_sum + intron_sum
            if denom > 0:
                metrics["rna_fraction_transcriptome_in_cells"] = f"{(transcriptome / denom) * 100.0:.2f}%"

        total_genes = _infer_total_genes_from_expression(sample_outdir, str(project))
        if total_genes is not None:
            metrics["rna_total_gene"] = fmt_int(total_genes)
    except Exception as e:
        logger.error(f"Error calculating summary metrics: {e}")
        
    return metrics

def get_omics_data(outdir, _sample, config):
    omics_data = {}
    outdir = Path(outdir)
    # In run_analysis_pipeline, outdir is passed as the XPRESS_PROCESSING directory.
    
    enabled = set(infer_enabled_sections(config))

    # RNA
    if 'rna' in enabled:
        try:
            # We calculate stats manually now
            project = config.get("project") or _sample or ""
            stat = calculate_summary_metrics(outdir, project=str(project))
            
            # Metadata
            rna_section = config.get('reference', {})
            stat['rna_species'] = _infer_transcriptome_label(rna_section)
            
            omics_data['rna'] = {'stat': stat, 'plot_dict': {}, 'table': None}
            
        except Exception as e:
            logger.warning(f"Failed to load RNA report data: {e}")

    return omics_data

# Helper functions for processing omics data
def _process_rna_data(omics_data, sample_outdir, combined_context):
    """Process RNA data and add to combined context"""
    if 'rna' not in omics_data:
        return
        
    # Add basic RNA statistics
    for key, value in omics_data['rna']['stat'].items():
        combined_context[key] = value
        
    # Process RNA-specific data
    _process_rna_stats_table_data(sample_outdir, combined_context)
    _process_barcode_report_data(sample_outdir, combined_context, combined_context.get("_run_config", {}))
    _process_rna_saturation_data(sample_outdir, combined_context)
    _process_rna_gene_body_coverage_data(sample_outdir, combined_context)
    _process_sequencing_quality_data(sample_outdir, combined_context)
    _process_rna_read_distribution_data(sample_outdir, combined_context)
    _process_rna_bead_count_data(sample_outdir, combined_context)
    _process_rna_cluster_assignment(sample_outdir, combined_context)
    _process_rna_rankplot_data(sample_outdir, combined_context)


def _process_sequencing_quality_data(sample_outdir, combined_context):
    combined_context["sequencing_quality_summary_data"] = "[]"
    try:
        config = combined_context.get("_run_config", {}) or {}
        project = config.get("project", "")
        q30_rows = _load_q30_rows(sample_outdir, project)
        read_stats = _load_read_stats_json(sample_outdir, project)

        total_reads, valid_bc_reads, unused_bc_reads, valid_bc_rate = _barcode_read_counts_from_read_stats(read_stats)
        if total_reads is None:
            r1_len = _configured_read_length(config, "file1")
            r1_total_bases = (q30_rows.get("R1") or {}).get("total_bases")
            total_reads = int(round(r1_total_bases / r1_len)) if r1_len and r1_total_bases else None
            valid_bc_reads = _sum_bcstats_reads(sample_outdir, project)
            # read_stats.json is missing, so Unused cannot be computed with a
            # matching metric semantics: raw-input minus BCstats would also count reads that
            # fqfilter dropped for quality, inflating the Unused number.
            unused_bc_reads = None
            valid_bc_rate = (valid_bc_reads / total_reads) if total_reads and valid_bc_reads is not None else None

        metric_map = [
            ("Total sequencing reads", total_reads, "int"),
            ("Valid barcode reads", valid_bc_reads, "int"),
            ("Unused barcode reads", unused_bc_reads, "int"),
            ("Valid barcode rate", valid_bc_rate, "pct"),
            ("Cell barcode Q30", (q30_rows.get("BC") or {}).get("q30_rate"), "pct"),
            ("UMI Q30", (q30_rows.get("UMI") or {}).get("q30_rate"), "pct"),
            ("Read1 cDNA Q30", (q30_rows.get("R1_cDNA") or {}).get("q30_rate"), "pct"),
            ("Read2 cDNA Q30", (q30_rows.get("R2_cDNA") or {}).get("q30_rate"), "pct"),
        ]
        payload = [
            {
                "metric": label,
                "value": _format_quality_value(value, value_type),
                "raw": value,
                "type": value_type,
            }
            for label, value, value_type in metric_map
        ]
        combined_context["sequencing_quality_summary_data"] = json.dumps(payload)
    except Exception as e:
        logger.warning(f"Failed to prepare sequencing quality summary data: {e}")

def _process_rna_stats_table_data(sample_outdir, combined_context):
    combined_context['rna_stats_table_data'] = '[]'
    combined_context['rna_stats_table_available'] = False
    try:
        project = (combined_context.get("_run_config", {}) or {}).get("project", "")
        stats_path = _resolve_project_artifact(sample_outdir, project, "stats", ".stats.tsv")
        if stats_path is None:
            return

        with open(stats_path, "r", encoding="utf-8", errors="ignore", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            records = [_normalize_stats_record(r) for r in reader if r]
        combined_context["rna_stats_table_data"] = json.dumps(records)
        combined_context["rna_stats_table_available"] = bool(len(records) > 0)
    except Exception as e:
        logger.warning(f"Failed to prepare RNA stats table data: {e}")

def _process_rna_gene_body_coverage_data(sample_outdir, combined_context):
    combined_context['rna_gene_body_percentile'] = '[]'
    combined_context['rna_gene_body_umi_maxnorm'] = '[]'
    combined_context['rna_gene_body_internal_maxnorm'] = '[]'
    combined_context['rna_gene_body_all_maxnorm'] = '[]'
    combined_context['rna_gene_body_coverage_available'] = False
    try:
        project = (combined_context.get("_run_config", {}) or {}).get("project", "")
        gb_path = _resolve_project_artifact(sample_outdir, project, "stats", ".geneBodyCoverage.txt")
        if gb_path is None:
            return

        x_vals = []
        umi_vals = []
        internal_vals = []
        all_vals = []
        with open(gb_path, "r", encoding="utf-8", errors="ignore", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                if not r:
                    continue
                try:
                    x = float(str(r.get("Percentile", "")).strip() or 0)
                    umi = float(str(r.get("UMI_MaxNorm", "")).strip() or 0)
                    internal = float(str(r.get("Internal_MaxNorm", "")).strip() or 0)
                    all_norm = float(str(r.get("All_MaxNorm", "")).strip() or 0)
                except Exception:
                    continue
                x_vals.append(x)
                umi_vals.append(umi)
                internal_vals.append(internal)
                all_vals.append(all_norm)

        combined_context["rna_gene_body_percentile"] = json.dumps(x_vals)
        combined_context["rna_gene_body_umi_maxnorm"] = json.dumps(umi_vals)
        combined_context["rna_gene_body_internal_maxnorm"] = json.dumps(internal_vals)
        combined_context["rna_gene_body_all_maxnorm"] = json.dumps(all_vals)
        combined_context["rna_gene_body_coverage_available"] = bool(len(x_vals) > 0 and len(x_vals) == len(umi_vals) == len(internal_vals) == len(all_vals))
    except Exception as e:
        logger.warning(f"Failed to prepare RNA gene body coverage data: {e}")

def _process_rna_read_distribution_data(sample_outdir, combined_context):
    combined_context['rna_read_distribution_bar_data'] = '{}'
    combined_context['rna_read_distribution_box_data'] = '{}'
    try:
        read_stats = {}
        project = (combined_context.get("_run_config", {}) or {}).get("project", "")
        stats_path = _resolve_project_artifact(sample_outdir, project, "stats", ".read_stats.json")
        if stats_path is not None:
            with open(stats_path, 'r', encoding='utf-8') as f:
                stats_json = json.load(f)
            read_stats = stats_json.get('read_stats', {}) or {}

        categories_order = ["Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped"]
        feat_colors = {
            "Exon": "#1A5084",
            "Intron": "#118730",
            "Intergenic": "#FFD700",
            "Ambiguity": "#FFA54F",
            "Unmapped": "#545454",
            "Unused BC": "#BABABA",
        }

        total_fracs_pct = {c: 0.0 for c in categories_order}
        unused_frac_pct = 0.0
        
        per_cell_pct = {c: [] for c in categories_order}
        
        # Calculate totals from read_stats.json
        if isinstance(read_stats, dict) and read_stats:
            total_barcodes = [k for k in read_stats.keys() if isinstance(k, str) and not k.startswith('__')]
            unused_total = 0
            try:
                unused_total = int((read_stats.get("__NO_CB__", {}) or {}).get("Unused BC", 0))
            except Exception:
                unused_total = 0

            totals = {c: 0 for c in categories_order}
            
            # Calculate per-cell percentages and total sums
            for bc in total_barcodes:
                counts = read_stats.get(bc, {}) or {}
                cell_total = 0
                for c in categories_order:
                    val = int(counts.get(c, 0) or 0)
                    totals[c] += val
                    cell_total += val
                
                if cell_total > 0:
                    for c in categories_order:
                        val = int(counts.get(c, 0) or 0)
                        per_cell_pct[c].append((val / cell_total) * 100.0)

            total_sum = sum(totals.values()) + unused_total
            total_fracs_pct = {c: (totals[c] / total_sum * 100.0) if total_sum > 0 else 0.0 for c in categories_order}
            unused_frac_pct = (unused_total / total_sum * 100.0) if total_sum > 0 else 0.0

        bar_payload = {
            "order": categories_order + ["Unused BC"],
            "values": {**total_fracs_pct, "Unused BC": unused_frac_pct},
            "colors": feat_colors,
        }
        box_payload = {
            "order": categories_order,
            "values": per_cell_pct,
            "colors": feat_colors,
        }

        combined_context['rna_read_distribution_bar_data'] = json.dumps(bar_payload)
        combined_context['rna_read_distribution_box_data'] = json.dumps(box_payload)
    except Exception as e:
        logger.warning(f"Failed to prepare RNA read distribution data: {e}")

def _process_rna_saturation_data(sample_outdir, combined_context):
    """Prepare RNA saturation arrays for plotting"""
    try:
        combined_context['rna_saturation_fraction'] = '[]'
        combined_context['rna_saturation_lib_pct'] = '[]'
        combined_context['rna_saturation_gene_pct'] = '[]'
        combined_context['rna_saturation_median_genes_umi'] = '[]'
        combined_context['rna_saturation_median_genes_read'] = '[]'
        combined_context['rna_saturation_available'] = False

        project = (combined_context.get("_run_config", {}) or {}).get("project", "")
        sat_path = _resolve_project_artifact(sample_outdir, project, "stats", ".saturation.tsv")
        if sat_path is not None:
            frac_vals = []
            sat_lib_pct = []
            sat_gene_pct = []
            med_umi_vals = []
            med_read_vals = []
            with open(sat_path, "r", encoding="utf-8", errors="ignore", newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for r in reader:
                    if not r:
                        continue
                    try:
                        frac = float(str(r.get("Fraction", "")).strip() or 0)
                        lib = float(str(r.get("Seq_Saturation_Library", "")).strip() or 0) * 100.0
                        gene = float(str(r.get("Seq_Saturation_Gene", "")).strip() or 0) * 100.0
                        med_umi = float(str(r.get("Median_Genes_UMI", "")).strip() or 0)
                        med_read = float(str(r.get("Median_Genes_Read", "")).strip() or 0)
                    except Exception:
                        continue
                    frac_vals.append(frac)
                    sat_lib_pct.append(lib)
                    sat_gene_pct.append(gene)
                    med_umi_vals.append(med_umi)
                    med_read_vals.append(med_read)

            if frac_vals and len(frac_vals) == len(sat_lib_pct) == len(sat_gene_pct):
                combined_context["rna_saturation_fraction"] = json.dumps(frac_vals)
                combined_context["rna_saturation_lib_pct"] = json.dumps(sat_lib_pct)
                combined_context["rna_saturation_gene_pct"] = json.dumps(sat_gene_pct)
                combined_context["rna_saturation_median_genes_umi"] = json.dumps(med_umi_vals)
                combined_context["rna_saturation_median_genes_read"] = json.dumps(med_read_vals)
                combined_context["rna_saturation_available"] = True

        combined_context['rna_saturation_sampling_fraction'] = combined_context.get('rna_saturation_fraction', '[]')
        combined_context['rna_saturation_seq_pct'] = combined_context.get('rna_saturation_lib_pct', '[]')
        combined_context['rna_saturation_median_genes'] = combined_context.get('rna_saturation_median_genes_umi', '[]')
        try:
            lib_pct = json.loads(combined_context.get("rna_saturation_lib_pct", "[]") or "[]")
            if isinstance(lib_pct, list) and lib_pct:
                combined_context["rna_saturation"] = f"{float(lib_pct[-1]):.1f}%"
        except Exception:
            pass
    except Exception as e:
        logger.warning(f"Failed to prepare RNA saturation arrays: {e}")


def generate_multi_report(name, outdir, config):
    """
    Generates the HTML analysis report.
    """
    logger.info(f"Generating HTML report for sample {name} in {outdir}")
    
    sample_outdir = Path(outdir)
    
    # Create necessary REPORT directories and outs directory based on config
    outs_dir = create_report_directories(sample_outdir, config)
    
    omics_data = get_omics_data(outdir, name, config)
    
    sample_type = str((config.get("sample") or {}).get("sample_type") or "").strip().lower()

    # Locate template. Auto and manual reports use separate templates because
    # 384-well plate QC and low-sample manual summaries have different layouts.
    template_dir = Path(__file__).parent / 'report_assets'
    template_path, report_mode = _select_report_template(sample_type, sample_outdir, template_dir, config)
    
    if not template_path.exists():
        raise FileNotFoundError(f"Report template not found: {template_path}")
    logger.info(f"Using report template: {template_path.name} (mode={report_mode}, configured={sample_type or 'unknown'})")

    with open(template_path, 'r', encoding='utf-8') as f:
        template_str = f.read()
    
    template = Template(template_str)
    
    combined_context = {}
    
    # Process each omics data type
    # Add sample name and version
    combined_context['samplename'] = name
    combined_context['version'] = __version__
    combined_context["sample_type"] = report_mode
    combined_context["configured_sample_type"] = sample_type
    combined_context["_run_config"] = config
    combined_context["vdj_t_target_enabled"] = "false"
    combined_context["vdj_b_target_enabled"] = "false"
    combined_context["fastq_display_html"] = ""
    combined_context["logo_data_uri"] = _load_logo_data_uri(template_dir)
    if config.get("_analysis_failed"):
        error_text = str(config.get("_analysis_error") or "Unknown pipeline error")
        combined_context["analysis_status_html"] = (
            '<section class="analysis-warning"><strong>Incomplete analysis:</strong> '
            'the pipeline failed before all stages finished. This partial report only '
            'contains outputs completed before the failure. '
            f'<br><strong>Last error:</strong> {html.escape(error_text)}</section>'
        )
    else:
        combined_context["analysis_status_html"] = ""

    # Add input CSV and config parameters
    combined_context['input_csv_data'] = config.get('csv_content', '')
    config_for_display = config.copy()
    config_for_display.pop('csv_content', None)
    combined_context['config_parameters'] = json.dumps(config_for_display, indent=4, default=str)

    # Handle Plotly JS loading
    plotly_candidates = [
        template_dir / 'plotly.js',
        template_dir / 'plotly-2.26.0.min.js',
        template_dir / 'plotly.min.js'
    ]
    raw_plotly_js = ''
    for js_path in plotly_candidates:
        if js_path.exists():
            try:
                with open(js_path, 'r', encoding='utf-8') as f:
                    raw_plotly_js = f.read()
                break
            except Exception:
                pass
    if raw_plotly_js:
        combined_context['plotly_loader_tag'] = '<script>' + raw_plotly_js + '</script>'
    else:
        # Fallback to CDN if local file not found
        combined_context['plotly_loader_tag'] = '<script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>'

    _process_rna_data(omics_data, sample_outdir, combined_context)

    # Reports can also be generated from a legacy or partial run config that
    # has no sequence_files entry.  The barcode report data is still
    # recoverable from stats.tsv and expect_id_barcode.tsv; populate it here so
    # the summary-card area does not silently render as an empty gap.
    if "barcode_mode_cards_data" not in combined_context:
        _process_rna_stats_table_data(sample_outdir, combined_context)
        _process_barcode_report_data(sample_outdir, combined_context, config)

    identifiers = set(re.findall(r"\$\{([a-zA-Z0-9_]+)\}", template_str))
    for key in identifiers:
        if key in _JS_TEMPLATE_PLACEHOLDERS:
            continue
        if key not in combined_context:
            combined_context[key] = _default_placeholder_value(key)

    report_html = template.safe_substitute(combined_context)

    sample_type = combined_context.get("sample_type") or "analysis"
    if config.get("_analysis_failed"):
        report_suffix = "partial_report"
    elif sample_type == "auto":
        report_suffix = "auto_plate_report"
    elif sample_type == "manual":
        report_suffix = "manual_report"
    elif sample_type == "custom":
        report_suffix = "custom_barcode_report"
    else:
        report_suffix = "analysis_report"
    out_file = outs_dir / f'{name}_{report_suffix}.html'
    tmp_file = out_file.with_suffix(out_file.suffix + ".tmp")
    try:
        with open(tmp_file, 'w', encoding='utf-8') as f:
            f.write(report_html)
            f.flush()
            os.fsync(f.fileno())

        if not config.get("_analysis_failed"):
            # run_manifest.json is the completion marker. Remove the previous
            # marker before touching exported data so an interrupted refresh
            # cannot be mistaken for a fully published run.
            success_manifest = outs_dir / "run_manifest.json"
            if success_manifest.exists():
                success_manifest.unlink()
            export_deliverables_to_outs(sample_outdir, outs_dir, name)
            export_diagnostics_to_outs(sample_outdir, outs_dir, name)
            os.replace(tmp_file, out_file)
            for stale_report in outs_dir.glob(f"{name}_*_report.html"):
                if stale_report != out_file:
                    stale_report.unlink()
            partial_manifest = outs_dir / "partial_run_manifest.json"
            if partial_manifest.exists():
                partial_manifest.unlink()
            write_run_manifest(outs_dir, name, config, "success", out_file)
        else:
            # A partial report is useful even if copying diagnostics later
            # fails, and it never replaces the last successful report.
            os.replace(tmp_file, out_file)
            export_diagnostics_to_outs(sample_outdir, outs_dir, name, include_config=True)
            write_run_manifest(outs_dir, name, config, "failed", out_file)
            logger.info("Analysis is incomplete; only diagnostic snapshots were exported.")
    finally:
        if tmp_file.exists():
            tmp_file.unlink()

    logger.info(f"HTML report saved to: {out_file}")
    logger.info("HTML report generation complete.")
    return out_file

#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

GA4GH_METRIC_BASE = "https://ga4gh.github.io/quality-control-wgs/metrics_definitions/"

# ---- Map JSON keys → GA4GH URL slugs ----
# Key names must match exactly as they appear in your input JSON.
METRIC_URLS = {
    "mean_autosome_coverage": "mean-autosome-coverage/",
    "mad_autosome_coverage": "mad-autosome-coverage/",
    "pct_autosomes_15x": "pct-autosomes-15x/",
    "pct_reads_mapped": "pct-reads-mapped/",
    "pct_reads_properly_paired": "pct-reads-properly-paired/",
    "insert_size_std_deviation": "insert-size-std-deviation/",
    "mean_insert_size": "mean-insert-size/",
    "cross_contamination_rate": "cross-contamination-rate/",
    "yield_bp_q30": "yield-bp-q30/",

    "count_snvs": "count-snvs/",
    "count_insertions": "count-insertions/",
    "count_deletions": "count-deletions/",
    "ratio_insertion_deletion": "ratio-insertion-deletion/",
    "ratio_transitions_transversions_snv": "ratio-transitions-transversions-snv/",
    "ratio_heterozygous_homzygous_snv": "ratio-heterozygous-homzygous-snv/",
    "ratio_heterozygous_homzygous_indel": "ratio-heterozygous-homzygous-indel/",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", required=True, help="Sample ID")
    parser.add_argument("--input_aln_metrics", required=True)
    parser.add_argument("--input_variants_metrics", required=True)
    parser.add_argument("--output_jsonld", default="./metrics.jsonld",
                        help="Output JSON-LD file path")
    parser.add_argument("--scratch_dir", default="./")
    return parser.parse_args()


def load_json(path):
    with open(path) as f:
        return json.load(f)


def inject_jsonld_context(merged_dict):
    """
    Converts metrics to JSON-LD by attaching @context with GA4GH URLs
    """
    context = {}

    # Populate @context for every recognized metric
    qc = merged_dict.get("wgs_qc_metrics", {})
    aln = qc.get("aln_metrics", {})
    var = qc.get("variant_metrics", {})

    for k in list(aln.keys()) + list(var.keys()):
        if k in METRIC_URLS:
            context[k] = GA4GH_METRIC_BASE + METRIC_URLS[k]

    # Add GA4GH schema context root (recommended)
    context["@vocab"] = "https://ga4gh.org/wgs_qc/"

    merged_dict["@context"] = context
    return merged_dict


def save_jsonld(aln, variants, outfile):
    # Merge metric blocks
    aln["wgs_qc_metrics"]["variant_metrics"].update(
        variants["wgs_qc_metrics"]["variant_metrics"]
    )

    # Convert to JSON-LD
    merged = inject_jsonld_context(aln)

    with open(outfile, "w") as f:
        json.dump(merged, f, indent=4, sort_keys=False)
        f.write("\n")


if __name__ == "__main__":
    args = parse_args()
    Path(args.scratch_dir).mkdir(parents=True, exist_ok=True)

    aln = load_json(args.input_aln_metrics)
    variants = load_json(args.input_variants_metrics)
    save_jsonld(aln, variants, args.output_jsonld)

#!/usr/bin/env python3

import logging
import argparse
import json
import metrics


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--multiqc_json", dest="multiqc_json", required=True,
                        default=None,
                        help="Path to multiqc_data.json.")
    parser.add_argument("--output_json", dest="output_json", required=False,
                        default="./metrics.json",
                        help="Path to output json. Default: ./metrics.json")
    parser.add_argument("--version_info", dest="version_info", required=False,
                        default=None,
                        help="Path to version_info file. Default: None")
    parser.add_argument("--loglevel", dest="loglevel", required=False,
                        default="INFO",
                        help="Set logging level to INFO (default), WARNING or DEBUG.")
    args = parser.parse_args()
    return args


def set_logging(loglevel):
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        level=numeric_level)


def load_multi_qc(multiqc_json):
    logging.info("Reading multiqc output: %s", multiqc_json)
    with open(multiqc_json) as json_file:
        d = json.load(json_file)
        # keep raw data to simplify parsing downstream
        mqc = d["report_saved_raw_data"]
    return mqc


def calculate_metrics(mqc):
    metrics_list = [
        # primary metrics
        "yield_raw_gb", "yield_pf_gb", "pct_q30_bases", "pct_q30_bases_read1", "pct_q30_bases_read2",
        "yield_raw_reads", "yield_pf_reads", "pct_pf_reads",
        # alignment
        "pct_aligned_bases", "pct_reads_aligned", "pct_reads_aligned_mapqge20",
        "pct_reads_unaligned", "pct_reads_aligned_in_pairs", "pct_reads_properly_paired",
        "total_alignments", "pct_singleton_alignments", "pct_alignments_diff_chrom",
        "pct_alignments_diff_chrom_mapqge5", "pct_chimeras", "pct_secondary_alignments",
        "pct_supplementary_alignments", "pct_duplicate_reads", "mismatch_rate", "mismatch_rate_mapqge20",
        # coverage
        "genome_territory", "mean_autosome_coverage", "sd_autosome_coverage", "median_autosome_coverage", "mad_autosome_coverage",
        "pct_autosomes_1x", "pct_autosomes_10x", "pct_autosomes_15x", "pct_autosomes_30x", "pct_autosomes_40x", "coverage_sg10k_062017",
        # insert size
        "mean_insert_size", "sd_insert_size", "median_insert_size", "mad_insert_size",
        "pct_adapters", "pct_overlapping_bases",
        # gc bias
        "at_dropout", "gc_dropout", "gc_nc_0_19", "gc_nc_20_39", "gc_nc_40_59",
        "gc_nc_60_79", "gc_nc_80_100",
        # variant calling
        "all_snps", "all_het_snps", "all_homalt_snps", "all_snp_het_hom",
        "pass_snps", "pass_het_snps", "pass_homalt_snps", "pass_snp_het_hom",
        "all_indels", "all_het_indels", "all_homalt_indels", "all_indel_het_hom",
        "pass_indels", "pass_het_indels", "pass_homalt_indels", "pass_indel_het_hom",
        "all_del", "all_ins", "pass_del", "pass_ins",
        "all_mnps", "pass_mnps",
        "all_complex_indels", "all_complex_ins", "all_complex_del",
        "pass_complex_indels", "pass_complex_ins", "pass_complex_del",
        "all_multiallelic_sites", "pass_multiallelic_sites",
        "all_snp_ts_tv", "pass_snp_ts_tv",
        # contamination and pst
        "pct_contamination", "pst_pct_concordance"
    ]
    result = dict()
    logging.info("Calculating %d metrics..." % len(metrics_list))
    for m in metrics_list:
        k, v = eval("metrics." + m + "(mqc)")
        result[k] = v

    return result


def add_version_info(parsed_metrics, version_info):
    version = "NA"
    if version_info is not None:
        with open(version_info) as f:
            version = f.read().strip()
    parsed_metrics["metrics_version"] = version
    return parsed_metrics


def save_output(d, outfile):
    with open(outfile, "w") as f:
        json.dump(d, f, sort_keys=True, indent=4)
        f.write("\n")


def done():
    logging.info("DONE")


if __name__ == "__main__":
    args = parse_args()
    set_logging(args.loglevel)

    mqc = load_multi_qc(args.multiqc_json)
    parsed_metrics = calculate_metrics(mqc)
    final_metrics = add_version_info(parsed_metrics, args.version_info)
    save_output(final_metrics, args.output_json)

    done()

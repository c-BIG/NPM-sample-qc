#!/usr/bin/env python3

import logging
import argparse
import json
import pprint as pp
import metrics


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--multiqc_json", dest="multiqc_json", required=True,
                        default=None,
                        help="Path to multiqc_data.json.")
    parser.add_argument("--output_json", dest="output_json", required=False,
                        default="./metrics.json",
                        help="Path to output json. Default: ./metrics.json")
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
        "yield_raw_reads", "yield_pf_reads", "pct_pf_reads",
        # alignment
        "pct_aligned_bases", "pct_reads_aligned", "pct_reads_aligned_mapqge20",
        "pct_reads_unaligned", "pct_reads_aligned_in_pairs", "pct_reads_properly_paired",
        "total_alignments", "pct_singleton_alignments", "pct_alignments_diff_chrom",
        "pct_alignments_diff_chrom_mapqge5", "pct_chimeras", "pct_secondary_alignments",
        "pct_supplementary_alignments", "pct_duplicate_reads", "mismatch_rate", "mismatch_rate_mapqge20",
        # coverage
        "genome_territory", "mean_coverage", "sd_coverage", "median_coverage", "mad_coverage",
        "pct_1x", "pct_15x", "pct_30x", "pct_40x",
        # insert size
        "mean_insert_size", "sd_insert_size", "median_insert_size", "mad_insert_size",
        "pct_adapters",
        # gc bias
        "at_dropout", "gc_dropout", "gc_nc_0_19", "gc_nc_20_39", "gc_nc_40_59",
        "gc_nc_60_79", "gc_nc_80_100",
        # variant calling
        "all_snps", "pass_snps", "all_indels", "pass_indels", "pass_sites_multiallelic",
        "pass_snps_multiallelic", "pass_mnps", "pass_complex_indels", "snp_ts_tv"
    ]
    result = dict()
    for m in metrics_list:
        r = eval("metrics." + m + "(mqc)")
        result.update(r)
    return result


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
    final_metrics = calculate_metrics(mqc)
    save_output(final_metrics, args.output_json)

    pp.pprint(final_metrics)
    done()

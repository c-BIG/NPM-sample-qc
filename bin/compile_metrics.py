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
    parser.add_argument("--sample_id", dest="sample_id", required=True,
                        default=None,
                        help="Sample ID. Default: None")
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
        "pct_q30_bases",
        # alignment
        "pct_reads_aligned_in_pairs", "pct_reads_properly_paired",
        # coverage
        "mean_autosome_coverage", "mad_autosome_coverage",
        "pct_autosomes_15x",
        # insert size
        "mean_insert_size"
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


def add_sample_id(version_metrics, sample_id):
    sampleid = "NA"
    if sample_id is not None:
        sampleid = sample_id
    version_metrics["sample"] = sampleid
    return version_metrics


def save_output(d, outfile):
    with open(outfile, "w") as f:
        d = {"biosample" : {"id" : args.sample_id}, "qc-metrics": d}
        json.dump(d, f, sort_keys=True, indent=4)
        f.write("\n")


def done():
    logging.info("DONE")


if __name__ == "__main__":
    args = parse_args()
    set_logging(args.loglevel)

    mqc = load_multi_qc(args.multiqc_json)
    parsed_metrics = calculate_metrics(mqc)
    version_metrics = add_version_info(parsed_metrics, args.version_info)
    final_metrics = add_sample_id(version_metrics, args.sample_id)
    save_output(final_metrics, args.output_json)

    done()

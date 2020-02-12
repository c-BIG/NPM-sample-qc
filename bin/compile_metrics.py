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
    result = dict()
    result["perc_mapped_reads"] = metrics.perc_mapped_reads(mqc)
    result["perc_properly_paired"] = metrics.perc_properly_paired(mqc)
    result["perc_diff_chrom_mapqge5"] = metrics.perc_diff_chrom_mapqge5(mqc)
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
    pp.pprint(final_metrics)

    save_output(final_metrics, args.output_json)
    done()

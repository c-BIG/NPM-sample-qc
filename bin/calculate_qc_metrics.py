#!/usr/bin/env python3

import logging
import argparse
import json
import pprint as pp
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--multiqc_json", dest="multiqc_json", required=True, 
                        default=None,
                        help="Path to multiqc_data.json.")
    parser.add_argument("--loglevel", dest="loglevel", required=False,
                        default="INFO",
                        help="Set logging level to INFO (default), WARNING or DEBUG.")
    args = parser.parse_args()
    return args


def set_logging(args):
    numeric_level = getattr(logging, args.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        level=numeric_level)


def load_multi_qc(args):
    logging.info("Reading multiqc output: %s", args.multiqc_json)
    with open(args.multiqc_json) as json_file:
        d = json.load(json_file)
        # keep raw data to simplify parsing downstream
        mqc = d["report_saved_raw_data"]
        return mqc


def calculate_metrics(mqc):
    result = dict()

    # samtools flagstat
    samtools_flagstat = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    # % mapped
    result["perc_mapped_reads"] = samtools_flagstat["mapped_passed_pct"]
    # % properly paired
    result["perc_properly_paired"] = samtools_flagstat["properly paired_passed_pct"]
    # % diff chrom
    x = 100 * np.divide(samtools_flagstat["with mate mapped to a different chr (mapQ >= 5)_passed"], samtools_flagstat["total_passed"])
    result["perc_diff_chrom_mapqge5"] = np.round(x, 2)

    # done
    return result


def done():
    logging.info("DONE")


if __name__ == "__main__":
    args = parse_args()
    set_logging(args)

    mqc = load_multi_qc(args)
    final_metrics = calculate_metrics(mqc)
    pp.pprint(final_metrics)

    done()

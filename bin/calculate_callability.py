#!/usr/bin/env python3

import logging
import argparse
import json
import pprint as pp
import subprocess
import numpy as np
import sys
import os
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_gvcf", dest="input_gvcf", required=True,
                        default=None,
                        help="Path to input gVCF.")
    parser.add_argument("--autosomes_bed", dest="autosomes_bed", required=True,
                        default=None,
                        help="Path to autosome BED.")
    parser.add_argument("--output_json", dest="output_json", required=False,
                        default="./callability.json",
                        help="Path to output file for callability metrics. Default: ./callability.json")
    parser.add_argument("--scratch_dir", dest="scratch_dir", required=False,
                        default="./",
                        help="Path to scratch dir. Default: ./")
    parser.add_argument("--keep_scratch", dest="keep_scratch", required=False,
                        default=False, action="store_true",
                        help="Keep scratch dir. Default: False.")
    parser.add_argument("--loglevel", dest="loglevel", required=False,
                        default="INFO",
                        help="Set logging level to INFO (default), WARNING or DEBUG.")
    args = parser.parse_args()

    # simple check for gvcf inputs
    if "gvcf" not in args.input_gvcf:
        sys.exit("Please provide a gvcf as input (gvcf in the file name)")

    # create scratch dir if it doesn't exist
    Path(args.scratch_dir).mkdir(parents=True, exist_ok=True)

    return args


def set_logging(loglevel):
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        level=numeric_level)


def calculate_callability(input_gvcf, autosomes_bed, scratch_dir):

    logging.info("Calculating pct_autosome_callability...")

    # initialise results dict
    r = dict()
    metrics_list = [
        "pct_autosome_callability"
    ]
    for m in metrics_list:
        r[m] = 0

    # autosome bases
    cmd = "cat %s | awk '{SUM += $3-$2} END {print SUM}'" % autosomes_bed
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    l = p.stdout.read().decode("utf-8").strip()
    total_bases = int(float(l))

    # called bases (pass)
    cmd = "bcftools view -f PASS %s" % input_gvcf
    cmd += " | bedtools intersect -a %s -b stdin > %s/called_bases.bed" % (autosomes_bed, scratch_dir)
    cmd += " ; cat %s/called_bases.bed" % scratch_dir
    cmd += " | awk '{SUM += $3-$2} END {print SUM}'"
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    l = p.stdout.read().decode("utf-8").strip()
    try:
        called_bases = int(l)
    except ValueError:
        called_bases = 0

    # autosome callability
    autosome_callability = np.divide(called_bases, total_bases)
    r["pct_autosome_callability"] = np.round(autosome_callability*100, 2)

    return r


def save_output(d, outfile):
    with open(outfile, "w") as f:
        json.dump(d, f, sort_keys=True, indent=4)
        f.write("\n")


def done(scratch_dir):
    if not args.keep_scratch:
        # clean up
        logging.info("Cleaning up...")
        os.remove("%s/called_bases.bed" % scratch_dir)
        if len(os.listdir(scratch_dir)) == 0:
            os.rmdir(scratch_dir)
        else:
            logging.warning("Couldn't remove scratch dir (%s): directory is not empty" % scratch_dir)
    logging.info("DONE")


if __name__ == "__main__":
    args = parse_args()
    set_logging(args.loglevel)

    callability = calculate_callability(args.input_gvcf, args.autosomes_bed, args.scratch_dir)
    save_output(callability, args.output_json)

    done(args.scratch_dir)

#!/usr/bin/env python3

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
    parser.add_argument("--sample_id", dest="sample_id", required=True,
                        default=None,
                        help="Sample ID")
    parser.add_argument("--input_aln_metrics", dest="input_aln_metrics", required=True,
                        default=None,
                        help="Path to input aln metrics list")
    parser.add_argument("--input_variants_metrics", dest="input_variants_metrics", required=True,
                        default=None,
                        help="Path to input variants metrics list")
    parser.add_argument("--output_json", dest="output_json", required=False,
                        default="./metrics.json",
                        help="Path to output file for variant metrics. Default: ./variant_counts.json")
    parser.add_argument("--scratch_dir", dest="scratch_dir", required=False,
                        default="./",
                        help="Path to scratch dir. Default: ./")
    args = parser.parse_args()

    # create scratch dir if it doesn't exist
    Path(args.scratch_dir).mkdir(parents=True, exist_ok=True)

    return args


def data1(input_aln_metrics):
    aln_input = {}
    with open(input_aln_metrics, 'r') as f1:
        aln_input = json.load(f1)
    return aln_input

def data2(input_variants_metrics):
    variants_input = {}
    with open(input_variants_metrics, 'r') as f2:
        variants_input = json.load(f2)
    return variants_input


def save_output(data1, data2, outfile):
    with open(outfile, "w") as f:
        data1['wgs_qc_metrics']['variant_metrics'].update(data2['wgs_qc_metrics']['variant_metrics'])
        print(data1)
        json.dump(data1, f, sort_keys=True, indent=4)
        #json.dump({"biosample": data1["biosample"], "wgs_qc_metrics": {**data1["wgs_qc_metrics"], **data2["wgs_qc_metrics"]}}, f, sort_keys=True, indent=4)
        f.write("\n")

if __name__ == "__main__":
    args = parse_args()

    aln = data1(args.input_aln_metrics)
    variants = data2(args.input_variants_metrics)
    save_output(aln, variants, args.output_json)

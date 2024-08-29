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
    parser.add_argument("--input_metrics", dest="input_metrics", required=True,
                        default=None,
                        help="Path to input aln metrics list")
    parser.add_argument("--output_json", dest="output_json", required=False,
                        default="./variant_counts.json",
                        help="Path to output file for variant metrics. Default: ./variant_counts.json")
    parser.add_argument("--scratch_dir", dest="scratch_dir", required=False,
                        default="./",
                        help="Path to scratch dir. Default: ./")
    args = parser.parse_args()

    # create scratch dir if it doesn't exist
    Path(args.scratch_dir).mkdir(parents=True, exist_ok=True)

    return args

def raw_data(input_metrics):
    d = {}
    with open(input_metrics) as f:
        for line in f:
            if not line.strip():
                continue
            row = line.split('\t')
            key = row[0]
            value_str = row[1]
            try:
                value = float(value_str.strip())
            except ValueError:
                value = value_str.strip()
            d[key] = value
    return d


def save_output(data_metrics, outfile):
    with open(outfile, "w") as f:
        data_metrics = {"sample" : {"id" : args.sample_id}, "wgs_metrics" : data_metrics}
        json.dump(data_metrics, f, sort_keys=True, indent=4)
        f.write("\n")

if __name__ == "__main__":
    args = parse_args()

    data_metrics = raw_data(args.input_metrics)
    save_output(data_metrics, args.output_json)

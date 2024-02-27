#!/usr/bin/env python3

import argparse
import logging


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", dest="input", required=True,
                        default=None,
                        help="Path to picard CollectWgsMetrics file.")
    args = parser.parse_args()
    return args


# Replace metrics class name
# This temporary fix is needed due to compatibility issue with GATK
# 4.1.2.0 and MultiQC 1.8. MultiQC searches for CollectWgsMetrics
# instead of WgsMetrics in the report file. Note that this should be
# resolved in the 1.9 release of MultiQC.
# From https://github.com/genialis/resolwe-bio/blob/master/resolwe_bio/processes/support_processors/wgs_metrics.py#L19
def replace_metrics_class(fname):
    with open(fname, 'r') as report:
        newlines = []
        for line in report.readlines():
            if line == '## METRICS CLASS\tpicard.analysis.WgsMetrics\n':
                line = '## METRICS CLASS\tCollectWgsMetrics$WgsMetrics\n'
                newlines.append(line)
            else:
                newlines.append(line)
    with open(fname, 'w') as report:
        for line in newlines:
            report.writelines(line)


if __name__ == "__main__":
    args = parse_args()
    replace_metrics_class(args.input)
    print("DONE: %s" % args.input)

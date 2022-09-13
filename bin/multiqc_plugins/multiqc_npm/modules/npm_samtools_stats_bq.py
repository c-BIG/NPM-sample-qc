#!/usr/bin/env python3

"""
Parser for samtools stats
"""

import logging
import pandas as pd
import numpy as np

from multiqc.utils import report

log = logging.getLogger(__name__)


def calc_pct_q30_bases(x):
    df = pd.DataFrame(x).apply(pd.to_numeric)
    df.columns = ["cycle_number"] + ["bq_%d" % i for i in range(0, len(df.columns)-1, 1)]

    bq_all = df.loc[:, 'bq_0':]
    bq_ge30 = df.loc[:, 'bq_30':]

    sum_bases_q30 = int(bq_ge30.values.sum())
    sum_bases_all = int(bq_all.values.sum())
    pct_q30_bases = np.round(np.divide(sum_bases_q30, sum_bases_all)*100, 2)

    return pct_q30_bases, sum_bases_q30, sum_bases_all


def parse_reports(self):

    # Set up vars
    self.samtools_stats_bq = dict()

    # Collect metrics
    for f in self.find_log_files('multiqc_npm/samtools_stats_bq', filehandles=True):

        # Collect relevant records and calculate metrics
        parsed_data = dict()
        ffq = list()
        lfq = list()

        for l in f['f']:
            if l.startswith("FFQ"):
                vals = l.strip("\n").split("\t")
                ffq.append(vals[1:])
            elif l.startswith("LFQ"):
                vals = l.strip("\n").split("\t")
                lfq.append(vals[1:])
            else:
                continue

        pct, total, q30 = calc_pct_q30_bases(ffq)
        parsed_data["pct_q30_bases_read1"] = pct
        parsed_data["total_bases_read1"] = total
        parsed_data["q30_bases_read1"] = q30

        pct, total, q30 = calc_pct_q30_bases(lfq)
        parsed_data["pct_q30_bases_read2"] = pct
        parsed_data["total_bases_read2"] = total
        parsed_data["q30_bases_read2"] = q30

        # Save results
        s_name = f["s_name"]
        self.samtools_stats_bq[s_name] = parsed_data

    # Write results
    if len(self.samtools_stats_bq) > 0:

        # Write parsed data to a file
        self.write_data_file(self.samtools_stats_bq, 'multiqc_npm_samtools_stats_bq')

    # Return the number of detected samples to the parent module
    return len(self.samtools_stats_bq)

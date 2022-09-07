#!/usr/bin/env python3

"""
Parser for run_mosdepth.sh
"""

import logging
import csv

from multiqc.utils import config, report

log = logging.getLogger(__name__)


def parse_reports(self):

    # Set up vars
    self.mosdepth = dict()

    # Collect metrics
    for f in self.find_log_files('multiqc_npm/mosdepth', filehandles=True):

        # Collect relevant records and calculate metrics
        vals = list()
        for l in f["f"]:
            v = l.strip("\n").split(",")
            vals.append(v)
        parsed_data = dict(zip(vals[0], vals[1]))

        # Save results
        s_name = f["s_name"]
        self.mosdepth[s_name] = parsed_data

    # Write results
    if len(self.mosdepth) > 0:

        # Write parsed data to a file
        self.write_data_file(self.mosdepth, 'multiqc_npm_mosdepth')

    # Return the number of detected samples to the parent module
    return len(self.mosdepth)

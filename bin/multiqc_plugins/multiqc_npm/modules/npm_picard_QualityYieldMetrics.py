#!/usr/bin/env python3

"""
Parser for picard QualityYieldMetrics
"""

import logging
import os
import re

from multiqc.utils import config, report

log = logging.getLogger(__name__)


def parse_reports(self):

    # Set up vars
    self.picard_quality_yield_metrics = dict()

    # Collect metrics
    for f in self.find_log_files('multiqc_npm/picard_quality_yield_metrics', filehandles=True):
        parsed_data = dict()
        s_name = None
        keys = None

        for l in f['f']:

            # Pull sample name from input
            if 'QualityYieldMetrics' in l and 'INPUT' in l:
                fn_search = re.search(r"INPUT(?:=|\s+)(\[?[^\s]+\]?)", l, flags=re.IGNORECASE)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1).strip('[]'))
                    s_name = self.clean_s_name(s_name, f['root'])

            # Parse metric names and values
            if s_name is not None:
                if 'QualityYieldMetrics' in l and '## METRICS CLASS' in l:
                    keys = f['f'].readline().strip("\n").split("\t")
                elif keys:
                    vals = l.strip("\n").split("\t")
                    if len(vals) == len(keys):
                        for i, k in enumerate(keys):
                            try:
                                parsed_data[k] = float(vals[i])
                            except ValueError:
                                parsed_data[k] = vals[i]

        # Save results
        self.picard_quality_yield_metrics[s_name] = parsed_data

    # Write results
    if len(self.picard_quality_yield_metrics) > 0:

        # Write parsed data to a file
        self.write_data_file(self.picard_quality_yield_metrics, 'multiqc_npm_picard_QualityYieldMetrics')

    # Return the number of detected samples to the parent module
    return len(self.picard_quality_yield_metrics)

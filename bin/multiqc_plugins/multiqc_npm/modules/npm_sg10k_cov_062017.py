#!/usr/bin/env python3

"""
Parser for sg10k-cov-062017.sh
"""

import logging

from multiqc.utils import config, report

log = logging.getLogger(__name__)


def parse_reports(self):

    # Set up vars
    self.sg10k_cov_062017 = dict()

    # Collect metrics
    for f in self.find_log_files("multiqc_npm/sg10k_cov_062017", filehandles=True):

        # Collect relevant records and calculate metrics
        parsed_data = dict()

        for l in f["f"]:
            # file contains only 1 line
            parsed_data["bases_sg10k_062017"] = int(l.strip("\n"))

        # Save results
        s_name = f["s_name"]
        self.sg10k_cov_062017[s_name] = parsed_data

    # Write results
    if len(self.sg10k_cov_062017) > 0:

        # Write parsed data to a file
        self.write_data_file(self.sg10k_cov_062017, "multiqc_npm_sg10k_cov_062017")

    # Return the number of detected samples to the parent module
    return len(self.sg10k_cov_062017)

#!/usr/bin/env python3

"""
Parser for count_variants.py
"""

import logging
import json

from multiqc.utils import report

log = logging.getLogger(__name__)


def parse_reports(self):

    # Set up vars
    self.count_variants = dict()

    # Collect metrics
    for f in self.find_log_files('multiqc_npm/count_variants'):
        parsed_data = json.loads(f["f"])
        # Save results
        s_name = f["s_name"]
        self.count_variants[s_name] = parsed_data

    # Write results
    if len(self.count_variants) > 0:

        # Write parsed data to a file
        self.write_data_file(self.count_variants, 'multiqc_npm_count_variants')

    # Return the number of detected samples to the parent module
    return len(self.count_variants)

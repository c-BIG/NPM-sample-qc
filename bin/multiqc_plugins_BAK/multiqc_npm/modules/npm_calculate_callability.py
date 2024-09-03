#!/usr/bin/env python3

"""
Parser for calculate_callability.py
"""

import logging
import json

from multiqc.utils import report

log = logging.getLogger(__name__)


def parse_reports(self):

    # Set up vars
    self.calculate_callability = dict()

    # Collect metrics
    for f in self.find_log_files('multiqc_npm/calculate_callability'):
        parsed_data = json.loads(f["f"])
        # Save results
        s_name = f["s_name"]
        self.calculate_callability[s_name] = parsed_data

    # Write results
    if len(self.calculate_callability) > 0:

        # Write parsed data to a file
        self.write_data_file(self.calculate_callability, 'multiqc_npm_calculate_callability')

    # Return the number of detected samples to the parent module
    return len(self.calculate_callability)

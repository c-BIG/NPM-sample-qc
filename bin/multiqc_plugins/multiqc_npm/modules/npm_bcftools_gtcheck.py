#!/usr/bin/env python3

"""
Parser for bcftools gtcheck
"""

import logging
import numpy as np

log = logging.getLogger(__name__)


def parse_reports(self):

    # Set up vars
    self.bcftools_gtcheck = dict()

    # Collect metrics
    for f in self.find_log_files('multiqc_npm/bcftools_gtcheck', filehandles=True):

        # Collect relevant records and calculate metrics
        parsed_data = dict()

        for l in f["f"]:
            if l.startswith("CN"):
                vals = l.strip("\n").split("\t")
            else:
                continue

        discordant_sites = int(float(vals[1]))
        # gtcheck considers snps + indels if available
        total_sites = int(float(vals[3]))
        frac_concordant_sites = 1 - np.divide(discordant_sites, total_sites)
        pct_concordance = np.round(frac_concordant_sites*100, 2)

        parsed_data["pst_pct_concordance"] = pct_concordance
        parsed_data["pst_sites_compared"] = total_sites

        # Save results
        s_name = f["s_name"].replace(".bcftools_gtcheck", "")
        self.bcftools_gtcheck[s_name] = parsed_data

    # Write results
    if len(self.bcftools_gtcheck) > 0:

        # Write parsed data to a file
        self.write_data_file(self.bcftools_gtcheck, 'multiqc_npm_bcftools_gtcheck')

    # Return the number of detected samples to the parent module
    return len(self.bcftools_gtcheck)

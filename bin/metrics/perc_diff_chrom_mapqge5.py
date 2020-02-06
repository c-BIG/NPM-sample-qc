#!/usr/bin/env python3

import numpy as np


def perc_diff_chrom_mapqge5(mqc):
    """
    Percent diff chrom mapQ >= 5
    """
    samtools_flagstat = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    m = 100 * np.divide(samtools_flagstat["with mate mapped to a different chr (mapQ >= 5)_passed"], samtools_flagstat["total_passed"])
    m = np.round(m, 2)
    return m

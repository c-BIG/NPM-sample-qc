#!/usr/bin/env python3

import numpy as np


def perc_properly_paired(mqc):
    """
    Percent properly paired reads
    """
    samtools_flagstat = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    m = samtools_flagstat["properly paired_passed_pct"]
    m = np.round(m, 2)
    return m

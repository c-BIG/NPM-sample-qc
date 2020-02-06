#!/usr/bin/env python3

import numpy as np


def perc_mapped_reads(mqc):
    """
    Percent mapped reads
    """
    samtools_flagstat = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    m = samtools_flagstat["mapped_passed_pct"]
    m = np.round(m, 2)
    return m

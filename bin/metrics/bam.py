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


def perc_properly_paired(mqc):
    """
    Percent properly paired reads
    """
    samtools_flagstat = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    m = samtools_flagstat["properly paired_passed_pct"]
    m = np.round(m, 2)
    return m


def perc_diff_chrom_mapqge5(mqc):
    """
    Percent diff chrom mapQ >= 5
    """
    samtools_flagstat = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    m = 100 * np.divide(samtools_flagstat["with mate mapped to a different chr (mapQ >= 5)_passed"], samtools_flagstat["total_passed"])
    m = np.round(m, 2)
    return m

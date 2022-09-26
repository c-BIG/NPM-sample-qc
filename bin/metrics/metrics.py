#!/usr/bin/env python3

import numpy as np
import inspect


DECIMALS = 5


def pct_q30_bases(mqc):
    """
    The percentage of PF bases with base quality >= 30.

    - Minimum mapping quality 0
    - Reads marked as duplicates not excluded
    - Secondary alignments excluded
    - Alignments are hard and soft clipped
    - Overlapping bases not excluded
    - Unique Molecular Identifiers not used to collapse reads

    Source: picard CollectMultipleMetrics(QualityYieldMetrics) - PF_Q30_BASES/PF_BASES
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_picard_QualityYieldMetrics"].values()))
        v = np.round(np.divide(d["PF_Q30_BASES"],
                               d["PF_BASES"])*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_reads_aligned_in_pairs(mqc):
    """
    The percentage of reads that have been aligned as pairs, i.e. percent aligned pairs * 2.

    - Genomic reference build GRCh38
    - Genomic filtering regions - excluded non autosomal regions
    - Minimum base quality not applied
    - Minimum mapping quality 0
    - Reads marked as duplicates not excluded
    - Secondary alignments excluded
    - Alignments are hard and soft clipped
    - Overlapping bases not excluded
    - Unique Molecular Identifiers not used to collapse reads

    Source: samtools stats (reads_mapped_and_paired / sequences)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_stats"].values()))
        v = np.divide(d["reads_mapped_and_paired"],
                      d["sequences"])
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_reads_properly_paired(mqc):
    """
    The percentage of reads that have been aligned as proper pairs.

    - Genomic reference build GRCh38
    - Genomic filtering regions - excluded non autosomal regions
    - Minimum base quality not applied
    - Minimum mapping quality 0
    - Reads marked as duplicates not excluded
    - Secondary alignments excluded
    - Alignments are hard and soft clipped
    - Overlapping bases included
    - Unique Molecular Identifiers not used to collapse reads

    Source: samtools stats (reads_properly_paired_percent)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_stats"].values()))
        v = d["reads_properly_paired_percent"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def mean_autosome_coverage(mqc):
    """
    The mean coverage in autosomes (The number of non-N bases in autosomes over which coverage will be evaluated).

    - Genomic reference build GRCh38
    - Genomic filtering regions - non autosomes and gap region (N bases) in the GRCh38 assembly
    - Minimum base quality not applied
    - Minimum mapping quality 20
    - Reads marked as duplicates not excluded
    - Secondary alignments excluded
    - Alignments are hard and soft clipped
    - Overlapping bases excluded
    - Unique Molecular Identifiers not used to collapse reads

    Source: in-house tool based on mosdepth v0.3.3 (mean_autosome_coverage)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = d["mean_autosome_coverage"]
        v = np.round(np.float(v), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def mad_autosome_coverage(mqc):
    """
    The median absolute deviation of coverage in autosomes, after coverage filters are applied (see mean_autosome_coverage for details).

    - Genomic reference build GRCh38
    - Genomic filtering regions - non autosomes and gap region (N bases) in the GRCh38 assembly
    - Minimum base quality not applied
    - Minimum mapping quality 20
    - Reads marked as duplicates not excluded
    - Secondary alignments excluded
    - Alignments are hard and soft clipped
    - Overlapping bases excluded
    - Unique Molecular Identifiers not used to collapse reads

    Source: in-house tool based on mosdepth v0.3.3 (mad_autosome_coverage)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = d["mad_autosome_coverage"]
        v = np.round(np.float(v), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_autosomes_15x(mqc):
    """
    The percentage of bases that attained at least 15X sequence coverage in autosomes, after coverage filters are applied (see mean_autosome_coverage for details).

    - Genomic reference build GRCh38
    - Genomic filtering regions - non autosomes and gap region (N bases) in the GRCh38 assembly
    - Minimum base quality not applied
    - Minimum mapping quality 20
    - Reads marked as duplicates not excluded
    - Secondary alignments excluded
    - Alignments are hard and soft clipped
    - Overlapping bases excluded
    - Unique Molecular Identifiers not used to collapse reads

    Source: in-house tool based on mosdepth v0.3.3 (ge_15x_autosome_bases/total_autosome_bases*100)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = np.divide(np.float(d["ge_15x_autosome_bases"]),
                      np.float(d["total_autosome_bases"]))
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def mean_insert_size(mqc):
    """
    The mean insert size over the "core" of the distribution.

    Note: Artefactual outliers in the distribution often cause calculation of nonsensical mean and
    stdev values. To avoid this the distribution is first trimmed to a "core" distribution of +/- N
    median absolute deviations around the median insert size.

    - Reads marked as duplicates are excluded
    - Minimum base quality not applied
    - Minimum mapping quality 0 
    - Secondary alignments excluded
    - Alignments hard and soft clipped
    - Overlapping bases not excluded
    - Unique Molecular Identifiers not used to collapse reads

    Source: picard CollectMultipleMetrics(InsertSizeMetrics) - MEAN_INSERT_SIZE
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_insertSize"].values()))
        v = d["MEAN_INSERT_SIZE"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

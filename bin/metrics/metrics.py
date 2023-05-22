#!/usr/bin/env python3

import numpy as np
import inspect


DECIMALS = 5


def yield_bp_q30(mqc):
    """
    The number of bases that pass filter (PF) and with base quality ≥ 30 (BQ).

    Source: picard QualityYieldMetrics (PF_Q30_BASES)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_picard_QualityYieldMetrics"].values()))
        v = d["PF_Q30_BASES"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pct_reads_mapped(mqc):
    """
    The percentage of primary reads, paired or single, that are mappable to the REF sequence with MAPQ > 0 after alignment.

    Source: samtools stats (reads_mapped / sequences)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_stats"].values()))
        v = np.divide(d["reads_mapped"],
                      d["sequences"])
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_reads_properly_paired(mqc):
    """
    The percentage of reads that have been aligned as proper pairs.

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
    The mean coverage in autosomes (as defined in genome_territory). Excludes:

    - bases in reads with low mapping quality (mapq < 20)
    - bases in reads marked as duplicates
    - overlapping bases in read pairs

    Source: run_mosdepth.sh (mean_autosome_coverage)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = d["mean_autosome_coverage"]
        v = np.round(np.float64(v), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def mad_autosome_coverage(mqc):
    """
    The median absolute deviation of coverage in autosomes, after coverage filters are applied (see mean_autosome_coverage for details).

    Source: run_mosdepth.sh (mad_autosome_coverage)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = d["mad_autosome_coverage"]
        v = np.round(np.float64(v), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_autosomes_15x(mqc):
    """
    The percentage of bases with at least 15X coverage in autosomes, after coverage filters are applied (see mean_autosome_coverage for details).

    Source: run_mosdepth.sh (ge_15x_autosome_bases * 100 / total_autosome_bases)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = np.divide(np.float64(d["ge_15x_autosome_bases"]),
                      np.float64(d["total_autosome_bases"]))
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def mean_insert_size(mqc):
    """
    The average absolute template length for paired and mapped reads.

    Source: samtools stats (insert_size_average)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_stats"].values()))
        v = d["insert_size_average"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def cross_contamination_rate(mqc):
    """
    Estimation of cross-individual contamination rate.
    Source: VerifyBamID2 (FREEMIX)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_verifybamid"].values()))
        v = d["FREEMIX"]
    except KeyError:
        v = "NA"

    return k, v

#!/usr/bin/env python3

import numpy as np
import inspect


def yield_raw_reads(mqc):
    """
    The total number of reads including all PF and non-PF reads.

    Source: picard AlignmentSummaryMetrics (TOTAL_READS)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
    v = d["TOTAL_READS"]
    v = int(v)

    r = {k:v}
    return r


def yield_pf_reads(mqc):
    """
    The number of PF reads where PF is defined as passing Illumina's filter.

    Source: picard AlignmentSummaryMetrics (PF_READS)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
    v = d["PF_READS"]
    v = int(v)

    r = {k:v}
    return r


def pct_pf_reads(mqc):
    """
    The percentage of reads that are PF, i.e. yield_pf_reads / yield_raw_reads.

    Source: picard AlignmentSummaryMetrics (PCT_PF_READS)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
    v = d["PCT_PF_READS"]
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def pct_aligned_bases(mqc):
    """
    The percentage of bases that are aligned to the reference.

    Source: samtools stats (bases mapped cigar / total length)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_samtools_stats"].values()))
    v = np.divide(d["bases_mapped_(cigar)"],
                  d["total_length"])
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def pct_reads_aligned(mqc):
    """
    The percentage of PF reads that aligned to the reference.

    Source: picard AlignmentSummaryMetrics (PCT_PF_READS_ALIGNED)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
    v = d["PCT_PF_READS_ALIGNED"]
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def pct_reads_aligned_mapqge20(mqc):
    """
    The percentage of PF reads that were aligned to the reference with a mapping quality of Q20 or higher.

    Source: picard AlignmentSummaryMetrics (PF_HQ_ALIGNED_READS / TOTAL_READS)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
    v = np.divide(d["PF_HQ_ALIGNED_READS"],
                  d["TOTAL_READS"])
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def pct_reads_unaligned(mqc):
    """
    The percentage of unmapped reads.

    Source: samtools stats (reads_unmapped / sequences)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_samtools_stats"].values()))
    v = np.divide(d["reads_unmapped"],
                  d["sequences"])
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def pct_reads_aligned_in_pairs(mqc):
    """
    The percentage of reads that have been aligned as pairs, i.e. percent aligned pairs * 2.

    Source: samtools stats (reads_mapped_and_paired / sequences)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_samtools_stats"].values()))
    v = np.divide(d["reads_mapped_and_paired"],
                  d["sequences"])
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def pct_reads_properly_paired(mqc):
    """
    The percentage of reads that have been aligned as proper pairs.

    Source: samtools stats (reads_properly_paired_percent)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_samtools_stats"].values()))
    v = d["reads_properly_paired_percent"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def total_alignments(mqc):
    """
    The total number of alignments in the BAM/CRAM file.
    Includes primary, secondary and supplementary alignments; excludes duplicates.

    Source: samtools flagstat (total_passed)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    v = int(d["total_passed"])

    r = {k:v}
    return r


def pct_singleton_alignments(mqc):
    """
    The percentage of alignments that are singletons.

    Source: samtools flagstat (singletons_passed_pct)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    v = d["singletons_passed_pct"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def pct_alignments_diff_chrom(mqc):
    """
    The percentage of alignments with mates aligned to different chromosomes.

    Source: samtools flagstat (with mate mapped to a different chr_passed / total_passed)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    v = np.divide(d["with mate mapped to a different chr_passed"],
                  d["total_passed"])
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def pct_alignments_diff_chrom_mapqge5(mqc):
    """
    The percentage of alignments with MAPQ>=5 and with mates aligned to different chromosomes.

    Source: samtools flagstat (with mate mapped to a different chr (mapQ >= 5)_passed / total_passed)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    v = np.divide(d["with mate mapped to a different chr (mapQ >= 5)_passed"],
                  d["total_passed"])
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def pct_chimeras(mqc):
    """
    The fraction of reads that map outside of a maximum insert size (usually 100kb)
    or that have the two ends mapping to different chromosomes.

    Source: picard AlignmentSummaryMetrics (PCT_CHIMERAS)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
    v = d["PCT_CHIMERAS"]
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def pct_secondary_alignments(mqc):
    """
    The percentage of alignments marked as secondary.

    Source: samtools flagstat (secondary_passed / total_passed)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    v = np.divide(d["secondary_passed"],
                  d["total_passed"])
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def pct_supplementary_alignments(mqc):
    """
    The percentage of alignments marked as supplementary.

    Source: samtools flagstat (supplementary_passed / total_passed)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
    v = np.divide(d["supplementary_passed"],
                  d["total_passed"])
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def pct_duplicate_reads(mqc):
    """
    The percentage of reads marked as duplicates.

    Source: samtools stats (reads_duplicated_percent)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_samtools_stats"].values()))
    v = d["reads_duplicated_percent"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def mismatch_rate(mqc):
    """
    The fraction of bases mismatching the reference for all bases aligned to the reference sequence.

    Source: picard AlignmentSummaryMetrics (PF_MISMATCH_RATE)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
    v = d["PF_MISMATCH_RATE"]
    v = np.round(v, 4)

    r = {k:v}
    return r


def mismatch_rate_mapqge20(mqc):
    """
    The fraction of bases that mismatch the reference in reads with MAPQ>=20.

    Source: picard AlignmentSummaryMetrics (PF_HQ_ERROR_RATE)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
    v = d["PF_HQ_ERROR_RATE"]
    v = np.round(v, 4)

    r = {k:v}
    return r


def mean_insert_size(mqc):
    """
    The mean insert size over the "core" of the distribution.
    Note: Artefactual outliers in the distribution often cause calculation of nonsensical mean and
    stdev values. To avoid this the distribution is first trimmed to a "core" distribution of +/- N
    median absolute deviations around the median insert size.

    Source: picard InsertSizeMetrics (MEAN_INSERT_SIZE)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_insertSize"].values()))
    v = d["MEAN_INSERT_SIZE"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def sd_insert_size(mqc):
    """
    The standard deviation of insert sizes over the "core" of the distribution.

    Source: picard InsertSizeMetrics (STANDARD_DEVIATION)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_insertSize"].values()))
    v = d["STANDARD_DEVIATION"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def median_insert_size(mqc):
    """
    The median insert size over the "core" of the distribution.

    Source: picard InsertSizeMetrics (MEDIAN_INSERT_SIZE)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_insertSize"].values()))
    v = d["MEDIAN_INSERT_SIZE"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def mad_insert_size(mqc):
    """
    The median absolute deviation of insert sizes over the "core" of the distribution.

    Source: picard InsertSizeMetrics (MEDIAN_ABSOLUTE_DEVIATION)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_insertSize"].values()))
    v = d["MEDIAN_ABSOLUTE_DEVIATION"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def pct_adapters(mqc):
    """
    The fraction of PF reads that are unaligned and match to a known adapter sequence right
    from the start of the read.

    Source: picard AlignmentSummaryMetrics (PCT_ADAPTER)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
    v = d["PCT_ADAPTER"]
    v = np.round(v*100, 2)

    r = {k:v}
    return r


def at_dropout(mqc):
    """
    Illumina-style AT dropout metric. Calculated by taking each GC bin independently
    and calculating (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[0..50].

    Source: picard GcBiasSummaryMetrics (AT_DROPOUT)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_gcbias"].values()))
    v = d["AT_DROPOUT"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def gc_dropout(mqc):
    """
    Illumina-style GC dropout metric. Calculated by taking each GC bin independently
    and calculating (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[50..100].

    Source: picard GcBiasSummaryMetrics (GC_DROPOUT)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_gcbias"].values()))
    v = d["GC_DROPOUT"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def gc_nc_0_19(mqc):
    """
    Normalized coverage over quintile of GC content ranging from 0 - 19.

    Source:  picard GcBiasSummaryMetrics (GC_NC_0_19)
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_gcbias"].values()))
    v = d["GC_NC_0_19"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def gc_nc_20_39(mqc):
    """
    Analogous to gc_nc_0_19.
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_gcbias"].values()))
    v = d["GC_NC_20_39"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def gc_nc_40_59(mqc):
    """
    Analogous to gc_nc_0_19.
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_gcbias"].values()))
    v = d["GC_NC_40_59"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def gc_nc_60_79(mqc):
    """
    Analogous to gc_nc_0_19.
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_gcbias"].values()))
    v = d["GC_NC_60_79"]
    v = np.round(v, 2)

    r = {k:v}
    return r


def gc_nc_80_100(mqc):
    """
    Analogous to gc_nc_0_19.
    """
    k = inspect.currentframe().f_code.co_name

    d = next(iter(mqc["multiqc_picard_gcbias"].values()))
    v = d["GC_NC_80_100"]
    v = np.round(v, 2)

    r = {k:v}
    return r

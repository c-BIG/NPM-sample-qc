#!/usr/bin/env python3

import numpy as np
import inspect


DECIMALS = 5


def yield_raw_gb(mqc):
    """
    The total number of bases in all reads.

    Source: picard QualityYieldMetrics (TOTAL_BASES)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_picard_QualityYieldMetrics"].values()))
        v = np.round(np.divide(d["TOTAL_BASES"], 1e9), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def yield_pf_gb(mqc):
    """
    The total number of bases in all PF reads.

    Source: picard QualityYieldMetrics (PF_BASES)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_picard_QualityYieldMetrics"].values()))
        v = np.round(np.divide(d["PF_BASES"], 1e9), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_q30_bases(mqc):
    """
    The percentage of PF bases with base quality >= 30.

    Source: picard QualityYieldMetrics (PF_Q30_BASES/PF_BASES)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_picard_QualityYieldMetrics"].values()))
        v = np.round(np.divide(d["PF_Q30_BASES"],
                               d["PF_BASES"])*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_q30_bases_read1(mqc):
    """
    The percentage of PF bases in read 1 with base quality >= 30.

    Source: samtools stats (FFQ) + MultiQC_NPM
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_samtools_stats_bq"].values()))
        v = np.round(d["pct_q30_bases_read1"], DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_q30_bases_read2(mqc):
    """
    Analogous to pct_q30_bases_read1.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_samtools_stats_bq"].values()))
        v = np.round(d["pct_q30_bases_read2"], DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def yield_raw_reads(mqc):
    """
    The total number of reads including all PF and non-PF reads.

    Source: picard AlignmentSummaryMetrics (TOTAL_READS)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
        v = d["TOTAL_READS"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def yield_pf_reads(mqc):
    """
    The number of PF reads where PF is defined as passing Illumina's filter.

    Source: picard AlignmentSummaryMetrics (PF_READS)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
        v = d["PF_READS"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pct_pf_reads(mqc):
    """
    The percentage of reads that are PF, i.e. yield_pf_reads / yield_raw_reads.

    Source: picard AlignmentSummaryMetrics (PCT_PF_READS)

    Note: picard reports a fraction, re-mapped to a percentage here.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
        v = d["PCT_PF_READS"]
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_aligned_bases(mqc):
    """
    The percentage of bases that are aligned to the reference.

    Source: samtools stats (bases mapped cigar / total length)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_stats"].values()))
        v = np.divide(d["bases_mapped_(cigar)"],
                      d["total_length"])
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_reads_aligned(mqc):
    """
    The percentage of PF reads that aligned to the reference.

    Source: picard AlignmentSummaryMetrics (PCT_PF_READS_ALIGNED)

    Note: picard reports a fraction, re-mapped to a percentage here.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
        v = d["PCT_PF_READS_ALIGNED"]
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_reads_aligned_mapqge20(mqc):
    """
    The percentage of PF reads that were aligned to the reference with a mapping quality of Q20 or higher.

    Note: picard reports a fraction, re-mapped to a percentage here.

    Source: picard AlignmentSummaryMetrics (PF_HQ_ALIGNED_READS / TOTAL_READS)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
        v = np.divide(d["PF_HQ_ALIGNED_READS"],
                      d["TOTAL_READS"])
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_reads_unaligned(mqc):
    """
    The percentage of unmapped reads.

    Source: samtools stats (reads_unmapped / sequences)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_stats"].values()))
        v = np.divide(d["reads_unmapped"],
                      d["sequences"])
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_reads_aligned_in_pairs(mqc):
    """
    The percentage of reads that have been aligned as pairs, i.e. percent aligned pairs * 2.

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


def total_alignments(mqc):
    """
    The total number of alignments in the BAM/CRAM file.
    Includes primary, secondary and supplementary alignments; excludes duplicates.

    Source: samtools flagstat (total_passed)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
        v = int(d["total_passed"])
    except KeyError:
        v = "NA"

    return k, v


def pct_singleton_alignments(mqc):
    """
    The percentage of alignments that are singletons.

    Source: samtools flagstat (singletons_passed_pct)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
        v = d["singletons_passed_pct"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_alignments_diff_chrom(mqc):
    """
    The percentage of alignments with mates aligned to different chromosomes.

    Source: samtools flagstat (with mate mapped to a different chr_passed / total_passed)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
        v = np.divide(d["with mate mapped to a different chr_passed"],
                      d["total_passed"])
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_alignments_diff_chrom_mapqge5(mqc):
    """
    The percentage of alignments with MAPQ>=5 and with mates aligned to different chromosomes.

    Source: samtools flagstat (with mate mapped to a different chr (mapQ >= 5)_passed / total_passed)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
        v = np.divide(d["with mate mapped to a different chr (mapQ >= 5)_passed"],
                      d["total_passed"])
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_chimeras(mqc):
    """
    The percentage of reads that map outside of a maximum insert size (usually 100kb)
    or that have the two ends mapping to different chromosomes.

    Source: picard AlignmentSummaryMetrics (PCT_CHIMERAS)

    Note: picard reports a fraction, re-mapped to a percentage here.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
        v = d["PCT_CHIMERAS"]
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_secondary_alignments(mqc):
    """
    The percentage of alignments marked as secondary.

    Source: samtools flagstat (secondary_passed / total_passed)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
        v = np.divide(d["secondary_passed"],
                      d["total_passed"])
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_supplementary_alignments(mqc):
    """
    The percentage of alignments marked as supplementary.

    Source: samtools flagstat (supplementary_passed / total_passed)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_flagstat"].values()))
        v = np.divide(d["supplementary_passed"],
                      d["total_passed"])
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_duplicate_reads(mqc):
    """
    The percentage of reads marked as duplicates.

    Source: samtools stats (reads_duplicated_percent)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_stats"].values()))
        v = d["reads_duplicated_percent"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def mismatch_rate(mqc):
    """
    The fraction of bases mismatching the reference for all bases aligned to the reference sequence.

    Source: picard AlignmentSummaryMetrics (PF_MISMATCH_RATE)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
        v = d["PF_MISMATCH_RATE"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def mismatch_rate_mapqge20(mqc):
    """
    The fraction of bases that mismatch the reference in reads with MAPQ>=20.

    Source: picard AlignmentSummaryMetrics (PF_HQ_ERROR_RATE)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
        v = d["PF_HQ_ERROR_RATE"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def genome_territory(mqc):
    """
    The number of non-N bases in autosomes over which coverage will be evaluated.

    Source: run_mosdepth.sh (total_autosome_bases)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = d["total_autosome_bases"]
        v = int(np.float(v))
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
        v = np.round(np.float(v), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def sd_autosome_coverage(mqc):
    """
    The standard deviation of coverage in autosomes, after coverage filters are applied (see mean_autosome_coverage for details).

    Source: run_mosdepth.sh (sd_autosome_coverage)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = d["sd_autosome_coverage"]
        v = np.round(np.float(v), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def median_autosome_coverage(mqc):
    """
    The median coverage in autosomes, after coverage filters are applied (see mean_autosome_coverage for details).

    Source: run_mosdepth.sh (median_autosome_coverage)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = d["median_autosome_coverage"]
        v = np.round(np.float(v), DECIMALS)
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
        v = np.round(np.float(v), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_autosomes_1x(mqc):
    """
    The percentage of bases that attained at least 1X sequence coverage in autosomes, after coverage filters are applied (see mean_autosome_coverage for details).

    Source: run_mosdepth.sh (ge_1x_autosome_bases/total_autosome_bases*100)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = np.divide(np.float(d["ge_1x_autosome_bases"]),
                      np.float(d["total_autosome_bases"]))
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_autosomes_10x(mqc):
    """
    Analogous to pct_autosomes_1x.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = np.divide(np.float(d["ge_10x_autosome_bases"]),
                      np.float(d["total_autosome_bases"]))
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_autosomes_15x(mqc):
    """
    Analogous to pct_autosomes_1x.
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


def pct_autosomes_30x(mqc):
    """
    Analogous to pct_autosomes_1x.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = np.divide(np.float(d["ge_30x_autosome_bases"]),
                      np.float(d["total_autosome_bases"]))
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_autosomes_40x(mqc):
    """
    Analogous to pct_autosomes_1x.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_mosdepth"].values()))
        v = np.divide(np.float(d["ge_40x_autosome_bases"]),
                      np.float(d["total_autosome_bases"]))
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def coverage_sg10k_062017(mqc):
    """
    The mean coverage in bases of the entire genome (3.1e9), after SG10K filters are applied.

    Note: filters include the following

    - before mapping: remove reads with >50% of bases with base quality <= 10
    - after mapping: remove duplicate reads, clipped bases and bases with base quality < 5

    Source: sg10k-cov-062017.sh (bases_sg10k_062017)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_sg10k_cov_062017"].values()))
        covered_bases = d["bases_sg10k_062017"]
        total_bases = float(3.1e9)
        v = np.round(np.divide(covered_bases, total_bases), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def mean_insert_size(mqc):
    """
    The mean insert size over the "core" of the distribution.

    Note: Artefactual outliers in the distribution often cause calculation of nonsensical mean and
    stdev values. To avoid this the distribution is first trimmed to a "core" distribution of +/- N
    median absolute deviations around the median insert size.

    Source: picard InsertSizeMetrics (MEAN_INSERT_SIZE)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_insertSize"].values()))
        v = d["MEAN_INSERT_SIZE"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def sd_insert_size(mqc):
    """
    The standard deviation of insert sizes over the "core" of the distribution.

    Source: picard InsertSizeMetrics (STANDARD_DEVIATION)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_insertSize"].values()))
        v = d["STANDARD_DEVIATION"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def median_insert_size(mqc):
    """
    The median insert size over the "core" of the distribution.

    Source: picard InsertSizeMetrics (MEDIAN_INSERT_SIZE)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_insertSize"].values()))
        v = d["MEDIAN_INSERT_SIZE"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def mad_insert_size(mqc):
    """
    The median absolute deviation of insert sizes over the "core" of the distribution.

    Source: picard InsertSizeMetrics (MEDIAN_ABSOLUTE_DEVIATION)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_insertSize"].values()))
        v = d["MEDIAN_ABSOLUTE_DEVIATION"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_overlapping_bases(mqc):
    """
    The percentage of bases that correspond to the second observation from an insert with overlapping reads.

    Source: picard WgsMetrics (PCT_EXC_OVERLAP)

    Note: picard reports a fraction, re-mapped to a percentage here.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_wgsmetrics"].values()))
        v = d["PCT_EXC_OVERLAP"]
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_adapters(mqc):
    """
    The percentage of PF reads that are unaligned and match to a known adapter sequence right
    from the start of the read.

    Source: picard AlignmentSummaryMetrics (PCT_ADAPTER)

    Note: picard reports a fraction, re-mapped to a percentage here.

    See picard's source code for details on the adapter sequences
    considered: https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/util/IlluminaUtil.java#L130
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_AlignmentSummaryMetrics"].values()))
        v = d["PCT_ADAPTER"]
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def at_dropout(mqc):
    """
    Illumina-style AT dropout metric. Calculated by taking each GC bin independently
    and calculating (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[0..50].

    Source: picard GcBiasSummaryMetrics (AT_DROPOUT)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_gcbias"].values()))
        v = d["AT_DROPOUT"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def gc_dropout(mqc):
    """
    Illumina-style GC dropout metric. Calculated by taking each GC bin independently
    and calculating (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[50..100].

    Source: picard GcBiasSummaryMetrics (GC_DROPOUT)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_gcbias"].values()))
        v = d["GC_DROPOUT"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def gc_nc_0_19(mqc):
    """
    Normalized coverage over quintile of GC content ranging from 0 - 19.

    Source:  picard GcBiasSummaryMetrics (GC_NC_0_19)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_gcbias"].values()))
        v = d["GC_NC_0_19"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def gc_nc_20_39(mqc):
    """
    Analogous to gc_nc_0_19.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_gcbias"].values()))
        v = d["GC_NC_20_39"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def gc_nc_40_59(mqc):
    """
    Analogous to gc_nc_0_19.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_gcbias"].values()))
        v = d["GC_NC_40_59"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def gc_nc_60_79(mqc):
    """
    Analogous to gc_nc_0_19.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_gcbias"].values()))
        v = d["GC_NC_60_79"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def gc_nc_80_100(mqc):
    """
    Analogous to gc_nc_0_19.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_picard_gcbias"].values()))
        v = d["GC_NC_80_100"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def all_snps(mqc):
    """
    The total number of SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_snps"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_het_snps(mqc):
    """
    The total number of het SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_het_snps"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_homalt_snps(mqc):
    """
    The total number of hom alt SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_homalt_snps"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_snp_het_hom(mqc):
    """
    The het/hom ratio for all SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_snp_het_hom"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pass_snps(mqc):
    """
    The number of PASS SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_snps"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_het_snps(mqc):
    """
    The number of PASS het SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_het_snps"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_homalt_snps(mqc):
    """
    The number of PASS hom alt SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_homalt_snps"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_snp_het_hom(mqc):
    """
    The het/hom ratio for PASS SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_snp_het_hom"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def all_indels(mqc):
    """
    The total number of INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_indels"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_het_indels(mqc):
    """
    The total number of het INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_het_indels"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_homalt_indels(mqc):
    """
    The total number of homalt INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_homalt_indels"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_indel_het_hom(mqc):
    """
    The het/hom ratio for all INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_indel_het_hom"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pass_indels(mqc):
    """
    The number of PASS INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_indels"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_het_indels(mqc):
    """
    The number of PASS het INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_het_indels"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_homalt_indels(mqc):
    """
    The number of PASS homalt INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_homalt_indels"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_indel_het_hom(mqc):
    """
    The het/hom ratio for PASS INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_indel_het_hom"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v



def all_del(mqc):
    """
    The total number of deletions.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_del"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_ins(mqc):
    """
    The total number of insertions.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_ins"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_del(mqc):
    """
    The number of PASS deletions.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_del"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_ins(mqc):
    """
    The number of PASS insertions.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_ins"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_mnps(mqc):
    """
    The total number of MNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_mnps"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_mnps(mqc):
    """
    The number of PASS MNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_mnps"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_complex_indels(mqc):
    """
    The total number of complex INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_complex_indels"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_complex_ins(mqc):
    """
    The total number of complex insertions.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_complex_ins"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_complex_del(mqc):
    """
    The total number of complex deletions.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_complex_del"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_complex_indels(mqc):
    """
    The number of PASS complex INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_complex_indels"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_complex_ins(mqc):
    """
    The number of PASS complex insertions.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_complex_ins"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_complex_del(mqc):
    """
    The number of PASS complex deletions.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_complex_del"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_multiallelic_sites(mqc):
    """
    The total number of multiallelic sites.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_multiallelic_sites"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def pass_multiallelic_sites(mqc):
    """
    The number of PASS multiallelic sites.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_multiallelic_sites"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v


def all_snp_ts_tv(mqc):
    """
    The transition to transversion ratio of all bi-allelic SNPs.

    Source: count_variants.py (bcftools stats - tstv)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["all_snp_ts_tv"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pass_snp_ts_tv(mqc):
    """
    The transition to transversion ratio of passing bi-allelic SNPs.

    Source: count_variants.py (bcftools stats - tstv)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_count_variants"].values()))
        v = d["pass_snp_ts_tv"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pct_contamination(mqc):
    """
    The percentage of cross-individual contamination.

    Source: VerifyBamID2 (FREEMIX)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_verifybamid"].values()))
        v = d["FREEMIX"]
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def pst_pct_concordance(mqc):
    """
    The percentage of sites with concordant genotypes between a query VCF and a PST VCF (positive sample tracking).
    Includes SNPs and INDELs.

    Source: bcftools gtcheck + MultiQC_NPM (pst_pct_concordance)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_npm_bcftools_gtcheck"].values()))
        v = d["pst_pct_concordance"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

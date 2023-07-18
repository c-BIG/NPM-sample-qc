#!/usr/bin/env python3

import numpy as np
import inspect

DECIMALS = 5

def yield_bp_q30(mqc):
    """
    Description: The number of bases in short paired-end sequencing high quality reads, primary alignments, achieving a base quality score of 30 or greater (Phred scale). Duplicated reads and clipped bases are included. No minimum mapping quality is imposed.

    Implementation details: In the NPM-sample-QC reference implementation it is computed using GATK Picard CollectQualityYieldMetrics, reporting the PF_Q30_BASES field. Only high quality bases from primary alignments are considered. No filter on duplicated reads, clipped bases or mapping qualiy is applied.
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
    Description: The percentage of short paired-end sequencing high quality reads, primary alignments, mapped on GRCh38 assembly. Duplicated reads and clipped bases are included. No minimum mapping quality is imposed.
    
    Implementation details: In the NPM-sample-QC reference implementation it is computed using samtools stats, reporting the percentage of reads mapped on GRCh38 assembly. Duplicated reads are included. No mapping qualiy is applied.
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
    Description: The percentage of short paired-end sequencing high quality, properly paired reads, primary alignments, mapped on GRCh38 assembly. Duplicated reads are included. No minimum mapping quality is imposed.
    
    Implementation details: In the NPM-sample-QC reference implementation it is computed using samtools stats, reporting the percentage of properly paired reads mapped on GRCh38 assembly. Duplicated reads are included. No mapping qualiy is applied.
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
    Description: The mean sequencing coverage derived from short paired-end sequencing high quality, non duplicated reads, primary alignments, achieving a mapping quality of 20 or greater, in autosomes non gap regions of GRCh38 assembly. Clipped bases are excluded. Overlapping bases are counted only once. It is critical that the (BAM/CRAM) alignment files be readily marked for duplicated reads and clipped bases.
    
    Implementation details: In the NPM-sample-QC reference implementation, the genome-wide sequencing coverage of non duplicated reads, non clipped bases, non overlapping bases, primary alignments, achieving a mapping quality of 20 or greater is derived from mosdepth v0.3.2. It is further narrowed down to the non gap regions of GRCh38 assembly, autosomes only using bedtools intersect. The mean coverage is then computed on 1,000bp windows and averaged for the selected region using datamash.
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
    Description: The median absolute deviation of sequencing coverage derived from short paired-end sequencing high quality, non duplicated reads, primary alignments, achieving a mapping quality of 20 or greater, in autosomes non gap regions of GRCh38 assembly. Clipped bases are excluded. Overlapping bases are counted only once. It is critical that the (BAM/CRAM) alignment files be readily marked for duplicated reads and clipped bases.
    
    Implementation details: In the NPM-sample-QC reference implementation, the genome-wide sequencing coverage of non duplicated reads, non clipped bases, non overlapping bases, primary alignments, achieving a mapping quality of 20 or greater is derived from mosdepth v0.3.2. It is further narrowed down to the non gap regions of GRCh38 assembly, autosomes only using bedtools intersect. The median absolute deviation of the coverage is then calculated using datamash.
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
    Description: The percentage of bases attaining at least 15X sequencing coverage in short paired-end sequencing high quality, non duplicated reads, primary alignments, achieving a mapping quality of 20 or greater, in autosomes non gap regions of GRCh38 assembly. Clipped bases are excluded. Overlapping bases are counted only once. It is critical that the (BAM/CRAM) alignment files be readily marked for duplicated reads and clipped bases.
    
    Implementation details: In the NPM-sample-QC reference implementation, the genome-wide sequencing coverage of non duplicated reads, non clipped bases, non overlapping bases, primary alignments, achieving a mapping quality of 20 or greater is derived from mosdepth v0.3.2. It is further narrowed down to the non gap regions of GRCh38 assembly, autosomes only using bedtools intersect. The percentage of bases attaining at least 15X coverage is then calculated using datamash.
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
    Description: The mean insert size of short paired-end sequencing high quality reads, primary alignments, mapped on GRCh38 assembly. Duplicated reads and clipped bases are included. No minimum mapping quality is imposed.
    
    Implementation details: In the NPM-sample-QC reference implementation it is computed using samtools stats, reporting the insert_size_average field. Duplicated reads are included. No mapping qualiy is applied.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_stats"].values()))
        v = d["insert_size_average"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

def insert_size_std_deviation(mqc):
    """
    Description: The insert size standard deviation of short paired-end sequencing high quality reads, primary alignments, mapped on GRCh38 assembly. Duplicated reads and clipped bases are included. No minimum mapping quality is imposed.
    
    Implementation details: In the NPM-sample-QC reference implementation it is computed using samtools stats, reporting the insert_size_standard_deviation field. Duplicated reads are included. No mapping qualiy is applied.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_samtools_stats"].values()))
        v = d["insert_size_standard_deviation"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

def cross_contamination_rate(mqc):
    """
    Description: Estimation of inter-sample contamination rate of short paired-end sequencing high quality, non duplicated reads, primary alignments, mapped on GRCh38 assembly. No minimum mapping quality is imposed. It is critical that the (BAM/CRAM) alignment files be readily marked for duplicated reads and clipped bases.
    
    Implementation details: The estimation of inter-sample DNA contamination of short paired-end sequencing high quality, aligned sequence reads (BAM/CRAM) mapped on GRCh38 assembly with pre-calculated reference panel of 1000 Genome Project dataset from the VerifyBamID resource using VerifyBamID2 with NumPC “4” (# of Principal Components used in estimation), the key information “FREEMIX” in “.selfSM” in the results indicates the estimated contamination level.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = next(iter(mqc["multiqc_verifybamid"].values()))
        v = d["FREEMIX"]
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

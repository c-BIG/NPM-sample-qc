#!/usr/bin/env python3

import numpy as np
import inspect

DECIMALS = 5

def yield_bp_q30(mqc, biosample_id):
    """
    Description: The number of bases in short paired-end sequencing high quality reads, primary alignments, achieving a base quality score of 30 or greater (Phred scale). Duplicated reads and clipped bases are included. No minimum mapping quality is imposed.

    Implementation details: In the NPM-sample-QC reference implementation it is computed using GATK Picard CollectQualityYieldMetrics, reporting the PF_Q30_BASES field. Only high quality bases from primary alignments are considered. No filter on duplicated reads, clipped bases or mapping qualiy is applied.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_npm_picard_QualityYieldMetrics"][biosample_id]
        v = d["PF_Q30_BASES"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v

def pct_reads_mapped(mqc, biosample_id):
    """
    Description: The percentage of short paired-end sequencing high quality reads, primary alignments, mapped on GRCh38 assembly. Duplicated reads and clipped bases are included. No minimum mapping quality is imposed.
    
    Implementation details: In the NPM-sample-QC reference implementation it is computed using samtools stats, reporting the percentage of reads mapped on GRCh38 assembly. Duplicated reads are included. No mapping qualiy is applied.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_samtools_stats"][biosample_id]
        v = np.divide(d["reads_mapped"],
                      d["sequences"])
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

def pct_reads_properly_paired(mqc, biosample_id):
    """
    Description: The percentage of short paired-end sequencing high quality, properly paired reads, primary alignments, mapped on GRCh38 assembly. Duplicated reads are included. No minimum mapping quality is imposed.
    
    Implementation details: In the NPM-sample-QC reference implementation it is computed using samtools stats, reporting the percentage of properly paired reads mapped on GRCh38 assembly. Duplicated reads are included. No mapping qualiy is applied.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_samtools_stats"][biosample_id]
        v = d["reads_properly_paired_percent"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

def mean_autosome_coverage(mqc, biosample_id):
    """
    Description: The mean sequencing coverage derived from short paired-end sequencing high quality, non duplicated reads, primary alignments, achieving a base quality of 20 or greater and mapping quality of 20 or greater, in autosomes non gap regions of GRCh38 assembly. Overlapping bases are counted only once. It is critical that the (BAM/CRAM) alignment files be readily marked for duplicated reads.
    
    Implementation details: In the NPM-sample-QC reference implementation, the genome-wide sequencing mean coverage of the non gap regions of GRCh38 assembly, autosomes only, non duplicated reads, non overlapping bases, primary alignments, achieving a base quality of 20 or greater and mapping quality of 20 or greater is derived from picard (2.27.0) CollectWgsMetrics.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_picard_wgsmetrics"][biosample_id]
        v = d["MEAN_COVERAGE"]
        v = np.round(np.float64(v), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

def mad_autosome_coverage(mqc, biosample_id):
    """
    Description: The median absolute deviation of sequencing coverage derived from short paired-end sequencing high quality, non duplicated reads, primary alignments, achieving a base quality of 20 or greater and mapping quality of 20 or greater, in autosomes non gap regions of GRCh38 assembly. Overlapping bases are counted only once. It is critical that the (BAM/CRAM) alignment files be readily marked for duplicated reads.
    
    Implementation details: In the NPM-sample-QC reference implementation, the genome-wide sequencing median absolute deviation coverage of the non gap regions of GRCh38 assembly, autosomes only, non duplicated reads, non overlapping bases, primary alignments, achieving a base quality of 20 or greater and mapping quality of 20 or greater is derived from picard (2.27.0) CollectWgsMetrics.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_picard_wgsmetrics"][biosample_id]
        v = d["MAD_COVERAGE"]
        v = np.round(np.float64(v), DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

def pct_autosomes_15x(mqc, biosample_id):
    """
    Description: The percentage of bases attaining at least 15X sequencing coverage in short paired-end sequencing high quality, non duplicated reads, primary alignments, achieving a base quality of 20 or greater and mapping quality of 20 or greater, in autosomes non gap regions of GRCh38 assembly. Overlapping bases are counted only once. It is critical that the (BAM/CRAM) alignment files be readily marked for duplicated reads.
    
    Implementation details: In the NPM-sample-QC reference implementation, the genome-wide sequencing coverage percentage of bases attaining at least 15X of the non gap regions of GRCh38 assembly, autosomes only, non duplicated reads, non overlapping bases, primary alignments, achieving a base quality of 20 or greater and mapping quality of 20 or greater is derived from picard (2.27.0) CollectWgsMetrics.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_picard_wgsmetrics"][biosample_id]
        v = d["PCT_15X"]
        v = np.round(v*100, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

def mean_insert_size(mqc, biosample_id):
    """
    Description: The mean insert size of short paired-end sequencing high quality reads, primary alignments, mapped on GRCh38 assembly. Duplicated reads and clipped bases are included. No minimum mapping quality is imposed.
    
    Implementation details: In the NPM-sample-QC reference implementation it is computed using samtools stats, reporting the insert_size_average field. Duplicated reads are included. No mapping qualiy is applied.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_samtools_stats"][biosample_id]
        v = d["insert_size_average"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

def insert_size_std_deviation(mqc, biosample_id):
    """
    Description: The insert size standard deviation of short paired-end sequencing high quality reads, primary alignments, mapped on GRCh38 assembly. Duplicated reads and clipped bases are included. No minimum mapping quality is imposed.
    
    Implementation details: In the NPM-sample-QC reference implementation it is computed using samtools stats, reporting the insert_size_standard_deviation field. Duplicated reads are included. No mapping qualiy is applied.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_samtools_stats"][biosample_id]
        v = d["insert_size_standard_deviation"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

def cross_contamination_rate(mqc, biosample_id):
    """
    Description: Estimation of inter-sample contamination rate of short paired-end sequencing high quality, non duplicated reads, primary alignments, mapped on GRCh38 assembly. No minimum mapping quality is imposed. It is critical that the (BAM/CRAM) alignment files be readily marked for duplicated reads and clipped bases.
    
    Implementation details: The estimation of inter-sample DNA contamination of short paired-end sequencing high quality, aligned sequence reads (BAM/CRAM) mapped on GRCh38 assembly with pre-calculated reference panel of 1000 Genome Project dataset from the VerifyBamID resource using VerifyBamID2 with NumPC “4” (# of Principal Components used in estimation), the key information “FREEMIX” in “.selfSM” in the results indicates the estimated contamination level.
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_verifybamid"][biosample_id]
        v = d["FREEMIX"]
    except KeyError:
        v = "NA"

    return k, v

def count_snvs(mqc, biosample_id):
    """
    The number of PASS SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_npm_count_variants"][biosample_id + ".variant"]
        v = d["pass_snps"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v

def pass_het_snps(mqc, biosample_id):
    """
    The number of PASS het SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_npm_count_variants"][biosample_id + ".variant"]
        v = d["pass_het_snps"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v

def pass_homalt_snps(mqc, biosample_id):
    """
    The number of PASS hom alt SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_npm_count_variants"][biosample_id + ".variant"]
        v = d["pass_homalt_snps"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v

def ratio_heterozygous_homzygous_snv(mqc, biosample_id):
    """
    The het/hom ratio for PASS SNPs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_npm_count_variants"][biosample_id + ".variant"]
        v = d["pass_snp_het_hom"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

def ratio_heterozygous_homzygous_indel(mqc, biosample_id):
    """
    The het/hom ratio for PASS INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_npm_count_variants"][biosample_id + ".variant"]
        v = d["pass_indel_het_hom"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

def count_deletions(mqc, biosample_id):
    """
    The number of PASS deletions.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_npm_count_variants"][biosample_id + ".variant"]
        v = d["pass_del"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v

def count_insertions(mqc, biosample_id):
    """
    The number of PASS insertions.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_npm_count_variants"][biosample_id + ".variant"]
        v = d["pass_ins"]
        v = int(v)
    except KeyError:
        v = "NA"

    return k, v

def ratio_insertion_deletion(mqc, biosample_id):
    """
    The ins/del ratio for PASS INDELs.

    Source: count_variants.py (bcftools view)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_npm_count_variants"][biosample_id + ".variant"]
        v = d["pass_ins_del"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v


def ratio_transitions_transversions(mqc, biosample_id):
    """
    The transition to transversion ratio of passing bi-allelic SNPs.

    Source: count_variants.py (bcftools stats - tstv)
    """
    k = inspect.currentframe().f_code.co_name

    try:
        d = mqc["multiqc_npm_count_variants"][biosample_id + ".variant"]
        v = d["pass_snp_ts_tv"]
        v = np.round(v, DECIMALS)
    except KeyError:
        v = "NA"

    return k, v

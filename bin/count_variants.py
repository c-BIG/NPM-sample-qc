#!/usr/bin/env python3

import logging
import argparse
import json
import pprint as pp
import subprocess
import numpy as np
import sys
import os
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_vcf", dest="input_vcf", required=True,
                        default=None,
                        help="Path to input VCF/gVCF.")
    parser.add_argument("--regions", dest="regions", required=True,
                        default=None,
                        help="Path to input autosomes_non_gap_regions_bed.")
    parser.add_argument("--output_json", dest="output_json", required=False,
                        default="./variant_counts.json",
                        help="Path to output file for variant metrics. Default: ./variant_counts.json")
    parser.add_argument("--biosample_id", dest="biosample_id", required=True,
                        default=None,
                        help="Sample ID.")
    parser.add_argument("--scratch_dir", dest="scratch_dir", required=False,
                        default="./",
                        help="Path to scratch dir. Default: ./")
    parser.add_argument("--keep_scratch", dest="keep_scratch", required=False,
                        default=False, action="store_true",
                        help="Keep scratch dir. Default: False.")
    parser.add_argument("--loglevel", dest="loglevel", required=False,
                        default="INFO",
                        help="Set logging level to INFO (default), WARNING or DEBUG.")
    args = parser.parse_args()

    # create scratch dir if it doesn't exist
    Path(args.scratch_dir).mkdir(parents=True, exist_ok=True)

    return args


def set_logging(loglevel):
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        level=numeric_level)


def setup(input_vcf, scratch_dir):

    # subset input vcf to include only variant sites to speed up execution from gvcf
    logging.info("Subsetting variant sites from input vcf...")
    variants_vcf = "%s/variants.vcf.gz" % scratch_dir
    cmd = "bcftools view"
    cmd += " -e 'GT=\"./.\" | GT=\".|.\" | GT=\".\" | GT=\"0/0\" | GT=\"0|0\" | GT=\"0\"'"
    cmd += " -O z -o %s %s" % (variants_vcf, input_vcf)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.run(cmd, shell=True)

    return variants_vcf


def count_variants(input_vcf, scratch_dir, regions):

    # initialise results dict
    r = dict()
    metrics_list = [
        "count_snvs", "ratio_heterozygous_homzygous_snv",
        "ratio_heterozygous_homzygous_indel",
        "count_deletions", "count_insertions", "ratio_insertion_deletion",
        "ratio_transitions_transversions_snv"
    ]
    for m in metrics_list:
        r[m] = 0

    # run bcftools stats - used to obtain ts/tv and sanity checks
    logging.info("Running bcftools stats for all variants...")
    cmd = "tabix -p vcf %s" % input_vcf
    cmd1 = "bcftools stats -R %s %s > %s/all.bcftools_stats" % (regions, input_vcf, scratch_dir)
    logging.debug("CMD: %s" % cmd)
    logging.debug("CMD: %s" % cmd1)
    p = subprocess.run(cmd, shell=True)
    p = subprocess.run(cmd1, shell=True)

    logging.info("Running bcftools stats for PASS variants...")
    cmd = "bcftools stats -f PASS -R %s %s > %s/pass.bcftools_stats" % (regions, input_vcf, scratch_dir)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.run(cmd, shell=True)

    # calculate metrics
    logging.info("Counting all_snps...")
    cmd = "bcftools view -H -v snps -R %s %s | wc -l" % (regions, input_vcf)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    # r["all_snps"] = int(p.stdout.read())
    all_snps = int(p.stdout.read())

    logging.info("Counting count_snvs...")
    cmd = "bcftools view -H -v snps -f PASS -R %s %s | wc -l" % (regions, input_vcf)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["count_snvs"] = int(p.stdout.read())

    logging.info("Counting pass_het_snps...")
    cmd = "bcftools view -H -v snps -f PASS -g het -R %s %s | wc -l" % (regions, input_vcf)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    # r["pass_het_snps"] = int(p.stdout.read())
    pass_het_snps = int(p.stdout.read())

    logging.info("Counting pass_homalt_snps...")
    cmd = "bcftools view -H -v snps -f PASS -g hom -R %s %s | wc -l" % (regions, input_vcf)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    # r["pass_homalt_snps"] = int(p.stdout.read())
    pass_homalt_snps = int(p.stdout.read())


    logging.info("Calculating ratio_heterozygous_homzygous_snv...")
    # snp_het_hom = np.divide(r["pass_het_snps"], r["pass_homalt_snps"])
    snp_het_hom = np.divide(pass_het_snps, pass_homalt_snps)
    r["ratio_heterozygous_homzygous_snv"] = np.round(snp_het_hom, 2)

    logging.info("Counting all_indels...")
    cmd = "bcftools view -H -v indels -R %s %s | wc -l" % (regions, input_vcf)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    # r["all_indels"] = int(p.stdout.read())
    all_indels = int(p.stdout.read())

    logging.info("Counting pass_indels...")
    cmd = "bcftools view -H -v indels -f PASS -R %s %s | wc -l" % (regions, input_vcf)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    # r["pass_indels"] = int(p.stdout.read())
    pass_indels = int(p.stdout.read())

    logging.info("Counting pass_het_indels...")
    cmd = "bcftools view -H -v indels -f PASS -g het -R %s %s | wc -l" % (regions, input_vcf)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    # r["pass_het_indels"] = int(p.stdout.read())
    pass_het_indels = int(p.stdout.read())

    logging.info("Counting pass_homalt_indels...")
    cmd = "bcftools view -H -v indels -f PASS -g hom -R %s %s | wc -l" % (regions, input_vcf)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    # r["pass_homalt_indels"] = int(p.stdout.read())
    pass_homalt_indels = int(p.stdout.read())

    logging.info("Calculating ratio_heterozygous_homzygous_indel...")
    # indel_het_hom = np.divide(r["pass_het_indels"], r["pass_homalt_indels"])
    indel_het_hom = np.divide(pass_het_indels, pass_homalt_indels)
    r["ratio_heterozygous_homzygous_indel"] = np.round(indel_het_hom, 2)


    logging.info("Counting count_deletions, count_insertions...")
    cmd = "bcftools view -H -v indels -f PASS -R %s %s" % (regions, input_vcf)
    # filter dragen annotation of variant records
    cmd += " | sed 's/,<NON_REF>//g'"
    # exclude multi-allelic sites
    cmd += " | awk '$4 !~ /,/ && $5 !~ /,/'"
    # if length(ALT) > length(REF), label as INS
    #cmd += " | awk '{ if( length($5) > length($4) ) { print \"INS\" }"
    cmd += " | awk '{ if( length($5) > length($4) ) { print \"insertions\" }"
    # if length(ALT) < length(REF), label as DEL
    #cmd += " else if ( length($5) < length($4) ) { print \"DEL\" }"
    cmd += " else if ( length($5) < length($4) ) { print \"deletions\" }"
    # label everything else as UNK and report its location
    cmd += " else { print \"UNK\" } }'"
    # count by label
    cmd += " | sort | uniq -c"
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    x = p.stdout.read()
    for l in x.splitlines():
        l = l.decode("utf-8").strip().split(" ")
        k = "count_%s" % l[1].lower()
        v = int(l[0])
        r[k] = v

    logging.info("Calculating ratio_insertion_deletion...")
    ins_del = np.divide(r["count_insertions"], r["count_deletions"])
    r["ratio_insertion_deletion"] = np.round(ins_del, 2)

    logging.info("Calculating pass_snp_ts_tv...")
    cmd = "cat %s/pass.bcftools_stats | grep '^TSTV'" % scratch_dir
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    l = p.stdout.read().decode("utf-8").strip().split("\t")
    r["ratio_transitions_transversions_snv"] = np.round(float(l[4]), 2)

    # sanity checks
    logging.info("Running sanity checks for variant counts...")

    error_count = 0
    cmd = "cat %s/all.bcftools_stats | grep '^SN'" % scratch_dir
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

    for l in p.stdout.read().splitlines():
        l = l.decode("utf-8").strip().split("\t")
        k = str(l[2]).lower()
        v = int(l[3])
        if "snps" in k and v != all_snps:
        #if "snps" in k and v != r["all_snps"]:
            # logging.error("Mismatch in SNP counts: bcftools stats=%d, count_variants=%d" % (v, r["all_snps"]))
            logging.error("Mismatch in SNP counts: bcftools stats=%d, count_variants=%d" % (v, all_snps))
            error_count += 1
        # elif "indels" in k and v != r["all_indels"]:    
        elif "indels" in k and v != all_indels:
            # logging.error("Mismatch in INDEL counts: bcftools stats=%d, count_variants=%d" % (v, r["all_indels"]))
            logging.error("Mismatch in INDEL counts: bcftools stats=%d, count_variants=%d" % (v, all_indels))
            error_count += 1

    if error_count > 0:
        sys.exit("Sanity checks failed, aborting execution")

    return r


def save_output(data_metrics, outfile):
    with open(outfile, "w") as f:
        # json.dump(data_metrics, f, sort_keys=True, indent=4)
        #data_metrics = {"biosample" : {"id" : args.biosample_id}, "wgs_qc_metrics" : data_metrics}
        data_metrics = {"biosample" : {"id" : args.biosample_id}, "wgs_qc_metrics" : {"variant_metrics" : data_metrics}}
        json.dump(data_metrics, f, sort_keys=True, indent=4)
        f.write("\n")


def done(scratch_dir):
    if not args.keep_scratch:
        # clean up
        logging.info("Cleaning up...")
        os.remove("%s/variants.vcf.gz" % scratch_dir)
        os.remove("%s/variants.vcf.gz.tbi" % scratch_dir)
        os.remove("%s/all.bcftools_stats" % scratch_dir)
        os.remove("%s/pass.bcftools_stats" % scratch_dir)
        if len(os.listdir(scratch_dir)) == 0:
            os.rmdir(scratch_dir)
        else:
            logging.warning("Couldn't remove scratch dir (%s): directory is not empty" % scratch_dir)
    logging.info("DONE")


if __name__ == "__main__":
    args = parse_args()
    set_logging(args.loglevel)

    variants_vcf = setup(args.input_vcf, args.scratch_dir)

    vcf_metrics = count_variants(variants_vcf, args.scratch_dir, args.regions)
    save_output(vcf_metrics, args.output_json)

    done(args.scratch_dir)

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


def count_variants(input_vcf, scratch_dir):

    # initialise results dict
    r = dict()
    metrics_list = [
        "all_snps",
        "pass_snps", "pass_het_snps", "pass_homalt_snps", "pass_snp_het_hom",
        "all_indels",
        "pass_indels", "pass_het_indels", "pass_homalt_indels", "pass_indel_het_hom",
        "pass_del", "pass_ins", "pass_ins_del",
        "all_mnps", "pass_mnps",
        "pass_complex_indels", "pass_complex_ins", "pass_complex_del",
        "pass_multiallelic_sites",
        "pass_snp_ts_tv"
    ]
    for m in metrics_list:
        r[m] = 0

    # run bcftools stats - used to obtain ts/tv and sanity checks
    logging.info("Running bcftools stats for all variants...")
    cmd = "bcftools stats %s > %s/all.bcftools_stats" % (input_vcf, scratch_dir)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.run(cmd, shell=True)

    logging.info("Running bcftools stats for PASS variants...")
    cmd = "bcftools stats -f PASS %s > %s/pass.bcftools_stats" % (input_vcf, scratch_dir)
    logging.debug("CMD: %s" % cmd)
    p = subprocess.run(cmd, shell=True)

    # calculate metrics
    logging.info("Counting all_snps...")
    cmd = "bcftools view -H -v snps %s | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["all_snps"] = int(p.stdout.read())

    logging.info("Counting pass_snps...")
    cmd = "bcftools view -H -v snps -f PASS %s | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["pass_snps"] = int(p.stdout.read())

    logging.info("Counting pass_het_snps...")
    cmd = "bcftools view -H -v snps -f PASS -g het %s | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["pass_het_snps"] = int(p.stdout.read())

    logging.info("Counting pass_homalt_snps...")
    cmd = "bcftools view -H -v snps -f PASS -g hom %s | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["pass_homalt_snps"] = int(p.stdout.read())

    logging.info("Calculating pass_snp_het_hom...")
    snp_het_hom = np.divide(r["pass_het_snps"], r["pass_homalt_snps"])
    r["pass_snp_het_hom"] = np.round(snp_het_hom, 2)

    logging.info("Counting all_indels...")
    cmd = "bcftools view -H -v indels %s | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["all_indels"] = int(p.stdout.read())

    logging.info("Counting pass_indels...")
    cmd = "bcftools view -H -v indels -f PASS %s | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["pass_indels"] = int(p.stdout.read())

    logging.info("Counting pass_het_indels...")
    cmd = "bcftools view -H -v indels -f PASS -g het %s | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["pass_het_indels"] = int(p.stdout.read())

    logging.info("Counting pass_homalt_indels...")
    cmd = "bcftools view -H -v indels -f PASS -g hom %s | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["pass_homalt_indels"] = int(p.stdout.read())

    logging.info("Calculating pass_indel_het_hom...")
    indel_het_hom = np.divide(r["pass_het_indels"], r["pass_homalt_indels"])
    r["pass_indel_het_hom"] = np.round(indel_het_hom, 2)


    logging.info("Counting pass_del, pass_ins...")
    cmd = "bcftools view -H -v indels -f PASS %s" % input_vcf
    # filter dragen annotation of variant records
    cmd += " | sed 's/,<NON_REF>//g'"
    # exclude multi-allelic sites
    cmd += " | awk '$4 !~ /,/ && $5 !~ /,/'"
    # if length(ALT) > length(REF), label as INS
    cmd += " | awk '{ if( length($5) > length($4) ) { print \"INS\" }"
    # if length(ALT) < length(REF), label as DEL
    cmd += " else if ( length($5) < length($4) ) { print \"DEL\" }"
    # label everything else as UNK and report its location
    cmd += " else { print \"UNK\" } }'"
    # count by label
    cmd += " | sort | uniq -c"
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    x = p.stdout.read()
    for l in x.splitlines():
        l = l.decode("utf-8").strip().split(" ")
        k = "pass_%s" % l[1].lower()
        v = int(l[0])
        r[k] = v

    logging.info("Calculating pass_ins_del...")
    ins_del = np.divide(r["pass_ins"], r["pass_del"])
    r["pass_ins_del"] = np.round(ins_del, 2)

    logging.info("Counting all_mnps...")
    cmd = "bcftools view -H -v mnps %s | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["all_mnps"] = int(p.stdout.read())

    logging.info("Counting pass_mnps...")
    cmd = "bcftools view -H -v mnps -f PASS %s | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["pass_mnps"] = int(p.stdout.read())


    logging.info("Counting pass_complex_indels...")
    cmd = "bcftools view -H -v indels -f PASS %s" % input_vcf
    # filter dragen annotation of variant records
    cmd += " | sed 's/,<NON_REF>//g'"
    # exclude multi-allelic sites
    cmd += " | awk '$4 !~ /,/ && $5 !~ /,/'"
    # if length(ALT) > length(REF) and length(REF) > 1, label as COMPLEX_INS
    cmd += " | awk '{ if( length($5) > length($4) && length($4) > 1 ) { print \"COMPLEX_INS\" }"
    # if length(ALT) < length(REF) and length(ALT) > 1, label as COMPLEX_DEL
    cmd += " else if ( length($5) < length($4) && length($5) > 1 ) { print \"COMPLEX_DEL\" } }'"
    cmd += " | sort | uniq -c"
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    x = p.stdout.read()
    for l in x.splitlines():
        l = l.decode("utf-8").strip().split(" ")
        k = "pass_%s" % l[1].lower()
        v = int(l[0])
        r[k] = v
    r["pass_complex_indels"] = r["pass_complex_ins"] + r["pass_complex_del"]

    logging.info("Calculating all_multiallelic_sites...")
    cmd = "bcftools view -H %s | sed 's/,<NON_REF>//g' | awk '$5 ~ /,/' | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["all_multiallelic_sites"] = int(p.stdout.read())

    logging.info("Calculating pass_multiallelic_sites...")
    cmd = "bcftools view -H -f PASS %s | sed 's/,<NON_REF>//g' | awk '$5 ~ /,/' | wc -l" % input_vcf
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    r["pass_multiallelic_sites"] = int(p.stdout.read())

    logging.info("Calculating pass_snp_ts_tv...")
    cmd = "cat %s/pass.bcftools_stats | grep '^TSTV'" % scratch_dir
    logging.debug("CMD: %s" % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    l = p.stdout.read().decode("utf-8").strip().split("\t")
    r["pass_snp_ts_tv"] = np.round(float(l[4]), 2)

    # sanity checks
    logging.info("Running sanity checks for variant counts...")

    error_count = 0
    cmd = "cat %s/all.bcftools_stats | grep '^SN'" % scratch_dir
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

    for l in p.stdout.read().splitlines():
        l = l.decode("utf-8").strip().split("\t")
        k = str(l[2]).lower()
        v = int(l[3])
        if "snps" in k and v != r["all_snps"]:
            logging.error("Mismatch in SNP counts: bcftools stats=%d, count_variants=%d" % (v, r["all_snps"]))
            error_count += 1
        elif "indels" in k and v != r["all_indels"]:
            logging.error("Mismatch in INDEL counts: bcftools stats=%d, count_variants=%d" % (v, r["all_indels"]))
            error_count += 1
        elif "mnps" in k and v != r["all_mnps"]:
            logging.error("Mismatch in MNP counts: bcftools stats=%d, count_variants=%d" % (v, r["all_mnps"]))
            error_count += 1

    if error_count > 0:
        sys.exit("Sanity checks failed, aborting execution")

    return r


def save_output(d, outfile):
    with open(outfile, "w") as f:
        json.dump(d, f, sort_keys=True, indent=4)
        f.write("\n")


def done(scratch_dir):
    if not args.keep_scratch:
        # clean up
        logging.info("Cleaning up...")
        os.remove("%s/variants.vcf.gz" % scratch_dir)
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

    vcf_metrics = count_variants(variants_vcf, args.scratch_dir)
    save_output(vcf_metrics, args.output_json)

    done(args.scratch_dir)

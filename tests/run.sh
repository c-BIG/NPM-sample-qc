#!/bin/bash

export NXF_DEFAULT_DSL=1
export CAPSULE_LOG=none
NXF_ANSI_LOG=false

uid=`uuidgen`
# sampleid="NA12878"
sampleid="NNNNN"

nextflow run NPM-sample-qc/main.nf \
	-config NPM-sample-qc/conf/nextflow-docker.config \
	-params-file NNNNN.yml \
	-resume \
	--outdir s3://npm-grids/gnanakkan/npm-sample-qc/output/${sampleid} \
	-bucket-dir s3://npm-grids/gnanakkan/npm-sample-qc/output/${sampleid}-work \
	-profile awsbatch \
	-with-trace \
	-bg >nf-test.${uid}.log

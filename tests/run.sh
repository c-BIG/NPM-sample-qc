#!/bin/bash

export NXF_DEFAULT_DSL=1
export CAPSULE_LOG=none
NXF_ANSI_LOG=false

uid=`uuidgen`
sampleid="NA12878"

nextflow run NPM-sample-qc/main.nf \
	-config NPM-sample-qc/conf/nextflow-docker.config \
	-params-file NPM-sample-qc/tests/sample_params.aws.bam.yml \
	-resume \
	--outdir s3://npm-grids/gnanakkan/npm-sample-qc/output/${sampleid}-bam-awsbatch \
	-bucket-dir s3://npm-grids/gnanakkan/npm-sample-qc/output/${sampleid}-bam-awsbatch-work \
	-profile awsbatch \
	-bg >nf-test.${uid}.log

#!/bin/bash
REPO_PATH="NPM-sample-qc"

export NXF_DEFAULT_DSL=1
export CAPSULE_LOG=none
NXF_ANSI_LOG=false

uid=`uuidgen`
sampleid="NA12878"

nextflow run $REPO_PATH/main.nf \
	-config $REPO_PATH/nextflow.config \
	-params-file $REPO_PATH/tests/sample_params.aws.bam.yml \
	-resume \
	--outdir s3://npm-grids/${sampleid} \
	-bucket-dir s3://npm-grids/${sampleid}-work \
	-profile awsbatch \
	-bg >nf-test.${uid}.log

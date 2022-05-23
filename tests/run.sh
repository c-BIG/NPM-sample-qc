#!/bin/bash

export NXF_DEFAULT_DSL=1
export CAPSULE_LOG=none

uid=`uuidgen`

nextflow run NPM-sample-qc/main.nf \
	-config NPM-sample-qc/conf/nextflow.config \
	-params-file sample_params.yml \
	-resume \
	--outdir s3://npm-ica-staging-sse/ica-byob/npm-ica-staging-sse/qc-test/${uid}/output \
	-profile local \
	-bg >nf-test.${uid}.log

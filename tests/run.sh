#!/bin/bash
cd test

export NXF_DEFAULT_DSL=1
export CAPSULE_LOG=none
NXF_ANSI_LOG=false

nextflow run ../main.nf \
	-config ../nextflow.config \
	-params-file sample_params.yml \
	-resume \
	-work-dir ./work \
	--outdir ./ \
	-profile docker

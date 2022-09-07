#!/bin/bash

export CAPSULE_LOG=none
NXF_ANSI_LOG=false

nextflow run ../main.nf \
	-config ../nextflow.config \
	-params-file sample_params.yml \
	-resume \
	-work-dir ./work \
	--outdir ./ \
	-profile docker \
	-bg >run.log

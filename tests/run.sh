#!/bin/bash

export CAPSULE_LOG=none
NXF_ANSI_LOG=false

nextflow run ../main.nf \
	-config ../nextflow.config \
	-params-file sample_params.yml \
	-resume \
	--outdir . \
	-work-dir ./work \
	-profile docker \
	-bg > run.log

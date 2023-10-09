#!/bin/bash

# run the workflow on sample NA12878 with below command.
# after successful run completion do "diff ./output/results/NA12878.metrics.json ./output_certified/results/NA12878.metrics.json" to verify the results.

export CAPSULE_LOG=none
NXF_ANSI_LOG=false

nextflow run ../../main.nf \
	-config ../../nextflow.config \
	--inputs_list params.yaml \
	-work-dir ./work \
	--publish_dir ./output

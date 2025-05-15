#!/bin/bash

# run the workflow on sample NA06991 with below command. This unittest sample is selected specifically to show 'cross_contatmination_rate' 0.01%
# after successful run completion do "diff ./output/results/NA06991.metrics.json ./output_certified/results/NA06991.metrics.json" to verify the results.

export CAPSULE_LOG=none
NXF_ANSI_LOG=false

nextflow run ../../main.nf \
	-config ../../nextflow.config \
	--inputs_list params.yaml \
	-work-dir ./work \
	--publish_dir ./output

#!/bin/bash

# run the workflow on sample NA12878 AKT1 gene region with below command.
# after successful run completion do "diff ./output/results/NA12878-chr14-AKT1.metrics.json ./output_certified/results/NA12878-chr14-AKT1.metrics.json" to verify the results.

export CAPSULE_LOG=none
NXF_ANSI_LOG=false

nextflow run ../../main.nf \
	-config ../../nextflow.config -config ./NA12878-chr14-AKT1.config \
        --inputs_list inputs.yaml \
	-work-dir ./work \
	--publish_dir ./output \
        -resume

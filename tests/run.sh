#!/bin/bash

# check if nextflow is available
which nextflow >/dev/null || { echo "ERROR: nextflow is not available"; exit 1; }

# make sure we are in the test directory
ls ../main.nf ../tests  >& /dev/null || { echo "ERROR: must run from NPM-sample-qc/tests directory"; exit 1; }

# define directories
projectdir=$(realpath "$(pwd)/..")
outdir="$1"
if [ -z "$outdir" ]; then
    echo "ERROR: Missing outdir argument, e.g.:" 1>&2
    echo "  outdir='/data/13000026/pipeline/dev/NPM-sample-qc-aux/tests/<my_run>'" 1>&2
    echo "  ./run.sh \$outdir"
    exit 1
fi
workdir="${outdir}/work"

# define profile
profile="$2"
if [ -z "$profile" ]; then
    profile="nscc"
fi
echo "Running with profile \"${profile}\"."
echo "This option can be changed when calling the script:"
echo "  ./run.sh <outdir> <profile>"

# run nextflow
nextflow run ${projectdir}/main.nf \
    -config ${projectdir}/conf/nextflow.config \
    -params-file $(pwd)/sample_params.yml \
    -profile ${profile} \
    -resume \
    -work-dir ${workdir} \
    --outdir ${outdir} \
    --keep_workdir

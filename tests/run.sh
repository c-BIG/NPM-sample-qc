#!/bin/bash

# check if nextflow is available
which nextflow >/dev/null || { echo "ERROR: nextflow is not available"; exit 1; }

# make sure we are in the test directory
ls ../main.nf ../tests  >& /dev/null || { echo "ERROR: must run from NPM-sample-qc/tests directory"; exit 1; }

# define directories
testdir="$1"
if [ -z "$testdir" ]; then
    echo "ERROR: Missing testdir argument, e.g.:" 1>&2
    echo "  testdir='/data/13000026/pipeline/dev/NPM-sample-qc/tests/<my_run>'" 1>&2
    echo "  ./run.sh \$testdir"
    exit 1
fi
workdir=$testdir/work
publishdir=$testdir/results
projectdir=$(realpath "$(pwd)/..")

# define profile
profile="$2"
if [ -z "$profile" ]; then
    profile="standard"
fi
echo "Running with profile \"${profile}\"."
echo "This option can be changed when calling the script:"
echo "  ./run.sh <testdir> <profile>"

# run nextflow
nextflow run ${projectdir}/main.nf \
    -config ${projectdir}/nextflow.config \
    -params-file $(pwd)/sample_params.yml \
    -w ${workdir} \
    --publishdir ${publishdir} \
    -profile ${profile} -resume --keep_workdir

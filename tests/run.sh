#!/bin/bash

# check if nextflow is available
which nextflow >/dev/null || { echo "ERROR: nextflow is not available"; exit 1; }

# make sure we are in the test directory
ls ../main.nf ../tests  >& /dev/null || { echo "ERROR: must run from NPM-sample-qc/tests directory"; exit 1; }

# define directories
testdir="$1"
if [ -z "$testdir" ]; then
    echo "ERROR: Missing testdir argument. This is where results and workdir will be created." 1>&2
    echo "You can for example use: testdir='/data/13000026/pipeline/dev/NPM-sample-qc/tests/<my_run>'" 1>&2
    echo "Then run: ./run.sh \$basedir"
    exit 1
fi
workdir=$testdir/work
publishdir=$testdir/results
projectdir=$(realpath "$(pwd)/..")

# run nextflow
nextflow run ${projectdir}/main.nf \
    -config ${projectdir}/nextflow.config \
    -params-file ${projectdir}/tests/sample_params.yaml \
    -w ${workdir} \
    --publishdir ${publishdir} \
    -profile standard -resume --keep_workdir

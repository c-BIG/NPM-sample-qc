#!/bin/bash

# check if nextflow is available
which nextflow >/dev/null || { echo "ERROR: nextflow is not available"; exit 1; }

# make sure we are in the test directory
ls ../main.nf ../tests  >& /dev/null || { echo "ERROR: must run from NPM-sample-qc/tests directory"; exit 1; }

# define directories
basedir="$1"
if [ -z "$basedir" ]; then
    echo "ERROR: Missing basedir argument. This is where results and workdir will be created." 1>&2
    echo "You can for example use: basedir='/data/13000026/pipeline/dev/NPM-sample-qc/tests'" 1>&2
    exit 1
fi
suffix="run-$(date +%Y%m%d-%H%M)"
tmpdir=$(mktemp -d ${basedir}/${suffix}-XXXXXX) || exit  1
projectdir=$(realpath "$(pwd)/..")
workdir=$tmpdir/work
publishdir=$tmpdir/results

# run nextflow
nextflow run ${projectdir}/main.nf \
    -config ${projectdir}/nextflow.config \
    -params-file ${projectdir}/tests/sample_params.yaml \
    -w ${workdir} \
    --publishdir ${publishdir} \
    -profile standard -resume --keep_workdir

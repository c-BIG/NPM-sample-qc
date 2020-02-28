#!/bin/bash

# run from the project root
ls main.nf >& /dev/null || { echo "ERROR: must run from NPM-sample-qc project root"; exit 1; }

# define directories
projectdir=$(realpath "$(pwd)")
targetdir="$1"
if [ -z "$targetdir" ]; then
    echo "ERROR: Missing targetdir argument, e.g.:" 1>&2
    echo "  targetdir='/data/13000026/pipeline/dev/NPM-sample-qc'" 1>&2
    echo "  ./deploy.sh \$targetdir"
    exit 1
fi

# automatically resolve version
version=$(git describe --tags)
if [[ $version != v[0-9]* ]]
then
    version=$(git rev-parse --short HEAD)
fi

# define final target
finaltarget="$targetdir/NPM-sample-qc-$version"

# deploy
if [[ -d "$finaltarget" ]]
then
    echo "Deployment directory already exists, will not overwrite: $finaltarget"
else
    echo "Deploying NPM-sample-qc into $finaltarget"
    rsync -a --exclude "." $projectdir/* $finaltarget
fi
echo "DONE"

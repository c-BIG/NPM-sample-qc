#!/bin/bash

# AWS
if [[ ${USER} = "ubuntu" && ${HOSTNAME} =~ "ip" ]]
then
    sudo singularity build --force npm-sample-qc.simg singularity.def
    singularity exec npm-sample-qc.simg echo "Build complete"
    aws s3 cp npm-sample-qc.simg s3://npm-sample-qc
fi

# NSCC
if [ ${HOSTNAME} = "gis01" ]
then
    BASEDIR="/data/13000026/pipeline/dev/NPM-sample-qc/containers"
    aws s3 cp s3://npm-sample-qc/npm-sample-qc.simg ${BASEDIR}
    singularity exec ${BASEDIR}/npm-sample-qc.simg echo "Deployment complete"
fi

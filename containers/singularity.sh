#!/bin/bash

# AWS
if [[ ${USER} = "ubuntu" && ${HOSTNAME} =~ "ip" ]]
then
    sudo singularity build --force npm-sample-qc.simg singularity.def
    aws s3 cp npm-sample-qc.simg s3://npm-sample-qc
fi

# NSCC
if [ ${HOSTNAME} = "gis01" ]
then
    aws s3 cp s3://npm-sample-qc/npm-sample-qc.simg .
fi

#!/bin/bash

# AWS
if [[ ${USER} = "ubuntu" && ${HOSTNAME} =~ "ip" ]]
then
    sudo singularity build --force npm-sample-qc.simg singularity.def
    singularity exec npm-sample-qc.simg echo "Build complete"
fi
### cp npm-sample-qc.simg npm-sample-qc-20052022.simg
### aws s3 cp npm-sample-qc-20052022.simg s3://npm-sample-qc/

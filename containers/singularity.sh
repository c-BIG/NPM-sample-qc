#!/bin/bash

# AWS
if [[ ${USER} = "ubuntu" && ${HOSTNAME} =~ "ip" ]]
then
    sudo singularity build --force npm-sample-qc.simg singularity.def
    singularity exec npm-sample-qc.simg echo "Build complete"
fi

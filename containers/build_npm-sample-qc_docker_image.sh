#!/bin/bash

# Create a docker image named 'npm-sample-qc:latest'
# This will be used when specify profile 'docker'.

# Test OS
if [[ ${OSTYPE} = "darwin"* ]]
then # MacOS
    docker build --platform linux/amd64 -t npm-sample-qc .
else # Linux / WSL
  docker build -t npm-sample-qc .
fi

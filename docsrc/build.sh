#!/bin/bash

# rm existing docs if present
[ -d api ] && rm -r api
[ -d build ] && rm -r build

# auto-generate metrics docs
sphinx-apidoc --output-dir api --no-toc ../bin/metrics

# create sphinx docs
make github

# copy images
cp ../npm-sample-qc-overview.PNG ../docs/

# clean up docsrc
[ -d api ] && rm -r api
[ -d build ] && rm -r build

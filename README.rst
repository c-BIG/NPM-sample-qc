=============
NPM-sample-qc
=============

NPM-sample-qc is a Nextflow_ workflow to obtain QC metrics from single-sample WGS results. It has been created to support QC efforts within the National Precision Medicine programme in Singapore (NPM), a local initiative that intends to sequence the genomes of 100K individuals (SG100K), but can be easily extended to other large-scale sequencing projects.

.. _Nextflow: https://www.nextflow.io/

Requirements
============

* `Install Nextflow`_
* `Install Docker`_
* Install and configure `AWS CLI`_

.. _Install Nextflow: https://www.nextflow.io/docs/latest/getstarted.html#installation
.. _Install Docker: https://docs.docker.com/get-docker/
.. _AWS CLI: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html

Quick start
===========

Clone this repository ::

  git clone git@github.com:c-BIG/NPM-sample-qc.git

Build docker image locally ::

  # Move to containers
  cd NPM-sample-qc/containers
  # Build docker image locally
  sh build_npm-sample-qc_docker_image.sh
  # Move back to project root
  cd ../

Run workflow on sample NA12878 from the 1000 Genomes Phase 3 Reanalysis with DRAGEN 3.7 ::

  nextflow run main.nf \
    -config      nextflow.config \
    -profile     docker \
    -params-file tests/sample_params.yml \
    -work-dir    test-run/work \
    --outdir     test-run

Please refer to the workflow help for more information on its usage and access to additional options: ::

  nextflow run NPM-sample-qc/main.nf --help

Understanding the workflow
==========================

Resources
---------

The workflow requires the following resources given in the ``conf/resources.config``

- *N-regions reference file*, used as an input for mosdepth. This file can be downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz (``2019-03-11 09:51, 12K``).         

- *Human Reference Genome FASTA file*, used as an input for multiple tools. This file can be downloaded from ``s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta``.

- *FASTA file index*. This file can be downloaded from ``s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai`` and not required to be specified in the config. The workflow will look fasta index file in a folder the fasta file is present.

Inputs
------

Input requirements can be split into two categories:

- **Generic workflow settings** specify parameters that will not vary from run to run, e.g. Nextflow profile declarations, trace/timeline/report/dag options, output structure and paths to data resources. See ``nextflow.config`` for additional details.

- **Sample-specific settings** contain paths to WGS data for a given sample, namely BAM/CRAM. The workflow expects the BAM/CRAM index (bai/crai) to be present in the same location. See ``tests/sample_params.yml`` for an example.

.. _Nextflow configuration: https://www.nextflow.io/docs/latest/config.html

Outputs
-------

Upon completion, the workflow will create the following files in the ``outdir`` directory: ::

  /path/to/outdir/
      pipeline_info/    # dag, timeline and trace files
          dag.pdf
          timeline.html
          report.html
          trace.txt
      results/          # final metrics.json and intermediate outputs
          <sample_id>.metrics.json    
          samtools/
          picard/
          mosdepth/
          multiqc/

If ``keep_workdir`` has been specified, the contents of the Nextflow work directory (``work-dir``) will also be preserved.

Workflow logic
==============

We provide a schematic representation of the workflow in the figure below:
  
.. raw:: html

   <img src="npm-sample-qc-overview.PNG" width="500px"/>   

In a nutshell, this workflow generates QC metrics from single-sample WGS results in three stages: **metrics calculation**, **parsing of intermediate outputs** and **generation of a final report**. This makes it possible to take full advantage of the parallelisation capabilities of Nextflow, allows users to leverage third-party tools or add custom scripts, and enables auto-documentation of metrics from code comments.

**Metrics calculation**

The current workflow combines widely-used third-party tools (samtools, picard, mosdepth) and custom scripts. Full details on which processes are run/when can be found in the actual workflow definition (``main.nf``). We also provide an example dag for a more visual representation (``tests/example_dag.pdf``).

**Metrics parsing**

Next, output files from each individual tool are parsed and combined into a single json file. This is done by calling ``bin/multiqc_plugins/multiqc_npm/``, a MultiQC plugin that extends the base tool to support additional files.

**Metrics reporting**

Finally, the contents of the MultiQC json are formatted into a final metrics report, also in json format. The reporting logic lives in the ``bin/compile_metrics.py`` script, and whilst its contents are simple, it enables automatic documentation of metric definitions from code comments (see the **Metric definitions** section).

Metric definitions
==================

The full list of metrics reported by this workflow and details on how they've been calculated can be found here_.

.. _here: https://c-big.github.io/NPM-sample-qc/metrics.html

When needed, page contents can be updated by running the following command: ::

  # Install sphinx
  pip install sphinx sphinx_rtd_theme sphinx_automodapi
  # Move to doc source
  cd docsrc
  # Build the doc
  ./build.sh

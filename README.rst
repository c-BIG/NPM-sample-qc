=============
NPM-sample-qc
=============

NPM-sample-qc is a Nextflow_ workflow to obtain QC metrics from single-sample WGS results. It has been created to support QC efforts within the National Precision Medicine programme in Singapore (NPM), a local initiative that intends to sequence the genomes of 100K individuals (SG100K), but can be easily extended to other large-scale sequencing projects.

.. _Nextflow: https://www.nextflow.io/


Quick start
===========

Here an example command to launch the workflow: ::

  nextflow run main.nf \
  -config nextflow.config \
  -params-file tests/sample_params.yml \
  -work-dir ./work \
  --outdir ./ \
  --keep_workdir

Please refer to the workflow help for more information on its usage and access to additional options: ::

  nextflow run main.nf --help


Understanding the workflow
==========================

Required inputs
---------------

NPM-sample-qc input requirements can be split into two categories:

- **Generic workflow settings** specify parameters that will not vary from run to run, e.g. Nextflow profile declarations, trace/timeline/report/dag options, output structure and paths to data resources. See ``nextflow.config`` for additional details.

- **Sample-specific settings** contain paths to WGS results for a given sample, namely CRAM/BAM. See ``tests/sample_params*.yml`` for an example.

.. _Nextflow configuration: https://www.nextflow.io/docs/latest/config.html

Output files
------------

Upon completion, NPM-sample-qc will create the following files in the ``outdir`` directory: ::

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

Resources
=========

Other than than the sample specific BAM/CRAM file, it requires the follwing data resources.
Human reference genome GRCh38_ in fasta format, length of the autosomes and formatted gapped regions of the assembly in bed format. Which can be downloaded from UCSC_ See ``resources/Homo_sapiens_assembly38.autosomes.bed`` and ``resources/Homo_sapiens_assembly38.autosomes.n_regions.bed`` for an example.

.. _GRCh38: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
.. _UCSC: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz

AWS batch deployment
====================

Edit the following files to suit your AWS batch configuration  
* conf/awsbatch.config
*

Limitation with ``awsbatch`` 'Executor', write tracing & visualization files locally and upload it to s3. Issue_

.. _Issue: https://github.com/nextflow-io/nextflow/issues/2342

Workflow logic
--------------

We provide a schematic representation of the workflow in the figure below:
  
.. raw:: html

   <img src="docs/npm-sample-qc-overview.PNG" width="500px"/>   

In a nutshell, NPM-sample-qc generates QC metrics from single-sample WGS results in three stages: metrics calculation, parsing of intermediate outputs and generation of a final report. This makes it possible to take full advantage of the parallelisation capabilities of Nextflow, allows users to leverage third-party tools or add custom scripts, and enables auto-documentation of metrics from code comments.

**Metrics calculation**

The current workflow combines widely-used third-party tools (samtools, picard) and custom scripts. Full details on which processes are run/when can be found in the actual workflow definition (``main.nf``). We also provide an example dag for a more visual representation (``tests/example_dag.pdf``).


**Metrics parsing**

Next, output files from each individual tool are parsed and combined into a single json file. This is done by calling MultiQC_NPM_, a MultiQC_ plugin that extends the base tool to support additional files.

.. _MultiQC_NPM: https://github.com/c-BIG/MultiQC_NPM/
.. _MultiQC: https://github.com/ewels/MultiQC

**Metrics reporting**

Finally, the contents of the MultiQC json are formatted into a final metrics report, also in json format. The reporting logic lives in the compile_metrics.py script, and whilst its contents are simple, it enables automatic documentation of metric definitions from code comments (see the **Metric definitions** section).


Metric definitions
==================

The full list of metrics reported by the NPM-sample-qc workflow and details on how they've been calculated can be found here_.

.. _here: https://c-big.github.io/NPM-sample-qc/metrics.html

When needed, page contents can be updated by running the following command: ::

  cd docsrc; ./build.sh
  

AWS batch deployment
====================

Edit the following files to suit your AWS batch configuration  
* conf/awsbatch.config


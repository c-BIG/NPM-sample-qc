=============
NPM-sample-qc
=============

NPM-sample-qc is a Nextflow_ workflow to obtain QC metrics from single-sample WGS results. It has been created to support QC efforts within the National Precision Medicine programme in Singapore (NPM), a local initiative that intends to sequence the genomes of 100K individuals (SG100K), but can be easily extended to other large-scale sequencing initiatives.

.. _Nextflow: https://www.nextflow.io/


Quick start
===========

Here an example command to launch the workflow: ::

  nextflow run main.nf \
  -config conf/nextflow.config \
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

- **Generic workflow settings** (``conf/nextflow.config``) specify parameters that will not vary from run to run, e.g. nextflow profile declarations, trace/timeline/dag generation, output structure, paths to data resources...

- **Sample-specific settings** (e.g. ``tests/sample_params.yml``) contain paths to WGS results for a given sample, namely CRAM and VCF/gVCFs. Optionally, you can also provide a positive sample tracking VCF (``pst_vcf``) to calculate genotype concordance against your WGS VCF (see :ref:`Metric definitions`).

NPM users will only be required to update ``sample_params.yml`` to match their sample of interest. In addition, any user can update the global settings by following standard `Nextflow configuration`_ guidelines.

.. _Nextflow configuration: https://www.nextflow.io/docs/latest/config.html

Output files
------------

Upon completion, NPM-sample-qc will create the following files in the ``outdir`` directory: ::

  /path/to/outdir/
      pipeline_info/    # dag, timeline and trace files
          dag.pdf
          timeline.html
          trace.txt
      results/          # final metrics.json and intermediate outputs
          <sample_id>.metrics.json    
          bcftools/
          count_variants/
          verifybamid2/
          samtools/
          plot_bamstats/
          picard/
          multiqc/

If ``keep_workdir`` has been specified, the contents of the Nextflow work directory (``work-dir``) will also be preserved.

Workflow logic
--------------

We provide a schematic representation of the workflow in the figure below:

.. image:: npm-sample-qc-overview.PNG

In a nutshell, the NPM-sample-qc workflow generates QC metrics from single-sample WGS results in three stages: (1) metrics calculation, (2) parsing of intermediate outputs and (3) generation of a final report. This makes it possible to take full advantage of the parallelisation capabilities of Nextflow, allows users to leverage third-party tools or add custom scripts, and enables auto-documentation of metrics from code comments.

**Metrics calculation**

The current workflow combines well-established third-party tools (samtools, bcftools, picard, VerifyBamID2) and custom scripts (e.g. count_variants.py) to obtain a rich set of QC metrics. Full details on which processes are run/when can be found in the actual workflow definition (``main.nf``). We also provide an example dag (``tests/example_dag.pdf``) for a more visual representation.


**Metrics parsing**

Next, output files for each individual metrics tool are parsed and combined into a single json file. This is done by calling MultiQC_NPM_, a MultiQC_ plugin that extends the base tool to support additional metrics outputs (e.g. bcftools gtcheck, picard CollectQualityYieldMetrics, count_variants.py...).

.. _MultiQC_NPM: https://github.com/c-BIG/MultiQC_NPM/
.. _MultiQC: https://github.com/ewels/MultiQC

**Metrics reporting**

Finally, the contents of the MultiQC json are formatted into a final metrics report, also in json format. The reporting logic lives in the compile_metrics.py script, and whilst its contents are simple, it enables automatic documentation of metric definitions from code comments (see :ref:`Metric definitions`).

.. _Metric definitions:

Metric definitions
==================

The full list of metrics reported by the NPM-sample-qc workflow can be found here_. Each metric has been documented with its definition and details on how it has been calculated.

.. _here https://c-big.github.io/NPM-sample-qc/metrics.html

When needed, page contents can be updated by running the following command: ::

  cd docsrc ; ./build.sh

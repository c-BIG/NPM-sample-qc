#!/usr/bin/env python3

"""
NPMExtraMetrics module
"""

from __future__ import print_function
import logging

from multiqc.utils import config
from multiqc.modules.base_module import BaseMultiqcModule

from . import npm_picard_QualityYieldMetrics
from . import npm_samtools_stats_bq
from . import npm_bcftools_gtcheck
from . import npm_count_variants
from . import npm_mosdepth

log = logging.getLogger('multiqc')


class NPMExtraMetrics(BaseMultiqcModule):

    def __init__(self):

        # Halt execution if we've disabled the plugin
        if config.kwargs.get('enable_plugin', False) is True \
               or getattr(config, 'enable_plugin', False) is True:
            log.info("Running with MultiQC NPM plugin")
        else:
            log.info("Skipping MultiQC NPM plugin as not enabled")
            return None

        # Initialise the parent module Class object
        super(NPMExtraMetrics, self).__init__(
            name='NPMExtraMetrics',
            anchor='npm_extra_metrics'
        )

        # Call each submodule and report progress
        n = dict()

        n['npm_picard_QualityYieldMetrics'] = npm_picard_QualityYieldMetrics.parse_reports(self)
        if n['npm_picard_QualityYieldMetrics'] > 0:
            log.info("Found %d npm_picard_QualityYieldMetrics reports" % n['npm_picard_QualityYieldMetrics'])

        n['npm_samtools_stats_bq'] = npm_samtools_stats_bq.parse_reports(self)
        if n['npm_samtools_stats_bq'] > 0:
            log.info("Found %d npm_samtools_stats_bq reports" % n['npm_samtools_stats_bq'])

        n['npm_bcftools_gtcheck'] = npm_bcftools_gtcheck.parse_reports(self)
        if n['npm_bcftools_gtcheck'] > 0:
            log.info("Found %d npm_bcftools_gtcheck reports" % n['npm_bcftools_gtcheck'])

        n['npm_count_variants'] = npm_count_variants.parse_reports(self)
        if n['npm_count_variants'] > 0:
            log.info("Found %d npm_count_variants reports" % n['npm_count_variants'])

        n['npm_mosdepth'] = npm_mosdepth.parse_reports(self)
        if n['npm_mosdepth'] > 0:
            log.info("Found %d npm_mosdepth reports" % n['npm_mosdepth'])

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

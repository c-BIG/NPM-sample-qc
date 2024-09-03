#!/usr/bin/env python3

"""
Command line options

See https://multiqc.info/docs/#command-line-options
"""

import click

enable_plugin = click.option('--enable-npm-plugin',
                             'enable_plugin',
                             is_flag=True,
                             default=False,
                             help="Enable the MultiQC NPM plugin on this run"
                            )

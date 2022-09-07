#!/usr/bin/env python3

"""
A MultiQC plugin for downstream analysis of NPM data
"""

from setuptools import setup, find_packages

version = '0.1'

setup(
    name='multiqc_npm',
    version=version,
    author='Mar Gonzalez-Porta',
    author_email='Mar_Gonzalez_Porta@gis.a-star.edu.sg',
    description="MultiQC plugin for NPM",
    long_description=__doc__,
    keywords='bioinformatics',
    url='https://github.com/c-BIG/NPM-sample-qc',
    download_url='https://github.com/c-BIG/NPM-sample-qc/releases',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'multiqc',
        'pandas'
    ],
    entry_points={
        'multiqc.modules.v1': [
            'multiqc_npm = multiqc_npm.modules.npm_extra_metrics:NPMExtraMetrics'
        ],
        'multiqc.cli_options.v1': [
            'enable-npm-plugin = multiqc_npm.cli:enable_plugin'
        ],
        'multiqc.hooks.v1': [
            'before_config = multiqc_npm.hooks:before_config'
        ]
    }
)

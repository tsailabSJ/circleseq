#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

setup(
    name='circleseq',
    version='1.1',
    description="An easy to use bioinformatic pipeline for the CIRCLE-seq assay.",
    author="Shengdar Q Tsai, Martin Aryee, Ved V Topkar, Jose Malagon-Lopez",
    author_email='STSAI4@mgh.harvard.edu, Aryee.Martin@mgh.harvard.edu, vedtopkar@gmail.com, jose.lopez@mail.harvard.edu',
    url='https://github.com/tsailabSJ/circleseq',
    packages=[
        'circleseq',
    ],
    package_dir={'circleseq':
                 'circleseq'},
    include_package_data=True,
    install_requires=requirements,
    license="AGPL",
    keywords='circleseq',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: Unix',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7'
    ]
)
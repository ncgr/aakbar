# -*- coding: utf-8 -*-
#
# This file is part of aakbar.
# Copyright (C) 2018, NCGR.
#
# aakbar is free software; you can redistribute it and/or modify
# it under the terms of the 3-Clause BSD License; see LICENSE.txt
# file for more details.
#
''' aakbar -- amino-acid k-mer signature tools '''
#
# Developers, install with:
#    pip install -r requirements.txt
#    python setup.py develop
#
from setuptools import setup, find_packages
import sys
# restrict to python 3.4 or later
if sys.version_info < (3,4,0,'final',0):
    raise SystemExit('Python 3.4 or later is required!')

NAME = 'aakbar'

tests_require = [
    'check-manifest>=0.25',
    'coverage>=4.0',
    'isort>=4.2.2.2',
    'pydocstyle>=1.0.0',
    'pytest-cache>=1.0',
    'pytest-cov>=1.8.0',
    'pytest-pep8>=1.0.6',
    'pytest>=2.8.0'
]
extras_require = dict(docs=['Sphinx>=1.4.2'], tests=tests_require)
extras_require['all'] = []
for reqs in extras_require.values():
    extras_require['all'].extend(reqs)
packages = find_packages()

setup(
    description=__doc__,
    packages=packages,
    setup_requires=['setuptools>30.3.0',
                    'setuptools-scm>1.5'
                    ],
    install_requires=['Click',],
    entry_points={
        'console_scripts': [NAME + ' = ' + NAME + ':cli']
    },
    use_scm_version={
        'version_scheme': 'guess-next-dev',
        'local_scheme': 'dirty-tag',
        'write_to': NAME + '/version.py'
    },
    extras_require=extras_require,
    tests_require=tests_require,
)
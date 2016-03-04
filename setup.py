# -*- coding: utf-8 -*-
'''
aakbar -- amino-acid k-mer signature tools
'''

# Packagers: aakbar's version is auto-generated using setuptools-scm.
# The resulting file (version.py) is not managed by Git.
# As a consequence, packagers should not use the GitHub
# tarballs, but rather the PyPI ones.

# Developers:
# Install with
# pip install --editable .
# and execute as a module.


from setuptools import setup

setup(
    name='aakbar',
    use_scm_version={
        'write_to': 'aakbar/version.py'
    },
    setup_requires=['setuptools_scm'],
    url='http://github.com/GenerisBio/aakbar',
    keywords=['science', 'biology', 'bioinformatics', 'phylogenomics', 'peptide', 'signatures'],
    license='BSD',
    description='Amino-Acid k-mer phylogenetic signature tools',
    long_description=open('README.rst').read(),
    author='Joel Berendzen',
    author_email='joel@generisbio.com',
    packages=['aakbar'],
    include_package_data=True,
    zip_safe=False,
    install_requires=['bcbio-gff', 'biopython', 'bitarray',
                      'click>=5.0','click_plugins', 'networkx',
                      'numpy', 'pandas'],
    entry_points={
                 'console_scripts':['aakbar = aakbar:cli']
                },
    classifiers=[
                        'Development Status :: 3 - Alpha',
                        'Environment :: Console',
                        'Environment :: MacOS X',
                        'Environment :: Win32 (MS Windows)',
                        'Intended Audience :: Science/Research',
                        'License :: Other/Proprietary License ',
                        'Operating System :: OS Independent',
                        'Programming Language :: Python',
                        'Topic :: Scientific/Engineering :: Bio-Informatics',
                        ]
)
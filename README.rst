aakbar
======
Amino-Acid k-mer tools for creating, searching, and analyzing phylogenetic signatures from genomes or reads of DNA.

Prerequisites
-------------
A 64-bit Python 3.4 or greater is required.  8 GB or more of memory is recommended.

The python dependencies of aakbar are: biopython, click>=5.0, click_plugins numpy, pandas, pyfaidx,
and pyyaml.  Running the examples also requires the `pyfastaq  https://pypi.python.org/pypi/pyfastaq`
package.

If you don't have a python installed that meets these requirements, I recommend getting
`Anaconda Python <https://www.continuum.io/downloads>` on MacOSX and Windows for the smoothness
of installation and for the packages that come pre-installed.  Once Anaconda python is installed,
you can get the dependencies like this on MacOSX::

    export PATH=~/anaconda/bin:${PATH}    # you might want to put this in your .profile
    conda install click
    conda install --channel https://conda.anaconda.org/IOOS click-plugins
    conda install --channel https://conda.anaconda.org/bioconda pyfaidx
    conda install --channel https://conda.anaconda.org/bioconda pyfastaq  # only required for examples


Installation
------------
This package is tested under Linux and MacOS using Python 3.5 and is available from the PyPI.  To
install via pip (or pip3 under some distributions) : ::

     pip install aakbar

If you wish to develop aakbar,  download a `release <https://github.com/ncgr/aakbar/releases>`_
and in the top-level directory: ::

	pip install --editable .

If you wish to have pip install directly from git, use this command: ::

	pip install git+https://github.com/ncgr/aakbar.git



Usage
-----
Installation puts a single script called ``aakbar`` in your path.  The usage format is::

    aakbar [GLOBALOPTIONS] COMMAND [COMMANDOPTIONS] [ARGS]

A listing of commands is available via ``aakbar --help``.  Current available commands are:

============================= ====================================================
  calculate-peptide-terms     Write peptide terms and histograms.
  conserved-signature-stats   Stats on signatures found in all input genomes.
  define-set                  Define an identifier and directory for a set.
  define-summary              Define summary directory and label.
  demo-simplicity             Demo self-provided simplicity outputs.
  filter-peptide-terms        Remove high-simplicity terms.
  init-config-file            Initialize a configuration file.
  install-demo-scripts        Copy demo scripts to the current directory.
  intersect-peptide-terms     Find intersecting terms from multiple sets.
  label-set                   Define label associated with a set.
  peptide-simplicity-mask     Lower-case high-simplicity regions in FASTA.
  search-peptide-occurrances  Find signatures in peptide space.
  set-letterfreq-window       Define size of letterfreq window.
  set-plot-type               Define label associated with a set.
  set-simplicity-object       Select simplicity-calculation object.
  show-config                 Print location and contents of config file.
  show-context-object         Print the global context object.
  test-logging                Logs at different severity levels.
============================= ====================================================

Examples
--------

Bash scripts that implement examples for calculating and using signature sets for
Firmicutes and Streptococcus, complete with downloading data from GenBank, will
be created in the (empty) current working directory when you issue the command:

    aakbar install-demo-scripts

On linux and MacOS, follow the instructions to run the demos.  On Windows, you will
need ``bash`` installed for the scripts to work.


+-------------------+------------+------------+
| Latest Release    | |pypi|     | |akbar|    |
+-------------------+------------+            +
| GitHub            | |repo|     |            |
+-------------------+------------+            +
| License           | |license|  |            |
+-------------------+------------+            +
| Documentation     | |rtd|      |            |
+-------------------+------------+            +
| Travis Build      | |travis|   |            |
+-------------------+------------+            +
| Coverage          | |coverage| |            |
+-------------------+------------+            +
| Code Grade        | |codacy|   |            |
+-------------------+------------+            +
| Dependencies      | |pyup|     |            |
+-------------------+------------+            +
| Issues            | |issues|   |            |
+-------------------+------------+------------+


.. |akbar| image:: docs/akbar-the-great.jpg
     :target: https://en.wikipedia.org/wiki/Akbar
     :alt: Akbar the Great

.. |pypi| image:: https://img.shields.io/pypi/v/aakbar.svg
    :target: https://pypi.python.org/pypi/aakbar
    :alt: Python package

.. |repo| image:: https://img.shields.io/github/commits-since/ncgr/aakbar/0.1.svg
    :target: https://github.com/LegumeFederation/lorax
    :alt: GitHub repository

.. |license| image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
    :target: https://github.com/ncgr/aakbar/blob/master/LICENSE.txt
    :alt: License terms

.. |rtd| image:: https://readthedocs.org/projects/aakbar/badge/?version=latest
    :target: http://aakbar.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Server

.. |travis| image:: https://img.shields.io/travis/ncgr/aakbar.svg
    :target:  https://travis-ci.org/ncgr/aakbar
    :alt: Travis CI

.. |codacy| image:: https://api.codacy.com/project/badge/Grade/75ebc8405ee74953a555a51abe16d9fa
    :target: https://www.codacy.com/manual/joelb123/aakbar?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ncgr/aakbar&amp;utm_campaign=Badge_Grade
    :alt: Codacy.io grade

.. |coverage| image:: https://codecov.io/gh/ncgr/aakbar/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/ncgr/aakbar
    :alt: Codecov.io test coverage

.. |issues| image:: https://img.shields.io/github/issues/ncgr/aakbar.svg
    :target:  https://github.com/ncgr/aakbar/issues
    :alt: Issues reported

.. |pyup| image:: https://pyup.io/repos/github/ncgr/aakbar/shield.svg
     :target: https://pyup.io/repos/github/ncgr/aakbar/
     :alt: pyup.io dependencies


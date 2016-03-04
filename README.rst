|Travis| |PyPI| |Landscape| |Coveralls| |Depsy|
======
aakbar
======
Amino-Acid k-mer tools for creating, searching, and analyzing phylogenetic and functional signatures in DNA sequences.

Prerequisites
-------------
A 64-bit Python 3.4 or greater is required.  32 GB or more of memory is recommended.

The python dependencies of aakbar are: biopython, bcbio-gff, click>=5.0, click_plugins numpy, pandas, networkx, and setuptools_scm.

If you don't have a python installed that meets these requirements, I recommend getting `Anaconda Python <https://www.continuum.io/downloads>` on MacOSX and Windows for the smoothness of its installation and the large number of useful packages that come pre-installed.  Installation of Anaconda python and other prerequisites on MacOSX goes like this:

::
	bash Anaconda3-2.5.0-MacOSX-x86_64.sh # version number may change
	export PATH=~/anaconda/bin:${PATH}    # you might want to put this in your .bashrc
	#
        # install required stuff
	#
	conda install gcc setuptools_scm
	#
	# add some useful channels for items not in default repository
	#
	conda config --add channels anaconda
	conda config --add channels r
	conda config --add channels bioconda
	conda install click 
	#
	# click-plugins currently in IOOS channel.  If not, search for it thus:
	# anaconda search -t conda click-plugins
	#
	conda install --channel https://conda.anaconda.org/IOOS click-plugins
	#
	# bcbio-gff and pyfaidx are from pip (the pyfaidx in bioconda is too old)
	#
	pip install bcbio-gff pyfaidx
	#
	# the rest are NOT dependencies for aakbar, but useful for bioinformatics
	#
	conda install bamtools bedtools bioawk pysam sambamba samtools scikit-bio


Installation
------------

This package is tested under Linux and MacOS using Python 3.5 and is available from the PyPI:

::
	pip install aakbar

If you wish to develop aakbar,  download a `release <https:/github.com/generisbio/aakbar/releases>`_ and in the top-level directory:

::
	python install --editable .


Basic Use
---------
aakbar is implemented as a command-line program with subcommands.  To list these subcommands: ::

    aakbar --help

Documentation
-------------
- `Readthedocs documentation <https://aakbar.readthedocs.org/en/latest/index.html>`_


License
-------
aakbar is distributed under a `BSD License`_.  Proprietary modules that extend the capabilities of aakbar
are available from `GenerisBio, LLC <http://www.generisbio.com>`_.

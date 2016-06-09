aakbar
======
Amino-Acid k-mer tools for creating, searching, and analyzing phylogenetic signatures from genomes and from short reads of DNA.

Prerequisites
-------------
A 64-bit Python 3.4 or greater is required.  8 GB or more of memory is recommended.

The python dependencies of aakbar are: biopython, click>=5.0, click_plugins numpy, pandas, pyfaidx,
and pyyaml.  Running the examples also requires pyfastaq

If you don't have a python installed that meets these requirements, I recommend getting
`Anaconda Python <https://www.continuum.io/downloads>` on MacOSX and Windows for the smoothness
of installation and for the packages that come pre-installed.  Once Anaconda python is installed,
you can get the dependencies like this on MacOSX::

	export PATH=~/anaconda/bin:${PATH}    # you might want to put this in your .profile
	conda install --channel https://conda.anaconda.org/IOOS click-plugins
        conda install --channel https://conda.anaconda.org/bioconda pyfastaq


Installation
------------
This package is tested under Linux and MacOS using Python 3.5 and is available from the PyPI: ::

     pip install aakbar

If you wish to develop aakbar,  download a `release <https://github.com/ncgr/aakbar/releases>`_
and in the top-level directory: ::

	pip install --editable .

If you wish to have pip install directly from git, use this command: ::

	pip install git+https://github.com/ncgr/aakbar.git#egg=proj

 


Basic Use
---------
aakbar is implemented as a command-line program with subcommands.  To list these subcommands: ::

    aakbar --help

Documentation
-------------
- `Readthedocs documentation <https://aakbar.readthedocs.org/en/latest/index.html>`_


License
-------
aakbar is distributed under a `BSD License`.

This directory contains some scripts for downloading data from GenBank, then
calculating and evaluating sets of peptide signatures. The scripts have been
verified to work under linux and MacOSX.

These scripts (though not aakbar itself) depend on the "pyfastaq" package
to chop genomes up into chunks for testing sensitivity.  Pyfastaq is available
from PYPI: https://pypi.python.org/pypi/pyfastaq.  

There are two top-level example scripts, which have no arguments:
  firmicutes9.sh - calculates signatures based on 9 firmicutes.
  strep10.sh - signatures based on 10 Streptococci.

Firmicutes are a phylum of gram-positive bacteria that have been implicated
in causing obesity.  The firmicutes9 example runs quickly (about 10 minutes on
a modern MacBook Pro), but using such a small set on a broad phylogenetic
category results in poor sensitivity (about 0.25%).

The Streptococci are a genus within the Firmicutes that include medically-
important bacteria (the "strep" of strep throat).  The sensitivity of this
set ranges between 15 and 50%, high enough to be useful for searching. This
example takes about an hour to run, mostly in searching the non-included
genome.

The three helper scripts called by the top-level example scripts are:
   genbank_downloader.sh - downloads and prepares data using curl
   calculate_signatures.sh - shows how to calculate signatures using aakbar
   split.sh - splits a genome into non-overlapping non-random chunks using fastaq
Each of these helper scripts contains help text.

Run "./firmicutes9.sh" when you are ready to see a full demo.
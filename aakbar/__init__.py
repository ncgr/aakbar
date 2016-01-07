#!/usr/bin/env python3
"""aakbar -- amino-acid k-mer signature tools
"""

# standard library imports
import sys
import os
import logging
import logging.handlers
import warnings

# external package imports
import click


# global constants
__prog__ = 'aakbar'
__version__ = '0.1'
__usage__ = '%prog COMMAND [command-specific options]'
__author__ = 'Joel Berendzen'
__email__ = 'joel@generisbio.com'
__copyright__ = """Copyright (C) 2016, GenerisBio, LLC.  All rights reserved.
"""
FILE_LOGLEVEL = logging.DEBUG
TERMINAL_LOGLEVEL = logging.INFO
REPORT_INTERVAL = 1000.
PLOT_TYPE = '.png'
DEFAULT_CUTOFF = 3
TEST_SEQ_LIST = ['ABCDEFGHIK',
               'ABCCEFGHIK',
               'ABCCCFGHIK',
               'ABCDEFGHII',
               'ABCDEFGHIA',
               'AACDEFGHIK',
               'ABCABCDEFG',
               'ABCDEFGHIII',
                "MEKLQLHVPRVKLGSQGLEISRLGFGCVGLSGLYNAPLSHEAGCSIIKEAFNMGVTFFDTSDFYGLNHDNEIMIGKALKE" +\
                "LPREKVQLATKFGLVRSDGVFAGVKGTPEYVRQCCEASLKRLDVEYIDLYYQHRVDTSVPIEDTMGELKKLVNEGKIKYI" +\
                "GLSQASPDTIKRAHAVHPISALQMEYSLWTRDIEEEIIPLCRELGIGIVAYSPLGHGFFAGKAAVETLPSQSALAEDARF" +\
                "SGENLEKNKLFYNRIADLASKHSCTPSQLALAWFLHQGNDIVPIPGTTKIKNLENNVGSVAVKLTNAELSEISDAVPVYE" +\
                "VAGEAPGLGSLSQYTWKFATTPSK"]


@click.group()
@click.version_option(version=__version__, prog_name=__prog__)
def aakbar():
    """aakbar -- amino-acid k-mer signature tools

    Outputs:
        a log file will be written to the current working directory
    """
    # set up logging
    global logger
    logger = logging.getLogger(__prog__)
    logger.setLevel(FILE_LOGLEVEL)
    stderrHandler = logging.StreamHandler(sys.stderr)
    stderrFormatter = logging.Formatter('%(message)s')
    stderrHandler.setFormatter(stderrFormatter)
    stderrHandler.setLevel(TERMINAL_LOGLEVEL)
    logger.addHandler(stderrHandler)

    memoryHandler = logging.handlers.MemoryHandler(100, flushLevel = 50)
    memoryHandler.setLevel(FILE_LOGLEVEL)
    logger.addHandler(memoryHandler)

    operation = 'test'
    # configure logging to file in present working directory
    logfileName = __prog__ + '-'+ operation + '.log'
    logfilePath = os.path.join(os.path.realpath('.'), logfileName)
    try:
        os.remove(logfilePath)
    except OSError:
        pass
    logfileHandler = logging.FileHandler(logfilePath)
    logfileFormatter = logging.Formatter('%(levelname)s: %(message)s')
    logfileHandler.setFormatter(logfileFormatter)
    logfileHandler.setLevel(FILE_LOGLEVEL)
    memoryHandler.setTarget(logfileHandler)
    memoryHandler.flush()
    memoryHandler.close()
    logger.removeHandler(memoryHandler)
    logger.addHandler(logfileHandler)

    # handle runtime warnings (e.g., from pandas) as exceptions
    warnings.filterwarnings('error')


@aakbar.command()
@click.option('--cutoff', default=DEFAULT_CUTOFF, help='minimum cutoff value')
@click.pass_context
def test_context(ctx, cutoff):
    print('cutoff =%d' %cutoff)
    print('context=', ctx)
    print(dir(ctx))
    print(ctx.info_name)
    logger.info('testing with mask cutoff = %d' %cutoff)
    n_seq = 0
    for seq in TEST_SEQ_LIST:
        logger.info('  %d: %s' %(n_seq, seq))
        n_seq += 1


if __name__ == '__main__':
    aakbar()
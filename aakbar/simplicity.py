# -*- coding: utf-8 -*-
'''Simplicity masking and scoring classes.
'''

# external packages
import numpy as np
import pandas as pd

# module imports
from .common import *
from . import cli, get_user_context_obj, logger

class SimplicityObject(object):
    '''Define the interfaces needed to make simplicity calculations.

    '''
    def __init__(self, default_cutoff=DEFAULT_SIMPLICITY_CUTOFF):
        self.cutoff = default_cutoff
        self.k = None
        self.label = 'null'
        self.desc = 'no simplicity calculation'
        self.testcases = [('     non-repeated', 'ABCDEFGHIJ'),
                          ('Simple Repeated Letters', ''),
                          ('15 in a row', 'AAAAAAAAAAAAAAA')
                          ]


    def set_cutoff(self, cutoff):
        '''Set the lower bound on simplicity.

        :param cutoff:  Cutoff value, lower results in more sequence masked.
        :return: None
        '''
        if cutoff < 2:
            logger.error('Simplicity cutoff must be >=2.')
            sys.exit(1)
        else:
            self.cutoff = cutoff


    def set_k(self, k):
        '''Set the k value over which rolling scores are summed.

        :param k: Number of characters in window.
        :return:  None.
        '''
        if k < 2:
            logger.error('k must be >=2')
            sys.exit(1)
        else:
            self.k = k


    def mask(self,seq):
        '''Returns the input string.

        :param seq: Input string.
        :return: Masked copy of input string.  Must be the same length.
        '''
        return seq


    def score(self, seq):
        '''Count the number masked (by lower-case) over a window.

        :param seq: String of characters.
        :param window_size: Size of window over which to calculate.
        :return: Integer array of counts.
        '''
        is_lower = np.array([char.islower() for char in to_str(seq)], dtype=int)
        rolling_count = pd.rolling_sum(is_lower,window=self.k).astype(int)
        return rolling_count[self.k-1:]


class RunlengthSimplicity(SimplicityObject):
    '''Define simplicity by the number of repeated letters.

    '''
    def __init__(self, default_cutoff=DEFAULT_SIMPLICITY_CUTOFF):
        super().__init__(default_cutoff=default_cutoff)
        self.label = 'runlength'
        self.desc = 'runlength (repeated characters)'
        self.testcases += [
                 ('Simple Repeated Letters', ''),
                 ('double, beginning', 'AABCDEFGHI'),
                 (' double in middle', 'ABCDEEFGHI'),
                 ('double at end', 'ABCDEFGHII'),
                 ('double everywhere', 'AABCDDEFGG'),
                 ('triple, beginning', 'AAABCDEFGH'),
                 (' triple in middle', 'ABCDEEEFGH'),
                 ('triple at end', 'ABCDEFGIII'),
                 ('quad in middle', 'ABCDEEEEFG'),
                 ('longer string with insert',
                  'AAAAAAABCDEFGHIJKLLLLLLL'),
                 #
                 ('Mixed Patterns', ''),
                 ('mixed repeat', 'BCAABCABCA')
        ]


    def _runlength(self, s):
        return [all([s[i+j+1] == s[i] for j in range(self.cutoff-1)])
                for i in range(len(s)-self.cutoff+1)]


    def mask(self,seq):
        '''Mask high-simplicity positions in a string.

        :param s: Input string.
        :return: Input string with masked positions changed to lower-case.
        '''
        for pos in [i for i, masked in
                    enumerate(self._runlength(to_str(seq).upper()))
                    if masked]:
            if isinstance(seq, str): #strings need to have whole length set
                seq = seq[:pos] +seq[pos:pos+self.cutoff].lower() + seq[pos+self.cutoff:]
            else:
                seq[pos:pos+self.cutoff] = to_str(seq[pos:pos+self.cutoff]).lower()
        return seq


#
# Instantiations of classes.
#
NULL_SIMPLICITY = SimplicityObject()
RUNLENGTH_SIMPLICITY = RunlengthSimplicity()

@cli.command()
@click.option('--cutoff', default=DEFAULT_SIMPLICITY_CUTOFF, show_default=True,
              help='Maximum simplicity to keep.')
@click.option('-k', default=DEFAULT_K, show_default=True,
              help='k-mer size for score calculation.')
def demo_simplicity(cutoff, k):
    '''Demo self-provided simplicity outputs.

    :param cutoff: Simplicity value cutoff, lower is less complex.
    :param window_size: Window size for masking computation..
    :return:
    '''
    user_ctx = get_user_context_obj()
    simplicity_obj = user_ctx['simplicity_object']
    simplicity_obj.set_cutoff(cutoff)
    logger.info('Simplicity function is %s with cutoff of %d.',
                simplicity_obj.desc, cutoff)
    simplicity_obj.set_k(k)
    logger.info('           Mask window demo for %d-mers:', k)
    mask_test = 'AAAAAAAAAAaaaaAAAAAAAAAAAaaaaaAAAAAAAAAAAAaaaaaaAAAAAAAAAAAAAaaaaaaaAAAAAAAAAAAAAA'
    logger.info('      in: %s\n S-score: %s\n', mask_test,
                ''.join(['%X'%i for i in
                         simplicity_obj.score(mask_test)]))
    for desc, case in simplicity_obj.testcases:
        if case is '':
            logger.info('              %s', desc)
        else:
            masked_str = simplicity_obj.mask(case)
            logger.info('%s:\n      in: %s\n     out: %s\n S-score: %s\n',
                        desc,
                        case,
                        masked_str,
                        ''.join(['%X'%i for i in
                                 simplicity_obj.score(masked_str)
                                 ]
                                )
                        )


# -*- coding: utf-8 -*-
'''Simplicity masking and scoring classes.
'''
import os
import shutil
# 3rd-party packages
import pyfaidx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# module packages
from . import cli
from .common import *
#
# global constants
#
TERM_CHAR = '$'
NUM_HISTOGRAM_BINS = 25
DEFAULT_SMOOTH_WINDOW = 11 # best if this is odd
WINDOW_TYPES = ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
BWT_ONLY = False
DENSITY_CUTOFF = 0.1
#
# class definitions
#


class RunlengthSimplicity(SimplicityObject):
    '''Define simplicity by the number of repeated letters.

    '''

    def __init__(self, default_cutoff=DEFAULT_SIMPLICITY_CUTOFF):
        super().__init__(default_cutoff=default_cutoff)
        self.label = 'runlength'
        self.desc = 'runlength (repeated characters)'

    def _runlength(self, s):
        return [all([s[i + j + 1] == s[i] for j in range(self.cutoff - 1)])
                for i in range(len(s) - self.cutoff + 1)]

    def mask(self, seq, print_results=False):
        '''Mask high-simplicity positions in a string.

        :param s: Input string.
        :return: Input string with masked positions changed to lower-case.
        '''
        for pos in [i for i, masked in
                    enumerate(self._runlength(to_str(seq).upper()))
                    if masked]:
            if isinstance(seq, str):  # strings need to have whole length set
                seq = seq[:pos] + seq[pos:pos + self.cutoff].lower() + \
                    seq[pos + self.cutoff:]
            else:
                seq[pos:pos +
                    self.cutoff] = to_str(seq[pos:pos +
                                              self.cutoff]).lower()
        return super().mask(seq)


class LetterFrequencySimplicity(SimplicityObject):
    '''Define simplicity by the number of repeated letters.
    '''

    def __init__(self,
                 default_cutoff=DEFAULT_SIMPLICITY_CUTOFF,
                 window_size=None):
        global config_obj
        super().__init__(default_cutoff=default_cutoff)
        if window_size is None:
            try:
                self.window_size = config_obj.config_dict['simplicity_window']
            except KeyError:
                self.window_size = DEFAULT_SIMPLICITY_WINDOW
        else:
            self.window_size = window_size
        self.label = 'letterfreq%d' % self.window_size
        self.desc = 'letter frequency in window of %d residues' % self.window_size


    def mask(self, seq, print_results=False):
        '''Mask high-simplicity positions in a string.

        :param s: Input string.
        :return: Input string with masked positions changed to lower-case.
        '''
        out_str = to_str(seq)
        end_idx = len(out_str) - 1
        byte_arr = np.array([char for char in to_bytes(out_str.upper())])
        mask_positions = set()
        #
        # test character by character for number of occurrances over window
        #
        for char in set(byte_arr):  # test character by character
            char_positions = list(np.where(byte_arr == char)[0])
            while len(char_positions) >= self.cutoff:
                testpos = char_positions.pop(0)
                next_positions = char_positions[:self.cutoff - 1]
                if next_positions[-1] - testpos < self.window_size:
                    mask_positions = mask_positions.union(
                        set([testpos] + next_positions))
        #
        # mask everything
        #
        for pos in mask_positions:
            out_str = out_str[:pos] + out_str[pos].lower() + out_str[pos + 1:]

        if isinstance(seq, str):  # strings need to have whole length set
            seq = out_str
        else:  # may be MutableSeq that needs lengths
            seq[:end_idx] = out_str[:end_idx]
        return super().mask(seq)


class GenerisSimplicity(SimplicityObject):
    '''Define simplicity by the number of repeated letters.

    '''

    def __init__(self,
                 default_cutoff=DEFAULT_SIMPLICITY_CUTOFF,
                 window_size=None):
        super().__init__(default_cutoff=default_cutoff)
        if window_size is None:
            try:
                self.window_size = config_obj.config_dict['simplicity_window']
            except KeyError:
                self.window_size = DEFAULT_SIMPLICITY_WINDOW
        else:
            self.window_size = window_size
        self.label = 'generis%d' % self.window_size
        self.desc = 'pattern by BW xform in window of %d residues' % self.window_size


    def _runlength(self, s):
        return [all([s[i + j + 1] == s[i] for j in range(self.cutoff - 1)])
                for i in range(len(s) - self.cutoff + 1)]

    def _bwt(self, s):
        '''Burrows-Wheeler Transform.

        :param s: Input string.  Must not contain TERMCHAR.
        :return: Transformed string.
        '''
        s = s + TERM_CHAR
        return ''.join([x[-1] for x in sorted([s[i:] + s[:i]
                                               for i in range(len(s))])])

    def _ibwt(self, s):
        '''Inverse Burrows-Wheeler Transform on uppercase-only comparisons.

        :param s: Transformed string with mixed upper and lower case.
        :return: Untransformed string with original order.
        '''
        L = [''] * len(s)
        for i in range(len(s)):
            L = sorted([s[i] + L[i] for i in range(len(s))],
                       key=str.upper)
        return [x for x in L if x.endswith(TERM_CHAR)][0][:-1]

    def merge_mask_regions(self, mask, max_run):
        "merge regions separated by less than max_run"
        runs = self.run_lengths(mask)
        for i in range(len(runs)):
            if runs[i] <= max_run:
                mask[i] = True
        return mask

    def unset_small_regions(self, mask, min_run):
        "merge regions separated by less than max_run"
        runs = self.run_lengths([int(not i) for i in mask])
        for i in range(len(runs)):
            if mask[i] and (runs[i] < min_run-1):
                mask[i] = False
        return mask


    def mask(self, seq, print_results=False):
        '''Mask high-simplicity positions in a string.

        :param s: Input string, will be converted to all-uppercase.
        :return: Input string with masked positions changed to lower-case.
        '''
        out_str = to_str(seq)
        end_idx = len(out_str) - 1
        upper_str = out_str.upper()
        # run-length mask in direct space
        if not BWT_ONLY:
            dir_mask = self._runlength(upper_str)
            for pos in [i for i, masked in
                        enumerate(dir_mask)
                        if masked]:
                out_str = out_str[:pos] + out_str[pos:pos +\
                          self.cutoff].lower() +\
                          out_str[pos + self.cutoff:]
            if print_results:
                logger.info('     rlm: %s', colorize_string(out_str))
        #run-length mask in Burrows-Wheeler space
        bwts = self._bwt(upper_str)
        bwt_mask = self._runlength(bwts)
        #bwt_mask = self.merge_mask_regions(bwt_mask, 2)
        #bwt_mask = self.unset_small_regions(bwt_mask, self.cutoff)
        for pos in [i for i, masked in
                    enumerate(bwt_mask)
                    if masked]:
            bwts = bwts[:pos] + bwts[pos:pos + \
                self.cutoff].lower() + bwts[pos + self.cutoff:]
        ibwts = self._ibwt(bwts)
        if print_results:
            logger.info('      bwt: %s', colorize_string(bwts))
            logger.info(' inv. bwt: %s', colorize_string(ibwts))
            logger.info('runlength: %s', colorize_string(out_str))
        # add in mask from inverse-transformed string
        for pos in [i for i, char in
                    enumerate(ibwts) if char.islower()]:
            out_str = out_str[:pos] + out_str[pos].lower() + out_str[pos + 1:]
        if print_results:
            logger.info('  generis: %s', colorize_string(out_str))
        if isinstance(seq, str):  # strings need to have whole length set
            seq = out_str
        else:                    # may be MutableSeq that needs lengths
            seq[:end_idx] = out_str[:end_idx]
        return  super().mask(seq)

#
# Instantiations of classes.
#
NULL_SIMPLICITY = SimplicityObject()
RUNLENGTH_SIMPLICITY = RunlengthSimplicity()
LETTERFREQ_SIMPLICITY = LetterFrequencySimplicity()
GENERIS_SIMPLICITY = GenerisSimplicity()

@cli.command()
@click.option('--smooth/--no-smooth', default=True, is_flag=True,
              help='Finish with real-space smoothing.')
@click.option('--cutoff', default=DEFAULT_SIMPLICITY_CUTOFF, show_default=True,
              help='Maximum simplicity to keep.')
@click.option('-k', default=DEFAULT_K, show_default=True,
              help='k-mer size for score calculation.')
def demo_simplicity(smooth, cutoff, k):
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
    simplicity_obj.use_smoother(smooth)
    logger.info('           Mask window demo for %d-mers:', k)
    mask_test = 'AAAAAAAAAAaaaaAAAAAAAAAAAaaaaaAAAAAAAAAAAAaaaaaaAAAAAAAAAAAAAaaaaaaaAAAAAAAAAAAAAA'
    logger.info('      in: %s\n S-score: %s\n', mask_test,
                ''.join(['%X' % i for i in
                         simplicity_obj.score(mask_test)]))
    try:
        window_size = user_ctx['simplicity_window']
    except KeyError:
        window_size = DEFAULT_SIMPLICITY_WINDOW
    smoother_tests = [('L end', 'aAAAAAAAAAAAAA'),
                      ('R end', 'AAAAAAAAAAAAAa'),
                      ('singleton', 'AAAAAaAAAAAAAAaaaaAAAAA'),
                      ('non-windowed singleton', 'AaAAAAaAAAAaAAAAAaAAAAAaAAAAAA')]
    logger.info('Demo smoother with window size of %d' %simplicity_obj.k)
    for label, pattern in smoother_tests:
        logger.info('%s:  %s', label, colorize_string(pattern))
        smoothed = simplicity_obj.smoother(pattern)
        logger.info('smoothed: %s', colorize_string(smoothed))
    for desc, case in simplicity_obj.testcases:
        if case is '':
            logger.info('              %s', desc)
        else:
            logger.info('\n%s:', desc)
            logger.info('      in: %s', case)
            masked_str = simplicity_obj.mask(case, print_results=True)
            logger.info(' smoothed: %s', colorize_string(masked_str))

def num_masked(seq):
    """Count the number of lower-case characters in a sequence.

        :param seq: Sequence of characters.
        :type seq: str, bytes, or other convertible sequence type.
        :return: Count of lower-case characters.
        :rtype: int
    """
    gene = to_str(seq)
    mask = []
    [mask.append(gene[i].islower()) for i in range(len(gene))]
    masked = sum(mask)
    return masked

class Smoother(object):
    '''Smooth data using a window function'''
    def __init__(self,
                 window_len=DEFAULT_SMOOTH_WINDOW,
                 window_type='flat'):
        if window_len < 3:
            raise ValueError('Window length must be >3')
        self.window_len = window_len
        if window_type not in WINDOW_TYPES:
            raise ValueError('Window type "%s" unknown'%window_type)
        if window_type == 'flat': #moving average
            self.window = np.ones(window_len,'d')
        else:
            self.window = eval('np.' + window_type + '(window_len)')
        self.window = self.window/self.window.sum() # normalize

    def smooth(self, x, reflect=False):
        """convolve the window with the signal"""
        if reflect:
            signal = np.r_[x[self.window_len-1:0:-1],
                           x,
                           x[-2:-self.window_len-1:-1]]
        else:
            signal = x
        return np.convolve(self.window, signal ,mode='valid')


@cli.command()
@click.option('--cutoff',
              default=DEFAULT_SIMPLICITY_CUTOFF,
              help='Minimum simplicity level to unmask.')
@click.option('--smooth/--no-smooth',
              default=True,
              help='Smooth mask profile over window.')
@click.argument('infilename', type=str)
@click.argument('outfilestem', type=str)
@click.argument('setlist', nargs=-1, type=DATA_SET_VALIDATOR)
def peptide_simplicity_mask(cutoff, smooth, infilename, outfilestem, setlist):
    '''Lower-case high-simplicity regions in FASTA.

    :param infilename: Name of input FASTA files for every directory in setlist.
    :param outfilestem: Stem of output filenames.
    :param cutoff: Minimum simplicity level to unmask.
    :param plot: If specified, make a histogram of masked fraction.
    :param setlist: List of defined sets to iterate over.
    :return:

    Note that this calculation is single-threaded and may be time-consuming, so
    starting multiple processes may be a good idea.
    '''
    global config_obj
    user_ctx = get_user_context_obj()
    setlist = DATA_SET_VALIDATOR.multiple_or_empty_set(setlist)
    simplicity_obj = user_ctx['simplicity_object']
    simplicity_obj.set_cutoff(cutoff)
    simplicity_obj.use_smoother(smooth)
    logger.info('Simplicity function is %s with cutoff of %d.',
                simplicity_obj.desc, cutoff)
    if simplicity_obj.smooth:
        logger.info('Mask will be smoothed over window of %d residues',
                    simplicity_obj.k)
    logger.debug('Reading from FASTA file "%s".', infilename)
    instem, ext = os.path.splitext(infilename)
    outfilename = outfilestem + ext
    logger.debug('Output FASTA file name is "%s".', outfilename)
    for calc_set in setlist:
        dir = config_obj.config_dict[calc_set]['dir']
        inpath = os.path.join(dir, infilename)
        outpath = os.path.join(dir, outfilename)
        shutil.copy(inpath, outpath)
        fasta = pyfaidx.Fasta(outpath, mutable=True)
        if user_ctx['first_n']:
            keys = list(fasta.keys())[:user_ctx['first_n']]
        else:
            keys = fasta.keys()
        if user_ctx['progress']:
            with click.progressbar(keys, label='%s genes processed' % calc_set,
                                   length=len(keys)) as bar:
                for key in bar:
                    masked_gene = simplicity_obj.mask(fasta[key])
        else:
            for key in keys:
                masked_gene = simplicity_obj.mask(fasta[key])
        fasta.close()


@cli.command()
@click.option('--window',
              default=21,
              help='Window size for calculations.')
@click.option('--histmax',
              default=10.,
              help='Maximum value for histogram.')
@click.option('--outname', '-o',
              default='marker_histogram.png',
              help='Maximum value for histogram.')
@click.argument('filelist', nargs=-1)
def plot_mask_stats(window, histmax, outname, filelist):
    '''Compute stats on masked sequences'''
    if not len(filelist):
        logger.error('filelist is required')
        sys.exit(1)
    logger.info('Computing marker stats over window of size %d', window)
    global config_obj
    first_n = get_user_context_obj()['first_n']
    mask_dict = {}
    smoother = Smoother(window_len=window,
                        window_type='hanning')
    namelist = []
    for filename in filelist:
        filepath = Path(filename)
        name = os.path.splitext(filepath.name)[0]
        namelist.append(name)
        fasta = pyfaidx.Fasta(str(filepath))
        smooth_list = []
        if first_n:
            keys = list(fasta.keys())[:first_n]
        else:
            keys = fasta.keys()
        logger.info('  Read %d sequences from file "%s".',
                    len(keys), filename)
        for key in keys:
            seq = to_str(fasta[key])
            if len(seq) > window:
                mask = np.array([int(c.islower()) for c in seq])
                smooth_list.append(smoother.smooth(mask))
        mask_dict[name] = {}
        densities =  np.concatenate(smooth_list)
        above_cutoff = (densities > DENSITY_CUTOFF).sum()*100/len(densities)
        print("%s: %.2f%% > %f"%(name, above_cutoff, DENSITY_CUTOFF))
        mask_dict[name]['Density'] = densities
        mask_dict[name]['cutoff'] = above_cutoff
    firstname = namelist[0]
    namelist.pop(0)
    ax = sns.distplot(mask_dict[firstname]['Density'],
                      bins=100,
                      hist=True,
                      kde_kws={'label': '%s (%d%%>%.1f)' % (firstname,
                                    int(mask_dict[firstname]['cutoff']+0.5),
                                                            DENSITY_CUTOFF)},
                      norm_hist=True,
                      rug=False,
                      )
    ax.set_xlabel('Mask Density')
    ax.set_ylabel('Percent in Bin')
    ax.set_xlim([0.01, 1.0])
    ax.set_ylim([0, histmax])
    for name in namelist:
        sns.distplot(mask_dict[name]['Density'],
                     bins=100,
                     hist=True,
                     kde_kws={'label': '%s (%d%%>%.1f)' %(name,
                                        int(mask_dict[name]['cutoff']+0.5),
                                                          DENSITY_CUTOFF)},
                     rug=False,
                     norm_hist=True,
                     ax=ax)
    plt.title('Mask Histogram + KDE Over Window of %d'%window)
    logger.info('Saving plot to %s', outname)
    #plt.yscale('log')
    plt.savefig(outname, dpi=200)
    #plt.show()

    #percent_masked = 100. * \
    #                 num_masked(masked_gene) / len(masked_gene)
    #percent_masked_list.append(percent_masked)
    #
    # histogram masked regions
    #
    #(hist, bins) = np.histogram(percent_masked_list,
    #                            bins=np.arange(0., 100., 100. / NUM_HISTOGRAM_BINS))
    #bin_centers = (bins[:-1] + bins[1:]) / 2.
    #hist = hist * 100. / len(percent_masked_list)
    #hist_filepath = os.path.join(dir, histfilename)
    #logger.debug('writing histogram to file "%s".', hist_filepath)
    #pd.Series(hist, index=bin_centers).to_csv(hist_filepath, sep='\t',
    #                                         header=True)
    #
    # plot histogram, if requested
    #
    #plotpath = os.path.join(dir, plotname)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(bin_centers, hist)
    #plt.title(
    #    'Peptide %s Simplicity Distribution with Cutoff %d' %
    #    (simplicity_obj.label.capitalize(), cutoff))
    #plt.xlabel('Percent of Peptide Sequence Masked')
    #plt.ylabel('Percent of Peptide Sequences')
    #plt.savefig(plotpath)


@cli.command()
@click.argument('window_size')
def set_simplicity_window(window_size):
    '''Define size of simplicity window.
    '''
    global config_obj
    if window_size == ():
        try:
            window_size = config_obj.config_dict['simplicity_window']
            default = ''
        except KeyError:
            window_size = DEFAULT_SIMPLICITY_WINDOW
            default = ' (default)'
        logger.info('Window size is %d residues%s',
                    window_size, default)
    try:
        window_size = int(window_size)
    except ValueError:
        logger.error('Window size must be an integer value.')
        sys.exit(1)
    if window_size < 3:
        logger.error('Window size must be >=3.')
        sys.exit(1)
    config_obj.config_dict['simplicity_window'] = window_size
    logger.info(
        'Window size for letter-frequency simiplicity calculation is now %d residues.',
        window_size)
    config_obj.write_config_dict()

@cli.command()
@click.argument('infile', type=str)
@click.argument('outfile', type=str)
def colorize_fasta(infile, outfile):
    "Color lower-case parts of sequence."
    infilepath = Path(infile)
    outfilepath = Path(outfile)
    with infilepath.open(mode='rU') as infh:
        line = infh.readline()
        with outfilepath.open(mode='wt') as outfh:
            for line in infh:
                if line.startswith('>'):
                    outfh.write(line)
                else:
                    outfh.write(colorize_string(line))

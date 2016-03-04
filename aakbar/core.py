# -*- coding: utf-8 -*-
'''Core commands for aakbar
'''

# standard library imports
import os
import shutil

# external packages
import click
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyfaidx
import bitarray

# module imports
from .common import *
from . import cli, get_user_context_obj, logger, log_elapsed_time

# global definitions
DEFAULT_K = 10
ALPHABETSIZE = 20 # number of amino acids
UNITNAME = 'Mbasepair'
UNITMULTIPLIER = 3.E6
AMBIGUOUS_RESIDUE = 'X'
DEFAULT_CUTOFF = 3
NUM_HISTOGRAM_BINS = 25

class DataSetValidator(click.ParamType):
    '''Validate that set names are defined.
    '''
    global config_obj
    name = 'set'
    def convert(self, argset, param, ctx):
        '''Verify that arguments refers to a valid set in the configuration dictionary.

        :param argset:
        :param param:
        :param ctx:
        :return:
        '''
        if argset not in config_obj.config_dict['sets'] and argset != 'all':
            logger.error('"%s" is not a recognized set', argset)
            sys.exit(1)
        return argset
DATA_SET_VALIDATOR = DataSetValidator()

def if_all_valid(seq):
    '''Check to see if all characters (residues) are unambiguous and non-masked.

    :param seq: Input sequence, possibly containing ambiguous positions ('X')
                and masked (lower-case) positions.
    :return: True if no ambiguous or masked positions are contained in seq.
    '''
    if any(filter(str.islower, seq)) or AMBIGUOUS_RESIDUE in seq:
        return False
    else:
        return True

def num_masked(seq):
    """Count the number of lower-case characters in a sequence.

        :param seq: Sequence of characters.
        :type seq: str, bytes, or other convertible sequence type.
        :return: Count of lower-case characters.
        :rtype: int
    """
    gene = to_str(seq)
    mask = bitarray.bitarray()
    [mask.append(gene[i].islower()) for i in range(len(gene))]
    masked = sum(mask)
    return masked

class RunlengthComplexity(object):
    '''Define complexity by the number of repeated letters.

    '''
    def __init__(self, default_cutoff=3):
        self.cutoff = default_cutoff
        self.label = 'runlength'
        self.desc = 'runlength (repeated characters)'
        self.__name__ = 'ComplexityObject'
        self.testcases = (('     non-repeated', 'ABCDEFGHIJ'),
                 #
                 ('Simple Repeated Letters', ''),
                 ('double, beginning', 'AABCDEFGHI'),
                 (' double in middle', 'ABCDEEFGHI'),
                 ('    double at end', 'ABCDEFGHII'),
                 ('double everywhere', 'AABCDDEFGG'),
                 ('triple, beginning', 'AAABCDEFGH'),
                 (' triple in middle', 'ABCDEEEFGH'),
                 ('    triple at end', 'ABCDEFGIII'),
                 ('   quad in middle', 'ABCDEEEEFG'),
                 #
                 ('Mixed Patterns', ''),
                 ('     mixed repeat', 'BCAABCABCA')
    )

    def set_cutoff(self, cutoff):
        if cutoff < 2:
            logger.error('Cutoff must be <=2.')
            sys.exit(1)
        else:
            self.cutoff = cutoff

    def _runlength(self, s):
        return [all([s[i+j+1] == s[i] for j in range(self.cutoff-1)])
                for i in range(len(s)-self.cutoff+1)]

    def below_cutoff(self, s):
        '''Check to see if string is low-complexity.

        :param s: Input string.
        :return: True if the complexity of s is lower than the cutoff.
        '''
        return any(self._runlength(s))


    def mask(self,seq):
        '''Mask low-complexity positions in a string.

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

RUNLENGTH_COMPLEXITY = RunlengthComplexity()

def calculate_frequency_histogram(freqs, filepath):
    '''Writes frequency histograms to a khmer-compatible file.

    :param freqs: Vector of frequencies
    :param filepath: Output file path.
    :return:
    '''
    # calculate frequency histogram
    binvals, freq_hist = np.unique(freqs, return_counts=True)
    max_freq = max(binvals)
    logger.debug('writing frequency histogram to %s', filepath)
    cumulative = np.cumsum(freq_hist)
    total = np.sum(freq_hist)
    pd.DataFrame({'abundance':binvals,
                  'count':freq_hist,
                  'cumulative':cumulative,
                  'cumulative_fraction':cumulative/total},
                  columns=('abundance',
                           'count',
                           'cumulative',
                           'cumulative_fraction')).to_csv(filepath,
                                                         index=False,
                                                         float_format='%.3f')
    logger.info('Maximum k-mer frequency is %d (%.1e/%s)',
                max_freq,
                max_freq*UNITMULTIPLIER/len(freqs),
                UNITNAME)

@cli.command()
@click.option('-k', default=DEFAULT_K, show_default=True, help='k-mer length')
@click.argument('infilename', type=str)
@click.argument('outfilestem', type=str)
@click.argument('setlist', nargs=-1, type=DATA_SET_VALIDATOR)
@log_elapsed_time()
def calculate_peptide_kmers(k, infilename, outfilestem, setlist):
    '''calculate peptide k-mers and histogram.


    '''
    user_ctx = get_user_context_obj()
    try:
        first_n = user_ctx['first_n']
        progress = user_ctx['progress']
    except KeyError:
        logger.debug('first_n and progress not found in user context.')
        first_n = None
        progress = None
    logger.debug('k-mer size is %d' %k)
    if first_n:
        logger.debug('Only first %d records will be used', first_n)
    if setlist[0] == 'all':
        setlist = list(config_obj.config_dict['sets'])
    if setlist == []:
        logger.error('Empty setlist, make sure sets are defined.')
        sys.exit(1)
    logger.info('Input file name is "%s".', infilename)
    logger.info('Output file stem is "%s".', outfilestem)
    for calc_set in setlist:
        logger.info('Doing set "%s"', calc_set)
        dir = config_obj.config_dict[calc_set]['dir']
        infilepath = os.path.join(dir, infilename)
        if not os.path.exists(infilepath):
            logger.error('input file "%s" does not exist.', infilepath)
            sys.exit(1)

        # read the first time just to get sizes
        logger.debug('Reading file %s', infilepath)
        nresidues = 0
        nkmers = 0
        nrecs = 0
        kmer_list_ends = []
        for record in SeqIO.parse(str(infilepath),'fasta'):
            seq = record.seq.rstrip('*')
            nresidues += len(seq)
            nkmers += len(seq) - k + 1
            nrecs += 1
            kmer_list_ends.append(nkmers-1)
            if nrecs == first_n:
                    break
        logger.info('%d records and %d residues will make %d non-unique k-mers' %(nrecs,
                                                                                  nresidues,
                                                                                  nkmers))
        # read the second time to populate k-mer array
        logger.debug('generating k-mer array')
        cur_rec = 0
        kmer_list = []
        if progress:
            with click.progressbar(SeqIO.parse(infilepath,'fasta'),
                                   length=nrecs,
                                   label='%s genes' %calc_set) as bar:
                for record in bar:
                    seq = record.seq.rstrip('*')
                    kmer_list.extend([str(seq[i:i+k]) for i in range(len(seq)-k+1)
                                     if if_all_valid(str(seq[i:i+k]))])
                    cur_rec += 1
                    if cur_rec == first_n:
                        break
        else:
                for record in SeqIO.parse(infilepath,'fasta'):
                    seq = record.seq.rstrip('*')
                    kmer_list.extend([str(seq[i:i+k]) for i in range(len(seq)-k+1)
                                     if if_all_valid(str(seq[i:i+k]))])
                    cur_rec += 1
                    if cur_rec == first_n:
                        break
        kmer_arr = np.array(kmer_list, dtype=np.dtype(('S%d'%(k))))
        logger.info('Ignoring %.0f%% of k-mers because they are masked.', 100.*(1.-len(kmer_list)/nkmers))
        del kmer_list

        # calculate unique k-mers
        unique_kmers, freqs = np.unique(kmer_arr, return_counts=True)
        logger.info('%d unique k-mers (%.6f%% of %d possible %d-mers)',
                    len(unique_kmers),
                    len(unique_kmers)*100./(ALPHABETSIZE**k),
                    ALPHABETSIZE**k,
                    k)

        #write k-mers and counts
        kmer_filepath = os.path.join(dir, outfilestem+'_terms.tsv')
        logger.debug('writing unique k-mers and counts to %s', kmer_filepath)
        sort_arr = np.argsort(freqs)
        pd.DataFrame({'term':[i.decode('UTF-8') for i in unique_kmers[sort_arr]],
                      'count':freqs[sort_arr]},
                     columns=('term', 'count')).to_csv(kmer_filepath, index=False, sep='\t')
        del sort_arr

        hist_filepath = os.path.join(dir, outfilestem+'_hist.csv')
        calculate_frequency_histogram(freqs, hist_filepath)


@cli.command()
@click.option('--runlength', default=DEFAULT_CUTOFF, show_default=True,
              help='minimum runlength to remove')
@click.argument('infilename', type=str)
@click.argument('outfilestem', type=str)
@click.argument('setlist', nargs=-1, type=DATA_SET_VALIDATOR)
@log_elapsed_time()
def filter_peptide_kmers(runlength, infilename, outfilestem, setlist):
    '''Removes low-complexity peptide k-mers.

    '''
    user_ctx = get_user_context_obj()
    complexity_obj = user_ctx['complexity_object']
    try:
        first_n = user_ctx['first_n']
        progress = user_ctx['progress']
    except KeyError:
        logger.debug('first_n and progress not found in user context.')
        first_n = None
        progress = None
    if first_n:
        logger.debug('Only first %d records will be used', first_n)
    if setlist[0] == 'all':
        setlist = list(config_obj.config_dict['sets'])
    if setlist == []:
        logger.error('Empty setlist, make sure sets are defined.')
        sys.exit(1)
    logger.info('Input file name is "%s".', infilename)
    logger.info('Output file stem is "%s".', outfilestem)
    logger.info('Maximum run length remaining is %d', runlength)
    complexity_obj.set_cutoff(runlength)
    for calc_set in setlist:
        logger.info('Doing set "%s"', calc_set)
        dir = config_obj.config_dict[calc_set]['dir']
        infilepath = os.path.join(dir, infilename)
        if not os.path.exists(infilepath):
            logger.error('input file "%s" does not exist.', infilepath)
            sys.exit(1)
        term_frame = pd.DataFrame().from_csv(infilepath, sep='\t')
        initial_term_count = len(term_frame)
        droplist = []
        k = len(term_frame.index[0])
        logger.info('%d %d-mer terms initially.', initial_term_count,
                    k)
        if progress:
            with click.progressbar(term_frame.index,
                                   length=initial_term_count,
                                   label='%s terms' %calc_set) as bar:
                for term in bar:
                    if complexity_obj.below_cutoff(term):
                        droplist.append(term)
        else:
            for term in term_frame.index:
                if complexity_obj.below_cutoff(term):
                    droplist.append(term)

        term_frame.drop(droplist, inplace=True)
        logger.info('Dropped %0.1f%% of terms.', 100.*(1.-len(term_frame)/initial_term_count))
        logger.info('%d k-mers remain (%.6f%% of %d possible %d-mers)',
                    len(term_frame),
                    len(term_frame)*100./(ALPHABETSIZE**k),
                    ALPHABETSIZE**k,
                    k)

        # write k-mers and counts
        kmer_filepath = os.path.join(dir, outfilestem+'_terms.tsv')
        term_frame.to_csv(kmer_filepath, sep='\t')

        hist_filepath = os.path.join(dir, outfilestem+'_hist.csv')
        calculate_frequency_histogram(term_frame, hist_filepath)


@cli.command()
@click.option('--cutoff', default=DEFAULT_CUTOFF, show_default=True,
              help='Minimum complexity to keep.')
def test_complexity_masking(cutoff):
    '''Tests the complexity function for some examples.

    :param runlength: Runlength value.
    :param pattern: If True, count repeating patterns.
    :return:
    '''
    user_ctx = get_user_context_obj()
    complexity_obj = user_ctx['complexity_object']
    logger.info('Complexity function is %s with cutoff of %d.',
                complexity_obj.desc, cutoff)
    complexity_obj.set_cutoff(cutoff)
    logger.info('    Description       Input     Masked   Cutoff')
    for desc, case in complexity_obj.testcases:
        if case is '':
            logger.info('           %s', desc)
        else:
            logger.info('%s: %s %s %s', desc, case,
                        complexity_obj.mask(case), complexity_obj.below_cutoff(case))


@cli.command()
@click.argument('filestem', type=str)
@click.argument('setlist', nargs=-1, type=DATA_SET_VALIDATOR)
@log_elapsed_time()
def merge_peptide_kmers(filestem, setlist):
    '''Merge kmers from a list of sets, discarding any kmers that only occur in one set.

    :param filestem: input and output filename less '_terms.tsv'
    :param setlist:
    :return:
    '''
    global config_obj
    user_ctx = get_user_context_obj()
    try:
        first_n = user_ctx['first_n']
        progress = user_ctx['progress']
    except KeyError:
        logger.debug('first_n and progress not found in user context.')
        first_n = None
        progress = None
    if first_n:
        logger.debug('Only first %d records will be used', first_n)
    if setlist[0] == 'all':
        setlist = list(config_obj.config_dict['sets'])
    if setlist == []:
        logger.error('Empty setlist, make sure sets are defined.')
        sys.exit(1)
    infilename = filestem+'_terms.tsv'
    logger.info('Input File name is "%s".', infilename)

    # get the output directory
    outdir = config_obj.config_dict['summary']['dir']
    if not os.path.exists(outdir) and not os.path.isdir(outdir):
        logger.info('Making summary directory %s', outdir)
        os.makedirs(outdir)

    # read and merge the kmer_lists
    merged_frame = None
    n_terms_in = 0
    n_sets = len(setlist)
    logger.info('Intersecting terms from %d sets.', n_sets)
    for calc_set in setlist:
        dir = config_obj.config_dict[calc_set]['dir']
        infilepath = os.path.join(dir, infilename)
        if not os.path.exists(infilepath):
            logger.error('input file "%s" does not exist.', infilepath)
            sys.exit(1)
        working_frame = pd.DataFrame.from_csv(infilepath,
                                              sep='\t', index_col=0)
        n_terms_in += len(working_frame)
        logger.info('%s: %d k-mers.' %(calc_set,len(working_frame)))
        if merged_frame is None:
            working_frame.columns = [calc_set]
            merged_frame = working_frame
        else:
            merged_frame[calc_set] = working_frame['count']

    # calculate the number of times a term appears and max count
    n_unique_terms = len(merged_frame)
    max_count = merged_frame.max(axis=1)
    intersections = n_sets - merged_frame.isnull().sum(axis=1)
    merged_frame['intersections'] = intersections
    merged_frame['max_count'] = max_count.astype(int)
    for calc_set in setlist: # delete original counts
        del merged_frame[calc_set]

    # drop terms that don't intersect in two sets
    merged_frame.drop(merged_frame[merged_frame['intersections'] == 1].index,
                      inplace=True)
    n_intersecting_terms = len(merged_frame)
    k = len(merged_frame.index[0])
    logger.info('%d shared terms (%.6f%% of possible %d-mers, %.2f%% of %d unique terms).',
                n_intersecting_terms,
                n_intersecting_terms*100./(ALPHABETSIZE**k),
                k,
                n_intersecting_terms*100./n_unique_terms,
                n_unique_terms)

    # write terms
    merged_frame.sort_values(by=['max_count', 'intersections'], inplace=True)
    kmer_filepath = os.path.join(outdir, filestem+'_terms.tsv')
    logger.debug('Writing merged k-mers and counts to "%s".', kmer_filepath)
    merged_frame.to_csv(kmer_filepath, sep='\t')

    # calculate frequency histogram to khmer-compatible histogram file
    freqs = merged_frame['max_count']
    binvals, freq_hist = np.unique(freqs, return_counts=True)
    max_freq = max(binvals)
    hist_filepath = os.path.join(outdir, filestem+'_hist.csv')
    logger.debug('Writing frequency histogram to "%s".', hist_filepath)
    cumulative = np.cumsum(freq_hist)
    total = np.sum(freq_hist)
    pd.DataFrame({'abundance':binvals,
                      'count':freq_hist,
                   'cumulative':cumulative,
                   'cumulative_fraction':cumulative/total},
                     columns=('abundance',
                              'count',
                              'cumulative',
                              'cumulative_fraction')).to_csv(hist_filepath,
                                                             index=False,
                                                             float_format='%.3f')

    # calculate histogram of intersections
    intersect_bins = []
    lastbin = 0
    nextbin = 1
    hists = {}
    sums = []
    intersect_filepath = os.path.join(outdir, filestem+'_intersect.tsv')
    logger.debug('Writing intersection frequency histograms to %s.', intersect_filepath)
    while lastbin < max_freq:
        inrange = merged_frame[merged_frame['max_count'].isin([lastbin,nextbin])]
        sums.append(len(inrange))
        ibins, ifreqs = np.unique(inrange['intersections'], return_counts=True)
        intersect_bins.append(nextbin)
        hists[nextbin] = pd.Series(ifreqs/ifreqs.sum(), index=ibins)
        lastbin = nextbin
        nextbin *= 2
    intersect_frame = pd.DataFrame(hists).transpose()
    intersect_frame['Number'] = sums
    intersect_frame.to_csv(intersect_filepath, sep='\t',
                           float_format='%.3f')

    # plot intersection histograms
    plot_type = config_obj.config_dict['plot_type']
    plot_filepath = os.path.join(outdir, filestem+'_plot.'+ plot_type)
    logger.debug('Plotting intersection histograms to %s.', plot_filepath)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x_vals = np.arange(2, n_sets+1)
    del intersect_frame['Number']
    for i in range(len(sums)):
        bin_edge = intersect_frame.index[i]
        data = intersect_frame.iloc[i]
        ax.plot(x_vals,
                data*100.,
                '-',
                label='>%d occurrances (%d)' %(bin_edge,
                                              sums[i]))
    ax.legend(loc=9)
    plt.title('Histogram of Signature Intersections Across %d Genomes' %n_sets)
    plt.xlabel('Number Intersecting')
    plt.ylabel('% of Signatures')
    #plt.show()
    plt.savefig(plot_filepath)

@cli.command()
@click.option('--cutoff', default=DEFAULT_CUTOFF, help='Minimum complexity level to unmask.')
@click.option('--plot/--no-plot', default=False, help='Plot histogram of mask fraction.')
@click.argument('infilename', type=str)
@click.argument('outfilestem', type=str)
@click.argument('setlist', nargs=-1, type=DATA_SET_VALIDATOR)
def protein_complexity_mask(cutoff, plot, infilename, outfilestem, setlist):
    '''Mask low-complexity regions by changing them to lower case

    :param infilename: Name of input files for every directory in setlist.
    :param outfilestem: Stem of output filenames.
    :param cutoff: Minimum complexity level to unmask.
    :param plot: If specified, make a histogram of masked fraction.
    :return:
    '''
    global config_obj
    user_ctx = get_user_context_obj()
    complexity_obj = user_ctx['complexity_object']
    try:
        first_n = user_ctx['first_n']
        progress = user_ctx['progress']
    except KeyError:
        logger.warn('first_n and progress not found in user context.')
        first_n = None
        progress = False
    if first_n:
        logger.info('Only first %d records will be used', first_n)
    if setlist[0] == 'all':
        setlist = list(config_obj.config_dict['sets'])
    if setlist == []:
        logger.error('Empty setlist, make sure sets are defined.')
        sys.exit(1)
    logger.info('Complexity function is %s with cutoff of %d.',
                complexity_obj.desc, cutoff)
    complexity_obj.set_cutoff(cutoff)
    logger.debug('Reading from FASTA file "%s".', infilename)
    instem, ext = os.path.splitext(infilename)
    outfilename = outfilestem + ext
    logger.debug('Output FASTA file name is "%s".', outfilename)
    histfilename = outfilestem + '-hist.tsv'
    logger.debug('Output histogram file is "%s".', histfilename)
    if plot:
        plotname = outfilestem + '.' + config_obj.config_dict['plot_type']
        logger.debug('Plot to file "%s".', plotname)
    for calc_set in setlist:
        dir = config_obj.config_dict[calc_set]['dir']
        inpath = os.path.join(dir, infilename)
        outpath = os.path.join(dir, outfilename)
        shutil.copy(inpath, outpath)
        fasta  = pyfaidx.Fasta(outpath, mutable=True)
        percent_masked_list = []
        n_genes = 0
        if first_n:
            keys = list(fasta.keys())[:first_n]
        else:
            keys = fasta.keys()
        if progress:
            with click.progressbar(fasta, label='%s genes processed' %calc_set,
                                   length=len(keys)) as bar:
                for gene in bar:
                    masked_gene = complexity_obj.mask(gene)
                    percent_masked = 100.*num_masked(masked_gene)/len(gene)
                    percent_masked_list.append(percent_masked)
                    n_genes += 1
        else:
            for gene in fasta:
                gene = complexity_obj.mask(gene)
                percent_masked = 100.*num_masked(gene)/len(gene)
                percent_masked_list.append(percent_masked)
                n_genes += 1
        fasta.close()

        # make histogram
        (hist, bins) = np.histogram(percent_masked_list, bins=np.arange(0.,100.,100./NUM_HISTOGRAM_BINS))
        bin_centers = (bins[:-1] + bins[1:])/2.
        hist = hist*100./len(percent_masked_list)
        hist_filepath = os.path.join(dir, histfilename)
        logger.debug('writing histogram to file "%s".', hist_filepath)
        pd.Series(hist, index=bin_centers).to_csv(hist_filepath, sep='\t',
                                                  float_format='%.3f')

        # plot histogram, if requested
        if plot:
            plotpath = os.path.join(dir, plotname)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(bin_centers, hist)
            plt.title('Protein %s Complexity Distribution with Cutoff %d'%(complexity_obj.label.capitalize()
                                                                           ,cutoff))
            plt.xlabel('Percent of Protein Sequence Masked')
            plt.ylabel('Percent of Protein Sequences')
            plt.savefig(plotpath)

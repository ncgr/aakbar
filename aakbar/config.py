# -*- coding: utf-8 -*-
'''Implements commands related to configuration parameters
'''

# third-party imports
import matplotlib.pyplot as plt

# module imports
from . import cli, get_user_context_obj, logger
from .common import *

# private context function
_ctx = click.get_current_context

@cli.command()
def show_config():
    '''Prints location and contents of configuration file

        Example:
            aakbar -v show_config
    '''
    global  config_obj

    if config_obj.config_dict == {}:
        logger.info('No configuration file was found.')
    else:
        logger.info('Configuration file path is "%s".', config_obj.path)
        for key in config_obj.config_dict.keys():
            logger.info('  %s: %s', key, config_obj.config_dict[key])


@cli.command()
@click.option('-f', '--force', is_flag=True, default=False,
              help='Overwrite existing definition.')
@click.argument('identifier', type=str)
@click.argument('dir', type=click.Path(writable=True))
def define_set(identifier, dir, force):
    '''Define an identifier and directory for a set.

    IDENTIFIER:
        a short code by which this data set will be known.
    DIR:
        absolute or relative path to files.  Must be writable.

    Example::
        aakbar define_set glymax ./Glycine_max

    '''
    global config_obj
    if identifier in config_obj.config_dict['sets'] and not force:
        logger.error('Set "%s" is already defined, use "--force" to override.',
                     identifier)
    else:
        path = dir
        if identifier in config_obj.config_dict['sets'] and force:
            logger.info('Overriding existing set "%s" definition.', identifier)
            [config_obj.config_dict['sets'].remove(identifier)
             for i in range(config_obj.config_dict['sets'].count(identifier))]
        config_obj.config_dict['sets'].append(identifier)
        config_obj.config_dict[identifier] = {'dir': str(path), 'label': identifier}
        config_obj.write_config_dict()
        logger.info('Set "%s" will refer to directory "%s".', identifier, dir)


@cli.command()
@click.argument('identifier', type=str)
@click.argument('label', type=str)
def label_set(identifier, label):
    '''Define label associated with a set.
    '''
    global config_obj
    if identifier not in config_obj.config_dict['sets']:
        logger.error('Set "%s" is not defined.', identifier)
        sys.exit(1)
    else:
        config_obj.config_dict[identifier]['label'] = label
        logger.info('Label for %s is now "%s".', identifier, label)
        config_obj.write_config_dict()


@cli.command()
@click.option('--plot_type', type=str, default=None)
def set_plot_type(plot_type=None):
    '''Define label associated with a set.
    '''
    global config_obj
    plot_types = plt.gcf().canvas.get_supported_filetypes()
    if plot_type is None:
        logger.info('Must supply a plot type.')
        logger.info('Supported plot types are:')
        for typename in plot_types:
            logger.info('   %s', typename)
    elif plot_type not in plot_types:
        logger.error('Plot type "%s" is not defined.', plot_type)
        logger.info('Supported plot types are:')
        for typename in plot_types:
            logger.info('   %s', typename)
        sys.exit(1)
    else:
        config_obj.config_dict['plot_type'] = plot_type
        logger.info('Plot type is now %s.', plot_type)
        config_obj.write_config_dict()


@cli.command()
@click.argument('dir', type=click.Path(writable=True))
@click.argument('label', type=str)
def define_summary(dir, label):
    '''Define summary directory and label.
    '''
    global config_obj
    path = str(Path(dir).expanduser())
    logger.info('Summary output will go in directory "%s".', path)
    config_obj.config_dict['summary']['dir'] = path
    config_obj.config_dict['summary']['label'] = label
    config_obj.write_config_dict()


@cli.command()
@click.argument('dir', type=str, default='')
def init_config_file(dir):
    '''Initialize a configuration file.

    :param dir: Optional directory in which to initialize the file.
    If not present, the system-dependent default application directory
    will be used.  If this argument is '.', then the current working
    directory will be used.  This argument accepts tilde expansions.
    '''
    global config_obj
    config_obj.write_config_dict(dir=dir, config_dict={})

@cli.command()
@click.option('--function', type=str, default=None)
def set_complexity_function(function):
    '''Defines the function to be used for complexity calculations.

    :param function: Function name of global scope that begins with "complexity_".
    :return:
    '''
    global config_obj
    known_complexity_objects = _ctx().obj['complexity_objects']
    if function is None:
        logger.info('Possible complexity objects:')
        for obj in known_complexity_objects:
            logger.info('   %s: "%s"', obj.label, obj.desc)
        try:
            current_complexity_object = config_obj.config_dict['complexity_object_label']
        except KeyError:
            current_complexity_object = 'undefined'
        logger.info('Current complexity object is %s.', current_complexity_object)
    else:
        for obj in known_complexity_objects:
            if obj.label == function:
                logger.info('Complexity function is now %s.', function)
                config_obj.config_dict['complexity_object_label'] = function
                config_obj.write_config_dict()
                return
        logger.error('Function "%s" is not a known complexity function.', function)
        sys.exit(1)


# -*- coding: utf-8 -*-

import sys
import os

from functools import wraps
from docopt import docopt
from schema import Schema, And, Use, SchemaError
import variants
import evaluate
from targetregions import get_target_regions_from_file, eval_target_regions

__version__ = '0.1.0'

Tools = {
    'sprites': variants.SpritesBedpeFile
}

Chromosomes = ['all'] + [str(i) for i in range(1, 23)] + ['X', 'Y']

def argparsed(func):
    @wraps(func)
    def wrapped(argv):
        args = docopt(func.__doc__, argv=argv)
        return func(args)
    return wrapped


@argparsed
def sim(args):
    """
Usage: svtools sim -r <file0> <tool1> <file1> <tool2> <file2>

Dump audio meta data of the <files>.

Options:
  -r <file0>      Give the file of true variants
    """
    print "here..."


@argparsed
def taregion(args):
    """
Usage: svtools taregion -r <file0> -c <chrom> <file1>

Dump audio meta data of the <files>.

Options:
  -r <file0>    The true variants file
  -c <chrom>    Specify a chromosome
  <file1>       The target regions file
    """

    schema = Schema({
        '-r': os.path.isfile,
        '-c': And(str, lambda s: s in Chromosomes),
        '<file1>': os.path.isfile,
        'taregion': True
    })

    try:
        args = schema.validate(args)
    except SchemaError as e:
        exit(e)

    truth = evaluate.get_variants(variants.SvsimBedpeFile, args['-r'], 5)
    taregions = get_target_regions_from_file(args['<file1>'])
    if args['<chrom>'] != 'all':
        truth = [x for x in truth if x.reference_name() == args['chrom']]
        taregions = [x for x in taregions if x.reference_name() == args['chrom']]
    eval_target_regions(taregions, truth)


@argparsed
def diff(args):
    """
Usage: svtools diff -r <file0> <tool1> <file1> <tool2> <file2>

Show the difference between results obtained by two tools.

Options:
  -r <file0>    The true variants file
  <tool1>       Tool name should be sprites, lumpy or pindel
  <file1>       The results file obtained by <tool1>
  <tool2>       Tool name should be sprites, lumpy or pindel
  <file2>       The results file obtained by <tool2>
    """

    schema = Schema({
        '-r': os.path.isfile,
        '<tool1>': And(Use(str.lower), lambda s: s in Tools.keys()),
        '<file1>': os.path.isfile,
        '<tool2>': And(Use(str.lower), lambda s: s in Tools.keys()),
        '<file2>': os.path.isfile,
        'diff': True
    })

    try:
        args = schema.validate(args)
    except SchemaError as e:
        exit(e)

    truth = evaluate.get_variants(variants.SvsimBedpeFile, args['-r'], 50)
    pred1 = evaluate.get_variants(Tools[args['<tool1>']], args['<file1>'])
    pred2 = evaluate.get_variants(Tools[args['<tool2>']], args['<file2>'])

    evaluate.compare_variants(
        truth, pred1, pred2, lambda v1, v2: v1.matched_with_both_overlaps(v2))


def help(argv):
    if len(argv) > 1:
        cmd = argv[-1]
        try:
            print(globals()[cmd].__doc__)
        except KeyError:
            exit("%r is not a svtools command. See 'svtools help'." % cmd)
    else:
        docopt(main.__doc__, argv='-h')


def main(argv=None):
    """
A toolbox for research on structural variation detection.

Usage: svtools <command> [<options>...]

General Options:
  -h, --help      Show help.
  --version       Show version and exit.

Commands:
  sim             Generate an artificial genome.
  diff            Show the difference between calls of two tools
  taregion        Evaluate target regions

See 'svtools help <command>' for more information on a specific command.
    """

    args = docopt(
        main.__doc__,
        version='svtools version %s' % __version__,
        options_first=True,
        argv=argv or sys.argv[1:]
    )

    cmd = args['<command>']
    try:
        method = globals()[cmd]
        assert callable(method)
    except (KeyError, AssertionError):
        exit("%r is not a svtools command. See 'svtools help'." % cmd)

    argv = [args['<command>']] + args['<options>']
    return method(argv)

if __name__ == "__main__":
    main()

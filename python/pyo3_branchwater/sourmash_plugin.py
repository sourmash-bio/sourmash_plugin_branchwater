#! /usr/bin/env python
import sys
import argparse
from sourmash.plugins import CommandLinePlugin
from sourmash.logging import notify

from . import pyo3_branchwater

class Branchwater_Manysearch(CommandLinePlugin):
    command = 'manysearch'
    description = 'massively parallel sketch search'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_paths',
                       help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('against_paths',
                       help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('-o', '--output', required=True)
        p.add_argument('-t', '--threshold', default=0.01, type=float,
                       help="containment threshold for matches")
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help="k-mer size for which to load sketches & do search")
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help="scaled value for which to load sketches & do search")

    def main(self, args):
        num_threads = pyo3_branchwater.get_num_threads()
        notify(f"searching all sketches in '{args.query_paths}' against '{args.against_paths}' using {num_threads} threads")
        super().main(args)
        status = pyo3_branchwater.do_search(args.query_paths,
                                            args.against_paths,
                                            args.threshold,
                                            args.ksize,
                                            args.scaled,
                                            args.output)
        if status == 0:
            notify(f"...manysearch is done! results in '{args.output}'")
        return status


class Branchwater_Fastgather(CommandLinePlugin):
    command = 'fastgather'
    description = 'massively parallel sketch gather'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_sig', help="metagenome sketch")
        p.add_argument('against_paths', help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('-o', '--output-gather', required=True,
                       help="save gather output (minimum metagenome cover) to this file")
        p.add_argument('--output-prefetch',
                       help="save prefetch output (all overlaps) to this file")
        p.add_argument('-t', '--threshold-bp', default=100000, type=float)
        p.add_argument('-k', '--ksize', default=31, type=int)
        p.add_argument('-s', '--scaled', default=1000, type=int)

    def main(self, args):
        num_threads = pyo3_branchwater.get_num_threads()
        notify(f"gathering all sketches in '{args.query_sig}' against '{args.against_paths}' using {num_threads} threads")
        super().main(args)
        status = pyo3_branchwater.do_countergather(args.query_sig,
                                                   args.against_paths,
                                                   int(args.threshold_bp),
                                                   args.ksize,
                                                   args.scaled,
                                                   args.output_gather,
                                                   args.output_prefetch)
        if status == 0:
            notify(f"...fastgather is done! gather results in '{args.output_gather}'")
            if args.output_prefetch:
                notify(f"prefetch results in '{args.output_prefetch}'")
        return status


class Branchwater_Fastmultigather(CommandLinePlugin):
    command = 'fastmultigather'
    description = 'massively parallel sketch multigather'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_paths', help="a text file containing paths to .sig/.sig.gz files to query")
        p.add_argument('against_paths', help="a text file containing paths to .sig/.sig.gz files to search against")
        p.add_argument('-t', '--threshold-bp', default=100000, type=float)
        p.add_argument('-k', '--ksize', default=31, type=int)
        p.add_argument('-s', '--scaled', default=1000, type=int)

    def main(self, args):
        num_threads = pyo3_branchwater.get_num_threads()
        notify(f"gathering all sketches in '{args.query_paths}' against '{args.against_paths}' using {num_threads} threads")
        super().main(args)
        status = pyo3_branchwater.do_multigather(args.query_paths,
                                                 args.against_paths,
                                                 int(args.threshold_bp),
                                                 args.ksize,
                                                 args.scaled)
        if status == 0:
            notify(f"...fastmultigather is done!")
        return status

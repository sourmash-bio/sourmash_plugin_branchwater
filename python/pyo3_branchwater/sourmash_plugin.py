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
        notify(f"searching all sketches in '{args.query_paths}' against '{args.against_paths}'")
        super().main(args)
        pyo3_branchwater.do_search(args.query_paths,
                                    args.against_paths,
                                    args.threshold,
                                    args.ksize,
                                    args.scaled,
                                    args.output)
        notify(f"...done! results in '{args.output}'")


class Branchwater_Manygather(CommandLinePlugin):
    command = 'manygather'
    description = 'massively parallel sketch gather'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_paths')
        p.add_argument('against_paths')
        p.add_argument('-o', '--output-gather', required=True)
        p.add_argument('--output-prefetch', required=True)
        p.add_argument('-t', '--threshold-bp', default=100000, type=float)
        p.add_argument('-k', '--ksize', default=31, type=int)
        p.add_argument('-s', '--scaled', default=1000, type=int)

    def main(self, args):
        notify(f"gathering all sketches in '{args.query_paths}' against '{args.against_paths}'")
        super().main(args)
        pyo3_branchwater.do_countergather(args.query_paths,
                                           args.against_paths,
                                           int(args.threshold_bp),
                                           args.ksize,
                                           args.scaled,
                                           args.output_gather,
                                           args.output_prefetch)
        notify(f"...done! gather results in '{args.output_gather}'")
        notify(f"prefetch results in '{args.output_prefetch}'")

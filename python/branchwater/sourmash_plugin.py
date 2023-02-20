#! /usr/bin/env python
import sys
import argparse
from sourmash.plugins import CommandLinePlugin
from sourmash.logging import notify

from . import branchwater

class Branchwater_Manysearch(CommandLinePlugin):
    command = 'manysearch'
    description = 'massively parallel sketch search'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_paths')
        p.add_argument('against_paths')
        p.add_argument('-o', '--output', required=True)
        p.add_argument('-t', '--threshold', default=0.01, type=float)
        p.add_argument('-k', '--ksize', default=31, type=int)
        p.add_argument('-s', '--scaled', default=1000, type=int)

    def main(self, args):
        notify(f"searching all sketches in '{args.query_paths}' against '{args.against_paths}'")
        super().main(args)
        branchwater.do_search(args.query_paths,
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
        p.add_argument('-t', '--threshold', default=0.01, type=float)
        p.add_argument('-k', '--ksize', default=31, type=int)
        p.add_argument('-s', '--scaled', default=1000, type=int)

    def main(self, args):
        print(args)
        super().main(args)
        branchwater.do_countergather(args.query_paths,
                                     args.against_paths,
                                     args.threshold,
                                     args.ksize,
                                     args.scaled,
                                     args.output)

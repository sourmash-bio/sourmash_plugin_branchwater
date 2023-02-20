#! /usr/bin/env python
import sys
import argparse
from sourmash.plugins import CommandLinePlugin

from . import pymagsearch

class Branchwater_Extension(CommandLinePlugin):
    command = 'branchwater'
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
        print(args)
        super().main(args)
        pymagsearch.do_search(args.query_paths,
                              args.against_paths,
                              args.threshold,
                              args.ksize,
                              args.scaled,
                              args.output)

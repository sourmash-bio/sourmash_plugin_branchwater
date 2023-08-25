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
        p.add_argument('-o', '--output', required=True,
                       help='CSV output file for matches')
        p.add_argument('-t', '--threshold', default=0.01, type=float,
                       help='containment threshold for reporting matches')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to select sketches')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons')

    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold: {args.threshold}")
        num_threads = pyo3_branchwater.get_num_threads()
        notify(f"searching all sketches in '{args.query_paths}' against '{args.against_paths}' using {num_threads} threads")
        super().main(args)
        status = pyo3_branchwater.do_manysearch(args.query_paths,
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
        p.add_argument('-t', '--threshold-bp', default=50000, type=float,
                       help='threshold in estimated base pairs, for reporting matches (default: 50kb)')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to do comparisons (default: 31)')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons (default: 1000)')

    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold bp: {args.threshold_bp}")
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
        p.add_argument('-t', '--threshold-bp', default=50000, type=float,
                       help='threshold in estimated base pairs, for reporting matches (default: 50kb)')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to do comparisons (default: 31)')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons (default: 1000)')

    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold bp: {args.threshold_bp}")
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


class Branchwater_Index(CommandLinePlugin):
    command = 'index'
    description = 'Build Branchwater RevIndex'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('siglist',
                       help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('-o', '--output', required=True,
                       help='output file for the index')
        p.add_argument('-t', '--threshold', default=0.01, type=float,
                       help='containment threshold for reporting matches')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to select sketches')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons')
        p.add_argument('--save-paths', action='store_true',
                       help='save paths to signatures into index. Default: save full sig into index')

    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold: {args.threshold}")
        num_threads = pyo3_branchwater.get_num_threads()
        notify(f"indexing all sketches in '{args.siglist}'")
        super().main(args)
        status = pyo3_branchwater.do_index(args.siglist,
                                                args.ksize,
                                                args.scaled,
                                                args.threshold,
                                                args.output,
                                                args.save_paths,
                                                False) # colors - currently must be false?
        if status == 0:
            notify(f"...index is done! results in '{args.output}'")
        return status

class Branchwater_Check(CommandLinePlugin):
    command = 'check'
    description = 'Check Branchwater RevIndex'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('index',
                       help='index file')
        p.add_argument('--quick', action='store_true')

    def main(self, args):
        notify(f"checking index '{args.index}'")
        super().main(args)
        status = pyo3_branchwater.do_check(args.index, args.quick)
        if status == 0:
            notify(f"...index is consistent!")
        return status

# class Branchwater_Update(CommandLinePlugin):
#     command = 'update'
#     description = 'Update Branchwater RevIndex'

#     def __init__(self, p):
#         super().__init__(p)
#         p.add_argument('siglist',
#                        help="a text file containing paths to .sig/.sig.gz files")
#         p.add_argument('-o', '--output', required=True,
#                        help='output file for the index')
#         p.add_argument('-t', '--threshold', default=0.01, type=float,
#                        help='containment threshold for reporting matches')
#         p.add_argument('-k', '--ksize', default=31, type=int,
#                        help='k-mer size at which to select sketches')
#         p.add_argument('-s', '--scaled', default=1000, type=int,
#                        help='scaled factor at which to do comparisons')
#         p.add_argument('--save-paths', action='store_true',
#                         help='save paths to signatures into index. Default: save full sig into index')

#     def main(self, args):
#         notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold: {args.threshold}")
#         num_threads = pyo3_branchwater.get_num_threads()
#         notify(f"updating index with all sketches in '{args.siglist}'")
#         super().main(args)
#         status = pyo3_branchwater.do_update(args.siglist,
#                                                 args.ksize,
#                                                 args.scaled,
#                                                 args.threshold,
#                                                 args.output,
#                                                 args.save_paths,
#                                                 False) # colors - currently must be false
#         if status == 0:
#             notify(f"...index is done! results in '{args.output}'")
#         return status

class Branchwater_Search(CommandLinePlugin):
    command = 'search'
    description = 'Search against Branchwater RevIndex'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_paths',
                       help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('index',
                       help="an index prepared with pyo3_branchwater's sourmash scripts index")
        p.add_argument('-o', '--output',
                       help='CSV output file for matches')
        p.add_argument('-t', '--threshold-bp', default=50000, type=float,
                       help='threshold in estimated base pairs, for reporting matches (default: 50kb)')
        p.add_argument('-c', '--containment-threshold', default=0.01, type=float,
                       help='containment threshold for reporting matches')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to select sketches')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons')

    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold_bp: {args.threshold_bp}")
        num_threads = pyo3_branchwater.get_num_threads()
        notify(f"searching all sketches in '{args.query_paths}' against '{args.index}' using {num_threads} threads")
        super().main(args)
        # status = pyo3_branchwater.do_search(args.query_paths,
        status = pyo3_branchwater.do_mastiffmanysearch(args.query_paths,
                                              args.index,
                                              args.ksize,
                                              args.scaled,
                                              int(args.threshold_bp),
                                              args.containment_threshold,
                                              args.output)
        if status == 0:
            if args.output:
                notify(f"...search is done! results in '{args.output}'")
            else:
                notify(f"...search is done!")
        return status

class Branchwater_Gather(CommandLinePlugin):
    command = 'gather'
    description = 'Gather against Branchwater RevIndex'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('query_paths',
                       help="a text file containing paths to .sig/.sig.gz files")
        p.add_argument('index',
                       help="an index prepared with pyo3_branchwater's sourmash scripts index")
        p.add_argument('-o', '--output',
                       help='CSV output file for matches')
        p.add_argument('-t', '--threshold-bp', default=50000, type=float,
                       help='threshold in estimated base pairs, for reporting matches (default: 50kb)')
        p.add_argument('-k', '--ksize', default=31, type=int,
                       help='k-mer size at which to select sketches')
        p.add_argument('-s', '--scaled', default=1000, type=int,
                       help='scaled factor at which to do comparisons')

    def main(self, args):
        notify(f"ksize: {args.ksize} / scaled: {args.scaled} / threshold_bp: {args.threshold_bp}")
        num_threads = pyo3_branchwater.get_num_threads()
        notify(f"gathering all sketches in '{args.query_paths}' against '{args.index}' using {num_threads} threads")
        super().main(args)
        status = pyo3_branchwater.do_gather(args.query_paths,
                                              args.index,
                                              args.ksize,
                                              args.scaled,
                                              int(args.threshold_bp),
                                              args.output)
        if status == 0:
            if args.output:
                notify(f"...gather is done! results in '{args.output}'")
            else:
                notify(f"...gather is done!")
        return status

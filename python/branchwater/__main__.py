#! /usr/bin/env python
import sys
import argparse
from . import pymagsearch


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query_paths')
    p.add_argument('against_paths')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('-t', '--threshold', default=0.01, type=float)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-s', '--scaled', default=1000, type=int)
    args = p.parse_args()

    pymagsearch.do_search(args.query_paths,
                          args.against_paths,
                          args.threshold,
                          args.ksize,
                          args.scaled,
                          args.output)


if __name__ == '__main__':
    sys.exit(main())

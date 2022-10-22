#! /usr/bin/env python
import sys
import argparse
import pymagsearch


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query')
    p.add_argument('against_paths')
    p.add_argument('-t', '--threshold-bp', default=50000, type=int)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-s', '--scaled', default=1000, type=int)
    args = p.parse_args()

    pymagsearch.do_countergather2(args.query,
                                  args.against_paths,
                                  args.threshold_bp,
                                  args.ksize,
                                  args.scaled)


if __name__ == '__main__':
    sys.exit(main())

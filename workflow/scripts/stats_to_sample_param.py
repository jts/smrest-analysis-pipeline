#! /usr/bin/env python

import argparse
import sys

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--genome-size', required=False, type=int, default=3000000000)
    parser.add_argument('--target-coverage', required=True, type=float)
    parser.add_argument('--stats-file', required=True, type=str)
    args = parser.parse_args()

    total_length = None
    with open(args.stats_file) as f:
        for line in f:
            fields = line.rstrip().split('\t')
            if fields[0] == "SN" and fields[1] == "total length:":
                total_length = int(fields[2])

    total_coverage = float(total_length) / args.genome_size
    proportion = args.target_coverage / total_coverage
    assert(proportion < 1.0)

    print(f"{proportion:.3f}")

if __name__ == "__main__":
    main()


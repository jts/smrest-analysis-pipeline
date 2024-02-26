import sys
import csv
from collections import namedtuple
from common import *

def calculate_span_for_file(filename):
    span = 0
    with open(filename) as f:
        for line in f:
            if line[0] != '#':
                fields = line.rstrip().split()
                span += int(fields[2]) - int(fields[1])
    return span

if __name__ == "__main__":
    print("tumour_sample\ttumour_depth\tnormal_sample\tnormal_depth\ttotal_depth\tpurity\ttechnology\tcaller\tcallset_region\tspan")
    for fn in sys.argv[1:]:

        # parse the filename
        short = fn.replace("analysis/beds/", "")
        metadata = parse_metadata_from_filename(short)

        span = calculate_span_for_file(fn)
        print(f"{metadata.s1}\t"
              f"{metadata.s1_depth}\t"
              f"{metadata.s2}\t"
              f"{metadata.s2_depth}\t"
              f"{metadata.s1_depth + metadata.s2_depth}\t"
              f"{metadata.expected_purity:.3}\t"
              f"{metadata.technology}\t"
              f"{metadata.caller}\t"
              f"{metadata.short_regions}\t"
              f"{span}")

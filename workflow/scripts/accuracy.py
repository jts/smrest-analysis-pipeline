import argparse
import sys
import csv
import re
from collections import namedtuple
from common import *

Stats = namedtuple('Stats', 'tp, fn, fp, sensitivity, precision, f1')

def calculate_stats_for_file(filename, min_qual, min_vaf):
    tp = 0
    fp = 0
    fn = 0

    with open(filename) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for record in reader:
            in_truth = int(record['in_truth'])
            in_called = float(record['obs_somatic_qual']) >= min_qual and record['filter_status'] == "PASS"
            
            if float(record['min_vaf']) < min_vaf:
                continue

            if int(record['in_high_confidence_region']) != 1:
                continue

            if in_truth and in_called:
                tp += 1
            if in_truth and not in_called:
                fn += 1
            if not in_truth and in_called:
                fp += 1

        sensitivity = 0.0
        precision = 0.0
        f1 = 0.0

        if tp + fn > 0:
            sensitivity = float(tp) / float(tp + fn)

        if tp + fp > 0:
            precision = float(tp) / float(tp + fp)
        
        if tp + fp + fn > 0:
            f1 = 2 * float(tp) / (2 * tp + fp + fn)
        return Stats(tp, fn, fp, sensitivity, precision, f1)
    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--min-vaf', required=True, type=float)
    parser.add_argument('--smrest-min-qual', required=False, default=10)
    args, files = parser.parse_known_args()

    print("tumour_sample\ttumour_depth\tnormal_sample\tnormal_depth\ttotal_depth\tpurity\ttechnology\tcaller\tcallset_region\ttp\tfn\tfp\tsensitivity\tprecision\tf1")
    for fn in files:

        # parse the filename
        short = fn.replace("analysis/annotated/", "")
        metadata = parse_metadata_from_filename(short)
        min_qual = 0.0
        if metadata.caller == "smrest":
            min_qual = float(args.smrest_min_qual)
        s = calculate_stats_for_file(fn, min_qual, args.min_vaf)
        print(f"{metadata.s1}\t"
              f"{metadata.s1_depth}\t"
              f"{metadata.s2}\t"
              f"{metadata.s2_depth}\t"
              f"{metadata.s1_depth + metadata.s2_depth}\t"
              f"{metadata.expected_purity:.3}\t"
              f"{metadata.technology}\t"
              f"{metadata.caller}\t"
              f"{metadata.short_regions}\t"
              f"{s.tp}\t"
              f"{s.fn}\t"
              f"{s.fp}\t"
              f"{s.sensitivity:.3}\t"
              f"{s.precision:.3}"
              f"\t{s.f1:.3}")

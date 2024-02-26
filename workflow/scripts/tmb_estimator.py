import sys
import csv
import pysam
import glob
import os
import argparse
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

def count_pass_calls(filename):
    vcf = pysam.VariantFile(filename)
    count = 0
    sum_somatic_post = 0.0
    for record in vcf:
        f = record.filter.items()
        if len(f) == 0 or f[0][0] == "PASS":
            count += 1
            p_true = 1.0 - pow(10.0, -(record.qual / 10.0))
            sum_somatic_post += p_true
    return (count, sum_somatic_post)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--qual', required=True, type=str)
    args = parser.parse_args()

    print("tumour_sample\ttumour_depth\tnormal_sample\tnormal_depth\ttotal_depth\tpurity\ttechnology\tcaller\tcallset_region\tcallset_region_mb\thard_calls_per_mb\tsoft_calls_per_mb")
    for path in glob.glob("smrest_calls/*"):
        sn = os.path.basename(path)
        pass_vcf_file = f"smrest_calls/{sn}/{sn}.whatshap_phased.final_q{args.qual}_pass_calls.vcf"
        final_bed_file = f"smrest_calls/{sn}/{sn}.whatshap_phased.final_call_regions.bed"

        if not os.path.exists(pass_vcf_file) or not os.path.exists(final_bed_file):
            continue

        # parse the filename
        metadata = parse_metadata_from_sample_name(sn)
        (hard_calls, soft_calls) = count_pass_calls(pass_vcf_file)
        span = calculate_span_for_file(final_bed_file)

        span_mb = float(span) / 1000000
        hard_calls_per_mb = "NA"
        sort_calls_per_mb = "NA"
        if span_mb > 0:
            hard_calls_per_mb = float(hard_calls) / span_mb
            soft_calls_per_mb = float(soft_calls) / span_mb
        print(f"{metadata.s1}\t"
              f"{metadata.s1_depth}\t"
              f"{metadata.s2}\t"
              f"{metadata.s2_depth}\t"
              f"{metadata.s1_depth + metadata.s2_depth}\t"
              f"{metadata.expected_purity:.3}\t"
              f"{metadata.technology}\t"
              f"{metadata.caller}\t"
              f"{metadata.short_regions}\t"
              f"{span_mb:.5}\t"
              f"{hard_calls_per_mb:.3}\t"
              f"{soft_calls_per_mb:.3}")


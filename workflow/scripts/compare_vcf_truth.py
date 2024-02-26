#! /usr/bin/env python

from intervaltree import Interval, IntervalTree
import argparse
import random
import pysam
import csv
import sys
import os

class MutationDesc:
    def __init__(self, chromosome, position, ref, alt):
        self.chromosome = chromosome
        self.position = position
        self.ref = ref
        self.alt = alt

    def key(self):
        return (self.chromosome, self.position, self.ref, self.alt)

    def __hash__(self):
        return hash( self.key() )

    def __eq__(self, other):
        return self.key() == other.key()
    
    def __lt__(self, other):
        return self.key() < other.key()

def load_bed_to_intervaltree(filename):
    trees = dict()
    
    with open(filename) as f:
        for (idx, line) in enumerate(f):
            if line[0] == '#':
                continue
            
            fields = line.rstrip().split()
            chromosome = fields[0]
            start = fields[1]
            end = fields[2]
            if chromosome not in trees:
                trees[chromosome] = IntervalTree()
            trees[chromosome].addi(int(start), int(end), idx)
    return trees

def is_in_region(regions, chromosome, position):
    if regions is not None:
        if chromosome in regions:
            return int(len(regions[chromosome].at(position)) > 0)
    return 0

def read_snvs_from_vcf(fn):
    out = dict()
    reader = pysam.VariantFile(fn)
    for v in reader:
        if len(v.ref) > 1 or len(v.alts[0]) > 1 or len(v.alts) > 1:
            continue

        m = MutationDesc(v.chrom, v.pos, v.ref, v.alts[0])
        out[m] = v
    return out

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--truth-vcf', required=True, type=str)
    parser.add_argument('--called-vcf', required=True, type=str)
    parser.add_argument('--mask-vcf', required=False, type=str)
    parser.add_argument('--truth-sample-name', required=False, type=str)
    parser.add_argument('--tumour-sample-name', required=False, type=str)
    parser.add_argument('--analysis-bed', type=str)
    args = parser.parse_args()
    
    hc_regions = None
    if args.analysis_bed:
        hc_regions = load_bed_to_intervaltree(args.analysis_bed)

    truth_mutations = dict()
    called_mutations = dict()
    mask_mutations = dict()

    truth_reader = pysam.VariantFile(args.truth_vcf)
    for v in truth_reader:
        if "PASS" in v.filter and len(v.ref) == 1 and len(v.alts[0]) == 1:
            m = MutationDesc(v.chrom, v.pos, v.ref, v.alts[0])
            truth_mutations[m] = v
    
    called_mutations = read_snvs_from_vcf(args.called_vcf)
    if args.mask_vcf is not None:
        mask_mutations = read_snvs_from_vcf(args.mask_vcf)

    print("chromosome\tposition\tref\talt\tin_high_confidence_region\tin_truth\tin_called\tfilter_status\ttruth_vaf\ttruth_alt_count\ttruth_ref_count\tobs_vaf\tobs_hvaf\tobs_strand_bias\tobs_somatic_qual\tmin_vaf")
    all_mut = set().union(truth_mutations, called_mutations)
    for m in sorted(all_mut):
        in_truth = int(m in truth_mutations)

        truth_vaf = "NA"
        truth_alt_count = "NA"
        truth_ref_count = "NA"
        if m in truth_mutations:
            tm = truth_mutations[m]

            # different mutation calls annotated vaf in different ways
            if len(tm.ref) == 1 and len(tm.alts[0]) == 1 and args.truth_sample_name is not None and args.truth_sample_name in tm.samples:
                # attempt to parse from genotype data
                gt_data = tm.samples[args.truth_sample_name]
                if "AD" in gt_data:
                    ad = tm.samples[args.truth_sample_name]['AD']
                    truth_ref_count = ad[0]
                    truth_alt_count = ad[1]
                elif "AU" in gt_data:
                    # strelka style        
                    truth_ref_count = gt_data[tm.ref + "U"][0]
                    truth_alt_count = gt_data[tm.alts[0] + "U"][0]
            else:
                if "bwaDP" in tm.info:
                    (truth_alt_count, total_depth) = tm.info['bwaDP']
                    # not strictly correct but close enough for Illumina SNVs
                    truth_ref_count = total_depth - truth_alt_count

            if truth_ref_count != "NA" and truth_alt_count != "NA":
                truth_depth = truth_ref_count + truth_alt_count
                if truth_depth > 0:
                    truth_vaf = f"{truth_alt_count / truth_depth:.3}"
            
        obs_filter = 'NA'
        obs_mutated_hvaf = 'NA'
        
        obs_alt_count = 'NA'
        obs_ref_count = 'NA'
        obs_vaf = 'NA'
        obs_strand_bias = 'NA'
        in_called = 0
        min_vaf = 0

        # we don't want NAs in the column we use to make PR curves
        # so we set qual to 0 for all variants that are not called, or filtered out
        obs_somatic_qual = 0

        if m in called_mutations:
            cm = called_mutations[m]
            if 'SomaticHaplotypeIndex' in cm.info:
                #gather smrest stats
                shi = int(cm.info['SomaticHaplotypeIndex'])

                obs_mutated_hvaf = float(cm.info['HaplotypeVAF'][shi])
                
                obs_depth = int(cm.info['HaplotypeDepth'][0]) + int(cm.info['HaplotypeDepth'][1])
                obs_alt_count = int(cm.info['HaplotypeAltCount'][0]) + int(cm.info['HaplotypeAltCount'][1])
                obs_ref_count = obs_depth - obs_alt_count
                obs_vaf = float(obs_alt_count) / obs_depth

                obs_strand_bias = float(cm.info['StrandBias'])
                obs_somatic_qual = f"{cm.qual:.3}"
            elif 'TLOD' in cm.info:
                # gather mutect2 stats
                obs_somatic_qual = int(cm.info['TLOD'][0])
                obs_vaf = float(cm.samples[args.tumour_sample_name]['AF'][0])
            else:
                # gather clairS stats
                obs_somatic_qual = f"{cm.qual:.3}"
                obs_vaf = float(cm.samples['SAMPLE']['AF'])

            filters = [str(x) for x in cm.filter]

            # some callers have no filters for pass, others list PASS. make consistent here
            if len(filters) == 1 and filters[0] == "PASS":
                filters =  []

            if m in mask_mutations:
                filters.append("PossibleNormalCellMutation")

            obs_filter = ",".join(filters)
            if len(obs_filter) == 0:
                obs_filter = "PASS"

            # apply hard-filtering to everything except LowQual/"weak_evidence" for mutect2
            if obs_filter != "PASS" and obs_filter != "LowQual" and obs_filter != "weak_evidence":
                obs_somatic_qual = 0
                in_called = 0
            else:
                in_called = 1

        if obs_vaf == "NA" and truth_vaf == "NA":
            continue
        elif obs_vaf == "NA":
            min_vaf = truth_vaf
        elif truth_vaf == "NA":
            min_vaf = obs_vaf
        else:
            min_vaf = obs_vaf

        in_hc_region = '-'
        if hc_regions is not None:
            in_hc_region = str(is_in_region(hc_regions, m.chromosome, m.position))

        print(f"{m.chromosome}\t{m.position}\t{m.ref}\t{m.alt}\t{in_hc_region}\t{in_truth}\t{in_called}\t{obs_filter}\t{truth_vaf:.4}\t{truth_alt_count}\t{truth_ref_count}\t{obs_vaf:.4}\t{obs_mutated_hvaf:.3}\t{obs_strand_bias:.4}\t{obs_somatic_qual}\t{min_vaf:.4}")
    
if __name__ == "__main__":
    main()


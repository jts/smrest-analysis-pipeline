import re
from collections import namedtuple

ExperimentMetadata = namedtuple('ExperimentMetadata', 's1, s1_depth, s2, s2_depth, replicate, technology, caller, expected_purity, short_regions')

def parse_metadata_from_sample_name(sn):
    fields = sn.split("_")

    s1 = fields[0]
    s1_depth = float(fields[1])

    f2_idx = None
    if len(fields) == 5:
        # mixture experiment
        s2 = fields[2]
        s2_depth = float(fields[3])
        f2_idx = 4
    else:
        s2 = "NA"
        s2_depth = 0
        f2_idx = 2
    
    (rep, tech) = fields[f2_idx].split(".")
    
    total_depth = s1_depth + s2_depth
    expected_purity = float(s1_depth) / float(total_depth)

    return ExperimentMetadata(s1, s1_depth, s2, s2_depth, rep, tech, "NA", expected_purity, "NA")

def parse_metadata_from_filename(fn):
    match = re.search(r'(\w+)_(\d+(?:\.\d+)?)_(\w+)_(\d+(?:\.\d+)?)_(\w+).(\w+).(\w+).(.*)', fn)

    s1 = match.group(1)
    s1_depth = float(match.group(2))

    s2 = match.group(3)
    s2_depth = float(match.group(4))

    rep = match.group(5)
    caller = match.group(6)
    tech = match.group(7)

    regions = match.group(8)
    short_regions = regions

    expected_purity = float(s1_depth) / (float(s1_depth) + float(s2_depth))
    total_depth = s1_depth + s2_depth

    regions = regions.replace(".bed", "")
    regions = regions.replace(".tsv", "")
    short_regions = regions.replace("annotated_", "")
    
    return ExperimentMetadata(s1, s1_depth, s2, s2_depth, rep, tech, caller, expected_purity, short_regions)

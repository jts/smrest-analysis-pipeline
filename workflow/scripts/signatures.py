from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from numpy import dot
from numpy.linalg import norm
import sys

path = sys.argv[1]
name = sys.argv[2]
matrices = matGen.SigProfilerMatrixGeneratorFunc(name, "GRCh38", path, plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)

samples = list(matrices['96'].keys())
print(samples)

for a in samples:
    for b in samples:
        sa = list(matrices['96'][a])
        sb = list(matrices['96'][b])
        # https://stackoverflow.com/questions/18424228/cosine-similarity-between-2-number-lists
        cos_sim = dot(sa, sb)/(norm(sa)*norm(sb))
        print(f"{a}\t{b}\t{cos_sim:.3}")

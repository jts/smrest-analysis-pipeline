def get_calling_windows(win_size_mb):
    MB=1000000
    window_size = win_size_mb * MB
    fai = config["reference"] + ".fai"
    chromosomes = [ "chr" + str(i) for i in range(0, 23) ]
    chromosomes.append("chrX")
    chromosomes.append("chrY")
    windows = list()
    with open(fai) as f:
        for line in f:
            fields = line.rstrip().split()
            chrom_name = fields[0]
            chrom_len = int(fields[1])

            if chrom_name not in chromosomes:
                continue

            for i in range(1, chrom_len, window_size):
                end = min(chrom_len, i + window_size)
                windows.append( f"{chrom_name}:{i}-{end}" )
    return windows

def get_reference_dictionary(wildcards):
    ref = config["reference"]
    (path, ext) = os.path.splitext(ref)
    return path + ".dict"

rule index_vcf:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    params:
        memory_per_thread="8G",
        extra_cluster_opt=""
    shell:
        "tabix {input}"

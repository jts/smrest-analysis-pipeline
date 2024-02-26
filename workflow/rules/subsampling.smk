#
# Subsample read data to certain depth
#

def get_input_bam_cram(wildcards):
    bam_path = f"full_bam/{wildcards.sample}.{wildcards.tech}.bam"
    cram_path = f"full_bam/{wildcards.sample}.{wildcards.tech}.cram"

    if wildcards.tech == "ont" and (wildcards.sample == "COLO829" or wildcards.sample == "COLO829BL" or wildcards.sample == "cliveome-d4.1.0-sup"):
        return cram_path
    else:
        return bam_path

rule make_stats_file:
    input:
        bam=ancient(get_input_bam_cram)
    output:
        stats_file="etc/{sample}.{tech}.samtools_stats.txt"
    params:
        memory_per_thread="32G",
        extra_cluster_opt=""
    shell:
        "samtools stats {input.bam} > {output.stats_file}"

rule make_tmp_subset_bam:
    input:
        bam=ancient(get_input_bam_cram),
        stats_file="etc/{sample}.{tech}.samtools_stats.txt"
    output:
        #temp("data/tmp/{sample}_{cov}_rep{rep}.{tech}.bam")
        "data/tmp/{sample}_{cov}_rep{rep}.{tech}.bam"
    params:
        sample_rate_script = srcdir("../scripts/stats_to_sample_param.py"),
        memory_per_thread="8G",
        extra_cluster_opt=""
    shell:
        "samtools view -b -s `python {params.sample_rate_script} --target {wildcards.cov} --stats {input.stats_file}` {input.bam} > {output}"

rule make_merged_bam:
    input:
        tbam="data/tmp/{tsample}_{tcov}_rep{rep}.{tech}.bam",
        nbam="data/tmp/{nsample}_{ncov}_rep{rep}.{tech}.bam"
    output:
        #temp("data/{tsample}_{tcov}_{nsample}_{ncov}_rep{rep}.{tech}.bam")
        "data/{tsample}_{tcov}_{nsample}_{ncov}_rep{rep}.{tech}.bam"
    params:
        memory_per_thread="8G",
        extra_cluster_opt=""
    shell:
        "samtools merge -o /dev/stdout {input.tbam} {input.nbam} | samtools addreplacerg -o {output} -r \"@RG\tID:{wildcards.tsample}_mixture\tSM:{wildcards.tsample}_mixture\" -"

rule make_single_sample_subset_bam:
    input:
        bam=ancient(get_input_bam_cram),
        stats_file="etc/{sample}.{tech}.samtools_stats.txt"
    output:
        temp("data/{sample}_{cov}_rep{rep}.{tech}.bam")
    params:
        sample_rate_script = srcdir("../scripts/stats_to_sample_param.py"),
        memory_per_thread="8G",
        extra_cluster_opt=""
    shell:
        "samtools view -b -s `python {params.sample_rate_script} --target {wildcards.cov} --stats {input.stats_file}` {input.bam} > {output}"

ruleorder: make_merged_bam > make_single_sample_subset_bam
ruleorder: make_tmp_subset_bam > make_single_sample_subset_bam

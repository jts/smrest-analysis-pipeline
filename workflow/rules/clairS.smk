
def get_ref_dir(wildcards):
    return os.path.dirname(os.path.realpath(config['reference']))

def get_ref_path(wildcards):
    return os.path.realpath(config['reference'])

def get_platform_arg(wildcards):
    if wildcards.tech == "ont":
        return "ont_r10_dorado_4khz"
    elif wildcards.tech == "pacbio":
        return "hifi_revio"
    assert(False)

rule run_clairS:
    input:
        tbam="data/tmp/{sample}_{tdepth}_rep{rep}.{tech}.bam",
        tbai="data/tmp/{sample}_{tdepth}_rep{rep}.{tech}.bam.bai",
        nbam="data/tmp/{sample}BL_{ndepth}_rep{rep}.{tech}.bam",
        nbai="data/tmp/{sample}BL_{ndepth}_rep{rep}.{tech}.bam.bai",
        sif="clairs_latest.sif"
    output:
        dir=directory("clairS_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.{tech}"),
        vcf="clairS_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.{tech}/output.vcf.gz"
    params:
        ref_dir=get_ref_dir,
        ref_path=get_ref_path,
        platform=get_platform_arg,
        wd=os.getcwd(),
        memory_per_thread="6G",
        extra_cluster_opt=""
    threads: 16
    shell:
        """
        singularity exec -B {params.ref_dir},{params.wd} {input.sif} /opt/bin/run_clairs --tumor_bam_fn {input.tbam} --normal_bam_fn {input.nbam} --ref_fn {params.ref_path} --threads {threads} --platform {params.platform} --output_dir {output.dir} --conda_prefix /opt/conda/envs/clairs
        """

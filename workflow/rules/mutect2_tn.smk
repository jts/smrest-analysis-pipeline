
def get_normal_name(wildcards):
    if wildcards.sample == "COLO829":
        if wildcards.direction == "tn":
            return "COLO829R"
        elif wildcards.direction == "nt": # reciprocal comparison
            return "COLO829T"
        #return "CLINVAL_0001_Pb_R_COLO_839BL_CRL-1980"
        #return "CLINVAL_0001_Sk_F_COLO_829_CRL-1974"
    elif wildcards.sample == "HCC1395":
        if wildcards.direction == "tn":
            return "WGS_IL_N_1"
        elif wildcards.direction == "nt": # reciprocal comparison
            return "WGS_IL_T_1"
    assert(false)

def get_mutect2_tumour_name(wildcards):
    if wildcards.sample == "COLO829":
            return "COLO829T"
    elif wildcards.sample == "HCC1395":
            return "WGS_IL_T_1"
    assert(false)

rule run_gatk_create_dictionary:
    input:
        "{prefix}.fna"
    output:
        "{prefix}.dict"
    params:
        memory_per_thread="16G",
        extra_cluster_opt=""
    shell:
        "{config[gatk]}/gatk CreateSequenceDictionary -R {input}"

rule run_mutect2_tn_region:
    input:
        tbam="data/tmp/{sample}_{tdepth}_rep{rep}.illumina.bam",
        tbai="data/tmp/{sample}_{tdepth}_rep{rep}.illumina.bam.bai",
        nbam="data/tmp/{sample}BL_{ndepth}_rep{rep}.illumina.bam",
        nbai="data/tmp/{sample}BL_{ndepth}_rep{rep}.illumina.bam.bai",
        gnomad_vcf=config["gnomad_af_vcf"],
        gnomad_vcf_tbi=config["gnomad_af_vcf"] + ".tbi",
        pon_vcf=config["mutect2_pon"],
        pon_vcf_tbi=config["mutect2_pon"] + ".tbi",
        pon_tbi="data/tmp/{sample}BL_{ndepth}_rep{rep}.illumina.bam.bai",
        dict=get_reference_dictionary
    output:
        vcf=temp("mutect2_{direction}_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina/per_region/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina.region.{region}.vcf.gz"),
        stats=temp("mutect2_{direction}_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina/per_region/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina.region.{region}.vcf.gz.stats")
    params:
        memory_per_thread="12G",
        extra_cluster_opt="",
        normal_sample_name=get_normal_name
    threads: 8
    shell:
        """
        {config[gatk]}/gatk Mutect2 -R {config[reference]} \
            -I {input.tbam} \
            -I {input.nbam} \
            -normal {params.normal_sample_name} \
            -O {output.vcf} \
            -L {wildcards.region} \
            --native-pair-hmm-threads {threads} \
            --pon {config[mutect2_pon]} \
            --germline-resource {config[gnomad_af_vcf]} \
            --max-mnp-distance 0
        """

rule mutect2_tn_pileup_table:
    input:
        bam="data/tmp/{S}_{D}_rep{R}.illumina.bam",
        bai="data/tmp/{S}_{D}_rep{R}.illumina.bam.bai",
        vcf=config["gnomad_common"],
        tbi=config["gnomad_common"] + ".tbi"
    output:
        "mutect2_{direction}_calls/pileups/{S}_{D}_rep{R}.illumina.pileup.table"
    params:
        memory_per_thread="96G",
        extra_cluster_opt=""
    shell:
        "{config[gatk]}/gatk GetPileupSummaries -I {input.bam} -V {input.vcf} -L {input.vcf} -O {output}"

rule mutect2_tn_calculate_contamination:
    input:
        T_table="mutect2_{direction}_calls/pileups/{sample}_{tdepth}_rep{rep}.illumina.pileup.table",
        N_table="mutect2_{direction}_calls/pileups/{sample}BL_{ndepth}_rep{rep}.illumina.pileup.table"
    output:
        contam="mutect2_{direction}_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina.contamination.table",
        segmentation="mutect2_{direction}_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina.segmentation.tsv",
    params:
        memory_per_thread="24G",
        extra_cluster_opt=""
    shell:
        "{config[gatk]}/gatk CalculateContamination -I {input.T_table} -matched {input.N_table} -O {output.contam} --tumor-segmentation {output.segmentation}"

rule filter_mutect2_tn_calls:
    input:
        vcf="mutect2_{direction}_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina.calls.raw.vcf.gz",
        tbi="mutect2_{direction}_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina.calls.raw.vcf.gz.tbi",
        stats="mutect2_{direction}_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina.calls.raw.vcf.gz.stats",
        contam="mutect2_{direction}_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina.contamination.table",
        segmentation="mutect2_{direction}_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina.segmentation.tsv",
    output:
        "mutect2_{direction}_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.illumina.calls.filtered.vcf.gz"
    params:
        memory_per_thread="8G",
        extra_cluster_opt=""
    shell:
        """
        {config[gatk]}/gatk FilterMutectCalls -R {config[reference]} \
                                               -V {input.vcf} \
                                               -O {output} \
                                               --contamination-table {input.contam} \
                                               --tumor-segmentation {input.segmentation}
        """

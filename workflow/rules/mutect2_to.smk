rule run_mutect2_region:
    input:
        bam="data/{sample}.bam",
        bai="data/{sample}.bam.bai",
        dict=get_reference_dictionary,
        gnomad_vcf=config["gnomad_af_vcf"],
        gnomad_vcf_tbi=config["gnomad_af_vcf"] + ".tbi",
        pon_vcf=config["mutect2_pon"],
        pon_vcf_tbi=config["mutect2_pon"] + ".tbi"
    output:
        vcf=temp("mutect2_to_calls/{sample}/per_region/{sample}.region.{region}.vcf.gz"),
        stats=temp("mutect2_to_calls/{sample}/per_region/{sample}.region.{region}.vcf.gz.stats")
    params:
        memory_per_thread="12G",
        extra_cluster_opt=""
    threads: 8
    shell:
        """
        {config[gatk]}/gatk Mutect2 -R {config[reference]} \
            -I {input.bam} \
            -O {output.vcf} \
            -L {wildcards.region} \
            --native-pair-hmm-threads {threads} \
            --pon {config[mutect2_pon]} \
            --germline-resource {config[gnomad_af_vcf]} \
            --max-mnp-distance 0
        """

rule merge_mutect2_region_calls:
    input:
        vcfs=expand("mutect2_{{type}}_calls/{{sample}}/per_region/{{sample}}.region.{r}.vcf.gz", r=get_calling_windows(mutect2_window_size))
    output:
        "mutect2_{type}_calls/{sample}/{sample}.calls.raw.vcf.gz"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "bcftools concat {input.vcfs} | bcftools sort | bgzip > {output}"

rule merge_mutect2_stats:
    input:
        stats=expand("mutect2_{{type}}_calls/{{sample}}/per_region/{{sample}}.region.{r}.vcf.gz.stats", r=get_calling_windows(mutect2_window_size))
    output:
        "mutect2_{type}_calls/{sample}/{sample}.calls.raw.vcf.gz.stats"
    params:
        memory_per_thread="8G",
        extra_cluster_opt=""
    run:
        statsarg = " ".join([f"--stats {i}" for i in input.stats ])
        shell("{config[gatk]}/gatk MergeMutectStats %s -O {output}" % (statsarg))

rule mutect_pileup_table:
    input:
        bam="data/{sample}.bam",
        bai="data/{sample}.bam.bai",
        vcf=config["gnomad_common"],
        tbi=config["gnomad_common"] + ".tbi"
    output:
        "mutect2_to_calls/{sample}/{sample}.pileup.table"
    params:
        memory_per_thread="96G",
        extra_cluster_opt=""
    shell:
        "{config[gatk]}/gatk GetPileupSummaries -I {input.bam} -V {input.vcf} -L {input.vcf} -O {output}"

rule mutect2_calculate_contamination:
    input:
        "mutect2_to_calls/{sample}/{sample}.pileup.table"
    output:
        contam="mutect2_to_calls/{sample}/{sample}.contamination.table",
        segmentation="mutect2_to_calls/{sample}/{sample}.segmentation.tsv",
    params:
        memory_per_thread="24G",
        extra_cluster_opt=""
    shell:
        "{config[gatk]}/gatk CalculateContamination -I {input} -O {output.contam} --tumor-segmentation {output.segmentation}"

rule filter_mutect2_calls:
    input:
        vcf="mutect2_to_calls/{sample}/{sample}.calls.raw.vcf.gz",
        tbi="mutect2_to_calls/{sample}/{sample}.calls.raw.vcf.gz.tbi",
        stats="mutect2_to_calls/{sample}/{sample}.calls.raw.vcf.gz.stats",
        contam="mutect2_to_calls/{sample}/{sample}.contamination.table",
        segmentation="mutect2_to_calls/{sample}/{sample}.segmentation.tsv",
    output:
        "mutect2_to_calls/{sample}/{sample}.calls.filtered.vcf.gz"
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

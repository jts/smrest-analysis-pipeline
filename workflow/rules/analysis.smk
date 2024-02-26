#
# Calculate post-run mutation annotations
#
rule make_called_bed:
    input:
        called_bed="smrest_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.ont/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.{tech}.whatshap_phased.called_regions.bed"
    output:
        bed="analysis/beds/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.HC-{tech}-called.bed",
    params:
        memory_per_thread="16G",
        extra_cluster_opt="",
        bed = lambda w: config['analysis_region_bed'][w.sample]
    shell:
        """bedtools intersect -a {params.bed} \
                              -b {input.called_bed} > {output.bed}"""

rule make_dual_called_bed:
    input:
        ont_called_bed="smrest_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.ont/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.ont.whatshap_phased.called_regions.bed",
        pacbio_called_bed="smrest_calls/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.pacbio/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.pacbio.whatshap_phased.called_regions.bed"
    output:
        bed="analysis/beds/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.HC-both-called.bed",
    params:
        memory_per_thread="16G",
        extra_cluster_opt="",
        bed = lambda w: config['analysis_region_bed'][w.sample]
    shell:
        """bedtools intersect -a {params.bed} \
                              -b {input.ont_called_bed} |\
           bedtools intersect -a {input.pacbio_called_bed} -b stdin > {output.bed}"""

rule make_best_practice_bed:
    output:
        "analysis/beds/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.HC.bed",
    params:
        memory_per_thread="4G",
        extra_cluster_opt="",
        bed = lambda w: config['analysis_region_bed'][w.sample]
    shell:
        "ln -s -r {params.bed} {output}"

rule copy_callset_smrest:
    input:
        vcf="smrest_calls/{sample_prefix}.{tech}/{sample_prefix}.{tech}.whatshap_phased.calls.vcf"
    output:
        vcf="analysis/callsets/{sample_prefix}/{sample_prefix}.smrest.{tech}.vcf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "cp {input} {output}"

rule copy_callset_clairS:
    input:
        vcf="clairS_calls/{sample_prefix}.{tech}/output.vcf.gz"
    output:
        vcf="analysis/callsets/{sample_prefix}/{sample_prefix}.clairS.{tech}.vcf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "zcat {input} > {output}"

rule copy_callset_mutect2:
    input:
        vcf="mutect2_{type}_calls/{sample_prefix}.illumina/{sample_prefix}.illumina.calls.filtered.vcf.gz"
    output:
        vcf="analysis/callsets/{sample_prefix}/{sample_prefix}.mutect2_{type}.illumina.vcf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "zcat {input} > {output}"

#
# derive a set of potential mutations in the normal cell line
# these are used to avoid calling TO false positives
#
rule make_nt_germline_bed:
    input:
        "mutect2_nt_calls/{sample}_40_{sample}BL_40_rep1.illumina/{sample}_40_{sample}BL_40_rep1.illumina.calls.filtered.vcf.gz"
    output:
        "etc/{sample}.normal_cellline_population_calls.bed"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        """
        bcftools filter -i \"INFO/POPAF < 5 && FILTER='PASS' && TYPE='snp'\" {input} |\
        awk '{{ if($1 !~ "#") {{ print $1 "\t" ($2 - 5000) "\t" ($2 + 5001) }} }}' > {output}
        """

rule make_nt_mutation_vcf:
    input:
        vcf="mutect2_nt_calls/{sample}_40_{sample}BL_40_rep1.illumina/{sample}_40_{sample}BL_40_rep1.illumina.calls.filtered.vcf.gz",
        bed="etc/{sample}.normal_cellline_population_calls.bed"
    output:
        "etc/{sample}.normal_cellline_mutation_calls.vcf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        """
        bcftools filter --soft-filter NEARHET --mask-file {input.bed} {input.vcf} | bcftools filter -i 'INFO/POPAF >= 5 && FILTER="PASS" && TYPE="snp"' > {output}
        """

# some truth VCFs use the genotype field and hence need to be provided with a sample name
def get_truth_sample_name_arg(wildcards):
    if wildcards.sample in config["truth_sample_name"]:
        return "--truth-sample-name " + config["truth_sample_name"][wildcards.sample]
    else:
        return ""

# needed to extract mutect2 AF
def get_tumour_sample_name_arg(wildcards):
    if wildcards.caller == "mutect2_tn":
        return "--tumour-sample-name " + get_mutect2_tumour_name(wildcards)
    elif wildcards.caller == "mutect2_to":
        return "--tumour-sample-name " + wildcards.sample + "_mixture"
    else:
        return ""

# we intentially always use the ONT bed file to define the callset regions
# this is so we are comparing the exact same set of truth mutations for each
# technology/program
rule annotate_calls:
    input:
        vcf="analysis/callsets/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.{caller}.{tech}.vcf",
        bed="analysis/beds/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.{bedtype}.bed",
        mask_vcf="etc/{sample}.normal_cellline_mutation_calls.vcf"
    output:
        tsv="analysis/annotated/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.{caller}.{tech}.annotated-{bedtype}.tsv"
    params:
        memory_per_thread="24G",
        extra_cluster_opt="",
        truth_vcf = lambda w: config['truth_vcf'][w.sample],
        truth_sample_name_arg = get_truth_sample_name_arg,
        tumour_sample_name_arg = get_tumour_sample_name_arg,
        compare_vcf=srcdir("../scripts/compare_vcf_truth.py")
    shell:
        """python {params.compare_vcf} --truth-vcf {params.truth_vcf} {params.truth_sample_name_arg} {params.tumour_sample_name_arg}\
                                       --called {input.vcf} \
                                       --mask-vcf {input.mask_vcf} \
                                       --analysis-bed {input.bed} > {output.tsv}"""

rule all_callset_stats_tsv:
    # ref: https://bioinformatics.stackexchange.com/questions/7184/mix-globbing-and-wildcards-when-specifying-rule-input
    #input:
    #    lambda wildcards: glob.glob('analysis/annotated/{sample}*.tsv'.format(sample=wildcards.sample))
    output:
        "analysis/summary/{sample}.all_callset_accuracy.tsv"
    params:
        accuracy_script = srcdir("../scripts/accuracy.py"),
        memory_per_thread="24G",
        extra_cluster_opt=""
    shell:
        "python {params.accuracy_script} --min-vaf {config[min_vaf]} --smrest-min-qual {config[min_somatic_qual]} analysis/annotated/{wildcards.sample}*.tsv > {output}"

rule all_tmb_stats:
    output:
        "analysis/summary/{sample}.q{qual}.tmb_estimates.tsv"
    params:
        tmb_script = srcdir("../scripts/tmb_estimator.py"),
        memory_per_thread="24G",
        extra_cluster_opt=""
    shell:
        "python {params.tmb_script} --qual {wildcards.qual} > {output}"

#
# Mutation signatures
#
rule prepare_sites_signatures:
    input:
        tsv="analysis/annotated/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.smrest.{tech}.annotated-HC-{tech}-called.tsv"
    output:
        sites="analysis/signatures/{sample}_{tdepth}_{sample}BL_{ndepth}_rep{rep}.smrest.{tech}.sites.txt"
    shell:
        """ awk '$5 == 1 && $16 >= {config[min_vaf]}' {input.tsv} > {output}"""

rule prepare_smrest_calls_signatures:
    input:
        vcf="analysis/callsets/{sample}_40_{sample}BL_40_rep1/{sample}_40_{sample}BL_40_rep1.smrest.ont.vcf",
        sites="analysis/signatures/{sample}_40_{sample}BL_40_rep1.smrest.ont.sites.txt"
    output:
        vcf="analysis/signatures/{sample}/{sample}-smrest.vcf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        """bcftools filter -T {input.sites} -i \"QUAL >= {config[min_somatic_qual]} && FILTER='.'\" {input.vcf} > {output.vcf}"""

rule prepare_truth_calls_signatures:
    input:
        sites="analysis/signatures/{sample}_40_{sample}BL_40_rep1.smrest.ont.sites.txt"
    output:
        vcf="analysis/signatures/{sample}/{sample}-truth.vcf"
    params:
        truth_vcf = lambda w: config['truth_vcf'][w.sample],
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        """bcftools filter -T {input.sites} {params.truth_vcf} > {output.vcf}"""

rule make_signature_plot:
    input:
        call_vcf="analysis/signatures/{sample}/{sample}-smrest.vcf",
        truth_vcf="analysis/signatures/{sample}/{sample}-truth.vcf"
    output:
        directory("analysis/signatures/{sample}/output")
    params:
        script = srcdir("../scripts/signatures.py"),
        memory_per_thread="24G",
        extra_cluster_opt=""
    shell:
        "python {params.script} analysis/signatures/{wildcards.sample} {wildcards.sample}"

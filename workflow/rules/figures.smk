rule purity_plot:
    input:
        "analysis/summary/{sample}.best_practice_called_regions.all_callset_accuracy.tsv"
    output:
        "plots/publication_figure_purity_{sample}_{metric}.pdf"
    params:
        script=srcdir("../scripts/plot_colo829_purity_sweep.R"),
    shell:
        "Rscript {params.script} -o {output} -f {input} --regions 'NA' -m {wildcards.metric}"

rule plot_comparisons:
    output:
        "plots/publication_{figure}_comparison_{sample}_{depth}_{region}_{type}.pdf"
    params:
        script=srcdir("../scripts/plot_comparison_pr.R"),
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "Rscript {params.script} --clean --sample {wildcards.sample} --min-vaf {config[min_vaf]} -o {output} -d {wildcards.depth} -r {wildcards.region} -t {wildcards.type}"

rule plot_VAFs:
    output:
        "plots/publication_{figure}_{sample}_truth_vaf.pdf"
    params:
        script=srcdir("../scripts/plot_vaf.R"),
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "Rscript {params.script} --sample {wildcards.sample} --min-vaf {config[min_vaf]} -o {output}"

rule plot_purity_sweep:
    output:
        "plots/publication_{figure}_{sample}_purity_{metric}.pdf"
    params:
        script=srcdir("../scripts/plot_purity_sweep.R"),
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        """Rscript {params.script} --sample {wildcards.sample} \
                                   --metric {wildcards.metric} \
                                   -o {output} \
                                   -f analysis/summary/{wildcards.sample}.all_callset_accuracy.tsv \
                                   --regions annotate-HC-ont-called"""

rule plot_tmb_blood:
    output:
        "plots/publication_figure5_tmb_blood.pdf"
    params:
        script=srcdir("../scripts/plot_TMB.R"),
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        """Rscript {params.script} -o {output} \
                                   --blood-only \
                                   -f analysis/summary/tmb_estimates.tsv"""

rule plot_tmb_colo829:
    output:
        "plots/publication_figure5_tmb_colo829.pdf"
    params:
        script=srcdir("../scripts/plot_TMB.R"),
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        """Rscript {params.script} -o {output} \
                                   -f analysis/summary/tmb_estimates.tsv"""

# generators
rule figure2_plots:
    input:
        expand("plots/publication_figure2_comparison_COLO829_{depth}_{region}_LR-vs-SR.pdf", depth = [ 20, 30, 40], region = [ "HC", "HC-ont-called" ])

rule figure3_plots:
    input:
        expand("plots/publication_figure3_comparison_{sample}_{depth}_{region}_ONT-and-PB.pdf", sample = [ "COLO829", "HCC1395" ], depth = [ 40 ], region = [ "HC", "HC-both-called" ]),
        expand("plots/publication_figure3_{sample}_truth_vaf.pdf", sample = [ "COLO829", "HCC1395" ])

rule figure4_plots:
    input:
        expand("plots/publication_figure4_COLO829_purity_{m}.pdf", m = [ "f1", "sens_prec" ])

rule figure5_plots:
    input:
        "plots/publication_figure5_tmb_blood.pdf", "plots/publication_figure5_tmb_colo829.pdf"

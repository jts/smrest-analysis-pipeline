rule all_simulation_purity_depth:
    input:
        expand("simulations/simulation_tmb{tmb}_error{er}_clonal{clonal_prop}/sim_output_purity{purity}_depth{depth}.tsv",
            tmb = [5.0],
            clonal_prop = [1.0, 0.75],
            er = [0.01],
            purity = np.linspace(0,1,21),
            depth = range(5, 150, 5))


rule all_simulation_purity_sweep:
    input:
        expand("simulations/purity_sweep/simulation_tmb{tmb}_error{er}_clonal{clonal_prop}/sim_output_purity{purity}_depth{depth}.tsv",
            tmb = [5.0],
            clonal_prop = [1.0, 0.75],
            er = [0.01, 0.05],
            purity = np.linspace(0,1,101),
            depth = [40, 80, 160])

rule all_plots:
    input:
        expand("plots/simulation_tmb{tmb}_error{er}.pdf",
            tmb = [1.0, 5.0, 10.0],
            er = [0.01, 0.05])

rule pileup_sim:
    output:
        "simulations/{set}/simulation_tmb{tmb}_error{er}_clonal{cp}/sim_output_purity{purity}_depth{depth}.tsv"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "{config[smrest]} sim-pileup --tmb {wildcards.tmb} --proportion-clonal {wildcards.cp} --error-rate {wildcards.er} --purity {wildcards.purity} --depth {wildcards.depth} -g {config[reference]} > {output}"

rule merge_sim:
    output:
        "simulations/merged_sim/simulation_tmb{tmb}_error{er}_clonal{cp}.tsv"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "cat simulations/simulation_tmb{wildcards.tmb}_error{wildcards.er}_clonal{wildcards.cp}/*.tsv | awk 'NR == 1 || $0 !~ /sens/' > {output}"

rule plot_sim:
    input:
        "simulations/merged_sim/simulation_tmb{tmb}_error{er}_clonal{cp}.tsv"
    output:
        "plots/simulation_tmb{tmb}_error{er}_clonal{cp}.pdf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt="",
        plot_r = srcdir("plot_simulations.R")
    shell:
        "Rscript {params.plot_r} -t {wildcards.tmb} -p {wildcards.cp} -e {wildcards.er} -s {input} -o {output}"


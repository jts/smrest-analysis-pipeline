# smrest analysis pipeline

This is a snakemake workflow that runs the experiments described in [Simpson, J.T., Detecting Somatic Mutations Without Matched Normal Samples Using Long Reads, BioRxiv](https://www.biorxiv.org/content/10.1101/2024.02.26.582089v1).

## Prerequisites

The workflow will try to download most software and resources it needs to run but some common software is assumed to already be installed:

- minimap2
- bwa
- samtools
- bcftools

## Configuration

Prior to running, the pipeline needs to be configured with paths to resources and software. The paths that need to be set are indicated at the top of `config.yaml`. 

## Usage

Once configured, the workflow can be executed to generate the results in the preprint. The main rules are `init` (to download datasets and resources), `comparison_experiments` and `tmb_experiments`. The `init` rule needs to be run from a server that has an internet connection, the others do not. The complete pipeline will run >100,000 tasks so a HPC environment is needed.

Example usage:

```
snakemake -s /path/to/smrest-analysis-pipeline/Snakefile comparison_experiments
```

## License

MIT license

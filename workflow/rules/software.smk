rule get_gatk:
    output:
        "software/gatk-4.4.0.0.zip"
    shell:
        "wget -O {output} https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip"

rule extract_gatk:
    input:
        "software/gatk-4.4.0.0.zip"
    output:
        directory("software/gatk-4.4.0.0")
    shell:
        "unzip -d software {input}"

rule get_sratools:
    output:
        "software/sratoolkit.3.0.7-ubuntu64.tar.gz"
    shell:
        "wget -O {output} https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.7/sratoolkit.3.0.7-ubuntu64.tar.gz"

rule extract_sratools:
    input:
        "software/sratoolkit.3.0.7-ubuntu64.tar.gz"
    output:
        "software/sratoolkit.3.0.7-ubuntu64/bin/fastq-dump"
    shell:
        "tar -xzf {input} -C software"

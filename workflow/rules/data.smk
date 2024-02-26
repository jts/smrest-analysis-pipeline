#
# ONT
#
def get_flowcell_crams_for_sample(wildcards):
    fc = list()
    if wildcards.sample == "COLO829":
        fc = [ "PAK76302", "PAO29420", "PAO32033" ]
    elif wildcards.sample == "COLO829BL":
        fc = [ "PAK76487", "PAO33946" ]
    elif wildcards.sample == "hg002":
        fc = [ "PAO83395", "PAO89685" ]
    elif wildcards.sample == "cliveome":
        fc = [ "PAM63114", "PAM63167", "PAM63974" ]
    elif wildcards.sample == "cliveome-d4.1.0-sup":
        fc = [ "ONLA29132", "ONLA29133", "ONLA29134" ]
    else:
        assert(false)

    return expand("ont_open_data/{sample}/{flowcell}.cram", sample=wildcards.sample, flowcell = fc)

rule merged_ont_cram:
    input:
        get_flowcell_crams_for_sample
    output:
        cram="full_bam/{sample}.ont.cram",
        crai="full_bam/{sample}.ont.cram.crai"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "samtools merge --write-index -o {output.cram} {input}"

rule get_colo829_flowcell_cram:
    output:
        "ont_open_data/COLO829/{flowcell}.cram"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829/{wildcards.flowcell}/cram/{wildcards.flowcell}.cram {output}"

rule get_colo829bl_flowcell_cram:
    output:
        "ont_open_data/COLO829BL/{flowcell}.cram"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "aws s3 cp --no-sign-request s3://ont-open-data/colo829_2023.04/COLO829BL/{wildcards.flowcell}/cram/{wildcards.flowcell}.cram {output}"

rule get_giab_flowcell_cram:
    output:
        "ont_open_data/hg002/{flowcell}.cram"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "aws s3 cp --no-sign-request s3://ont-open-data/giab_2023.05/analysis/hg002/sup/{wildcards.flowcell}.pass.cram {output}"

rule get_cliveome_flowcell_bam:
    output:
        temp("ont_open_data/cliveome/{flowcell}.bam")
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "aws s3 cp --no-sign-request s3://ont-open-data/cliveome_kit14_2022.05/gdna/basecalling/{wildcards.flowcell}/bonito_calls.bam {output}"

rule get_HCC1395_bam:
    output:
        "full_bam/HCC1395.ont.bam"
    threads: 24
    shell:
        "{config[fastq_dump]} --stdout SRR25005626 | minimap2 -t 24 -a -x map-ont {config[reference]} - | samtools sort -O bam -o {output}"

rule get_HCC1395BL_bam:
    output:
        "full_bam/HCC1395BL.ont.bam"
    threads: 24
    shell:
        "{config[fastq_dump]} --stdout SRR25005625 | minimap2 -t 24 -a -x map-ont {config[reference]} - | samtools sort -O bam -o {output}"

rule ont_bam2cram:
    input:
        "ont_open_data/{sample}/{flowcell}.bam"
    output:
        "ont_open_data/{sample}/{flowcell}.cram"
    shell:
        "samtools view -o {output} {input}"

#
# Illumina
#
rule get_COLO829_ilmn_bam:
    output:
        "full_bam/tmp/COLO829.illumina.GRCh37.bam"
    shell:
        "ascp -T -l 300m -P 33001 -i {config[aspera_id]} era-fasp@fasp.sra.ebi.ac.uk:vol1/run/ERR275/ERR2752450/COLO829T_dedup.realigned.bam {output}"

rule get_COLO829BL_ilmn_bam:
    output:
        "full_bam/tmp/COLO829BL.illumina.GRCh37.bam"
    shell:
        "ascp -T -l 300m -P 33001 -i {config[aspera_id]} era-fasp@fasp.sra.ebi.ac.uk:vol1/run/ERR275/ERR2752449/COLO829R_dedup.realigned.bam {output}"

rule remap_ilmn:
    input:
        "full_bam/tmp/{sample}.illumina.GRCh37.bam"
    output:
        "full_bam/{sample}.illumina.bam"
    threads: 16
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        """
        samtools collate -Oun128 {input} | samtools fastq -OT RG,BC - \
          | bwa mem -pt{threads} -CH <(samtools view -H {input}|grep ^@RG) {config[reference]} - \
          | samtools sort -@4 -m4g -o {output} -
        """

rule get_HCC1395_ilmn_bam:
    output:
        "full_bam/HCC1395.illumina.bam"
    shell:
        "wget -O {output} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_IL_T_1.bwa.dedup.bam"

rule get_HCC1395BL_ilmn_bam:
    output:
        "full_bam/HCC1395BL.illumina.bam"
    shell:
        "wget -O {output} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_IL_N_1.bwa.dedup.bam"

#
# Pacbio
#
def get_pacbio_run_bams_for_sample(wildcards):
    fc = list()
    if wildcards.sample == "COLO829":
        fc = [ "m84039_230312_025934_s1.hifi_reads.bc2026", "m84039_230328_000836_s3.hifi_reads.bc2025" ]
    elif wildcards.sample == "COLO829BL":
        fc = [ "m84039_230327_230708_s1.hifi_reads.bc2007", "m84039_230327_233730_s2.hifi_reads.bc2007" ]
    elif wildcards.sample == "HCC1395":
        fc = [ "m84039_230414_235240_s2.hifi_reads.default", "m84039_230415_002321_s3.hifi_reads.default" ]
    elif wildcards.sample == "HCC1395BL":
        fc = [ "m84039_230412_183900_s2.hifi_reads.bc2002", "m84039_230415_005427_s4.hifi_reads.bc2002" ]
    else:
        assert(false)

    return expand("pacbio_data/{sample}/{flowcell}.bam", sample=wildcards.sample, flowcell = fc)

def get_pacbio_sample_dir(wildcards):
    if wildcards.sample == "COLO829BL":
        return "COLO829-BL"
    elif wildcards.sample == "COLO829":
        return "COLO829"
    elif wildcards.sample == "HCC1395BL":
        return "HCC1395-BL"
    elif wildcards.sample == "HCC1395":
        return "HCC1395"
    return "unknown"

def get_pacbio_base_dir(wildcards):
    return wildcards.sample.replace("BL", "")

rule get_pacbio_run:
    output:
        "pacbio_data/{sample}/{flowcell}.bam"
    params:
        base_dir=get_pacbio_base_dir,
        sample_dir=get_pacbio_sample_dir
    shell:
        "wget -O {output} https://downloads.pacbcloud.com/public/revio/2023Q2/{params.base_dir}/{params.sample_dir}/{wildcards.flowcell}.bam"

rule merged_pacbio_bam:
    input:
        get_pacbio_run_bams_for_sample
    output:
        bam="full_bam/{sample}.pacbio.unmapped.bam"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "samtools merge -o {output.bam} {input}"

rule mapped_pacbio_bam:
    input:
        bam="full_bam/{sample}.pacbio.unmapped.bam"
    output:
        bam="full_bam/{sample}.pacbio.bam"
    threads: 16
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        """
        samtools fastq -T '*' {input} \
          | minimap2 -a -t {threads} -x map-hifi -y {config[reference]} - \
          | samtools sort -@4 -m4g -o {output} -
        """

#
# Truth data
#
rule get_COLO829_truth:
    output:
        "resources/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf"
    shell:
        "wget -O {output} https://bioinformatics.nygenome.org/wp-content/uploads/CancerCellLines/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf"

rule bcftools_norm_truth:
    input:
        "resources/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf"
    output:
        "resources/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.primatives.vcf"
    shell:
        "bcftools norm -a {input} > {output}"

rule get_HCC1395_truth:
    output:
        "resources/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz"
    shell:
        "wget -O {output} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/v1.2/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz"

rule unzip_HCC1395_truth:
    input:
        "resources/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz"
    output:
        "resources/HCC1395-high-confidence_sSNV_in_HC_regions_v1.2.vcf"
    shell:
        "gunzip -c {input} > {output}"

rule get_SEQC_bed:
    output:
        "resources/High-Confidence_Regions_v1.2.bed"
    shell:
        "wget -O {output} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/v1.2/High-Confidence_Regions_v1.2.bed"

rule get_GIAB_SEQC_bed:
    input:
        SEQC="resources/High-Confidence_Regions_v1.2.bed",
        GIAB="resources/GRCh38_notinalldifficultregions.bed"
    output:
        "resources/GIAB_SEQC_regions.bed"
    shell:
        "bedtools intersect -a {input.GIAB} -b {input.SEQC} > {output}"

# Mutect2 resources
rule get_mutect2_pon:
    output:
        "resources/1000g_pon.hg38.vcf.gz"
    shell:
        "wget -O {output} https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"

# clairS image
rule get_clairs_sif:
    output:
        "clairs_latest.sif"
    shell:
        "singularity pull docker://hkubal/clairs:latest"

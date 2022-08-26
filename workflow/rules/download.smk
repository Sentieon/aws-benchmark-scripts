#!/usr/bin/env snakemake

"""
Download reference and input files
"""

#### Output ####
FASTQ = "input_fq/{sample}/{read_idx}/{pair}/data.fq.gz"
CRAM = "input_cram/{sample}/ultima.cram"
REFERENCE = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
CRAM_REF = "reference/Homo_sapiens_assembly38.fasta"
AUTOSOMES = "autosomes_bed/autosomes.fa"
TRUTHSET = "truthset/v4.2.1/HG002/truthset"
SITE_VCF = "sites_vcfs/{vcf_name}/calls.vcf.gz"
DNASCOPE_MODEL = "dnascope_model/{platform}/dnascope.model"

DNASCOPE_LR_SCRIPT = "dnascope_lr_script/DNAscopeHiFiBeta0.4.pipeline/dnascope_HiFi.sh"
DNASCOPE_LR_MODEL = "dnascope_lr_script/DNAscopeHiFiBeta0.4.pipeline/DNAscopeHiFiBeta0.4.model"

#### Rules ####
rule download_samples:
  output:
    fq = FASTQ
  log:
    stdout = FASTQ + ".stdout",
    stderr = FASTQ + ".stderr",
  params:
    url = lambda wldc: config["input"]["samples"][wldc.sample][int(wldc.read_idx)]["fq"][int(wldc.pair)]
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.fq}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    curl -L -o "{output.fq}" "{params.url}"
    """

rule download_ultima:
  output:
    cram = CRAM,
    crai = CRAM + ".crai",
  log:
    stdout = CRAM + ".stdout",
    stderr = CRAM + ".stderr",
  params:
    url = lambda wldc: config["input"]["samples"][wldc.sample]["url"],
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.cram}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    curl -L -o "{output.cram}" "{params.url}"
    curl -L -o "{output.crai}" "{params.url}.crai"
    """

rule download_truthset:
  output:
    vcf = TRUTHSET + ".vcf.gz",
    vcf_idx = TRUTHSET + ".vcf.gz.tbi",
    bed = TRUTHSET + ".bed",
    lcr_bed = TRUTHSET + "_ultima-lcr.bed",
  log:
    stdout = TRUTHSET + ".stdout",
    stderr = TRUTHSET + ".stderr",
  params:
    vcf_url = "https://ftp.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
    bed_url = "https://ftp.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
    lcr_bed_url = "https://s3.amazonaws.com/ultima-selected-1k-genomes/beds/ug_lcr.bed",
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.vcf}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    curl -L -o "{output.vcf}" "{params.vcf_url}"
    curl -L -o "{output.vcf_idx}" "{params.vcf_url}.tbi"
    curl -L -o "{output.bed}" "{params.bed_url}"
    curl -L -o "{output.lcr_bed}" "{params.lcr_bed_url}"
    """

rule download_cram_ref:
  output:
    fa = CRAM_REF,
    fai = CRAM_REF + ".fai",
  log:
    stdout = CRAM_REF + ".stdout",
    stderr = CRAM_REF + ".stderr",
  params:
    url = "https://broad-references.s3.amazonaws.com/hg38/v0/Homo_sapiens_assembly38.fasta"
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.fa}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}"

    curl -L "{params.url}" > "{output.fa}"
    curl -L "{params.url}.fai" > "{output.fai}"
    """

rule download_ref:
  input:
    samtools = config["tools"]["samtools"],
    sentieon = config["tools"]["sentieon"],
    rtgtools = config["tools"]["rtgtools"],
  output:
    fa = REFERENCE,
    fai = REFERENCE + ".fai",
    amb = REFERENCE + ".amb",
    ann = REFERENCE + ".ann",
    bwt = REFERENCE + ".bwt",
    pac = REFERENCE + ".pac",
    sa = REFERENCE + ".sa",
    sdf = directory(REFERENCE+ ".sdf"),
  log:
    stdout = REFERENCE + ".stdout",
    stderr = REFERENCE + ".stderr",
  params:
    url = "https://giab.s3.amazonaws.com/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz"
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.fa}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}" 

    curl -L "{params.url}" | gzip -dc > "{output.fa}"
    {input.samtools} faidx "{output.fa}"
    {input.sentieon} bwa index "{output.fa}"
    {input.rtgtools} format -o "{output.sdf}" "{output.fa}"
    """

rule autosomes_bed:
  input:
    rules.download_ref.output.fai,
  output:
    AUTOSOMES,
  log:
    stdout = AUTOSOMES + ".stdout",
    stderr = AUTOSOMES + ".stderr",
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}"

    cat "{input}" | head -n 25 | awk -v OFS='\t' '{{print $1,0,$2}}' > "{output}"
    """

rule download_sites_vcf:
  input:
    sentieon = config["tools"]["sentieon"],
  output:
    vcf = SITE_VCF,
    tbi = SITE_VCF + ".tbi",
  log:
    stdout = SITE_VCF + ".stdout",
    stderr = SITE_VCF + ".stderr",
  params:
    url = lambda wldc: config["input"]["known_sites"][wldc.vcf_name],
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.vcf}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}"
    
    url="{params.url}"
    pat="*.vcf.gz"
    if [[ $url == $pat ]]; then
      curl -L -o "{output.vcf}" "{params.url}"
      "{input.sentieon}" util vcfindex "{output.vcf}"
    else
      curl -L "{params.url}" | \
        "{input.sentieon}" util vcfconvert - "{output.vcf}"
    fi
    """

rule download_model:
  output:
    DNASCOPE_MODEL
  log:
    DNASCOPE_MODEL + ".log",
  params:
    url = lambda wldc: config["input"]["dnascope_model"][wldc.platform],
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output}")
    mkdir -p "$outdir"
    exec &>"{log}"

    curl -L -o "{output}" "{params.url}"
    """

rule dnascope_lr_script:
  output:
    script = DNASCOPE_LR_SCRIPT,
    model = DNASCOPE_LR_MODEL,
  params:
    url = config["input"]["dnascope_LR_package"],
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.script}")
    outdir2=$(dirname "$outdir")
    mkdir -p "$outdir2"

    cd "$outdir2"
    curl -L "{params.url}" | tar -zxf -
    """

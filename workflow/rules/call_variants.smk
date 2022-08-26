#!/usr/bin/env snakemake

"""
Call variants relative to a reference genome
"""

#### Output ####
HAPLOTYPER = "haplotyper/{sample}/calls.vcf.gz"
DNASCOPE_TMP = "dnascope_tmp/{sample}/calls.vcf.gz"
DNASCOPE = "dnascope/{sample}/calls.vcf.gz"
DNASCOPE_LR = "dnascope_lr/{sample}/calls.vcf.gz"

#### Functions ####
# Get the DNAscope model for the data type
def get_dnascope_model(wldc):
  platform = "Illumina"
  if "Ultima" in wldc.sample:
    platform = "Ultima"
  elif "Element" in wldc.sample:
    platform = "Element"
  return DNASCOPE_MODEL.format(platform=platform)

# Get the read filter for Ultima data
def get_dnascope_rf(wldc, input):
  if "Ultima" not in wldc.sample:
    return ""
  return f"--read_filter UltimaReadFilter"

# Get the proper alignemnt file for calling
def get_alignemnts(get_idx=False):
  def _get_alignemnts(wldc):
    if "Ultima" in wldc.sample:
      fn = CRAM.format(sample=wldc.sample)
      return fn + ".crai" if get_idx else fn
    else:
      fn = DEDUP.format(sample=wldc.sample)
      return fn + ".bai" if get_idx else fn
  return _get_alignemnts

#### Rules ####
rule haplotyper:
  input:
    bam = rules.dedup.output.bam,
    bai = rules.dedup.output.bai,
    qcal = rules.qualcal.output,
    fa = rules.download_ref.output.fa,
    fai = rules.download_ref.output.fai,
    autosomes_bed = rules.autosomes_bed.output,
    dbsnp = SITE_VCF.format(vcf_name="dbsnp"),
    sentieon = config["tools"]["sentieon"],
  output:
    vcf = HAPLOTYPER,
    tbi = HAPLOTYPER + ".tbi",
  benchmark:
    HAPLOTYPER + ".benchmark.txt"
  log:
    stdout = HAPLOTYPER + ".stdout",
    stderr = HAPLOTYPER + ".stderr",
  params:
    xargs = "--pcr_indel_model none"
  threads:
    192
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.vcf}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    {input.sentieon} driver -i "{input.bam}" -r "{input.fa}" -q "{input.qcal}" \
      --interval "{input.autosomes_bed}" -t {threads} --algo Haplotyper -d "{input.dbsnp}" {params.xargs} "{output.vcf}"
    """

rule dnascope:
  input:
    bam = get_alignemnts(),
    bai = get_alignemnts(True),
    fa = get_qualcal_ref(),
    fai = get_qualcal_ref(".fai"),
    autosomes_bed = rules.autosomes_bed.output,
    model = get_dnascope_model,
    dbsnp = SITE_VCF.format(vcf_name="dbsnp"),
    sentieon = config["tools"]["sentieon"],
  output:
    vcf = DNASCOPE_TMP,
    tbi = DNASCOPE_TMP + ".tbi",
  benchmark:
    DNASCOPE_TMP + ".benchmark.txt"
  log:
    stdout = DNASCOPE_TMP + ".stdout",
    stderr = DNASCOPE_TMP + ".stderr",
  params:
    xargs = "--pcr_indel_model none",
    read_filter = get_dnascope_rf,
  threads:
    192
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.vcf}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    {input.sentieon} driver -i "{input.bam}" -r "{input.fa}" {params.read_filter} \
      --interval "{input.autosomes_bed}" -t {threads} --algo DNAscope {params.xargs} \
      --model "{input.model}" -d "{input.dbsnp}" "{output.vcf}"
    """

rule dnamodelapply:
  input:
    vcf = rules.dnascope.output.vcf,
    tbi = rules.dnascope.output.tbi,
    fa = get_qualcal_ref(),
    fai = get_qualcal_ref(".fai"),
    model = get_dnascope_model,
    sentieon = config["tools"]["sentieon"],
  output:
    vcf = DNASCOPE,
    tbi = DNASCOPE + ".tbi",
  benchmark:
    DNASCOPE + ".benchmark.txt"
  log:
    stdout = DNASCOPE + ".stdout",
    stderr = DNASCOPE + ".stderr",
  threads:
    192
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.vcf}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    {input.sentieon} driver -t {threads} -r "{input.fa}" --algo DNAModelApply \
      --model "{input.model}" -v "{input.vcf}" "{output.vcf}"
    """

rule dnascope_longread:
  input:
    bam = rules.merge_pb.output.bam,
    bai = rules.merge_pb.output.bai,
    fa = rules.download_ref.output.fa,
    fai = rules.download_ref.output.fai,
    autosomes_bed = rules.autosomes_bed.output,
    model = DNASCOPE_LR_MODEL,
    dbsnp = SITE_VCF.format(vcf_name="dbsnp"),
    dnascope_hifi = DNASCOPE_LR_SCRIPT,
    sentieon = config["tools"]["sentieon"],
    bedtools = config["tools"]["bedtools"],
    bcftools = config["tools"]["bcftools"],
  output:
    vcf = DNASCOPE_LR,
    tbi = DNASCOPE_LR + ".tbi",
  benchmark:
    DNASCOPE_LR + ".benchmark.txt"
  log:
    stdout = DNASCOPE_LR + ".stdout",
    stderr = DNASCOPE_LR + ".stderr",
  threads:
    192
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.vcf}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    sentieon_dir=$(dirname "{input.sentieon}")
    bedtools_dir=$(dirname "{input.bedtools}")
    bcftools_dir=$(dirname "{input.bcftools}")

    export PATH=$sentieon_dir:$bedtools_dir:$bcftools_dir:$PATH
    bash "{input.dnascope_hifi}" -r "{input.fa}" -i "{input.bam}" -m "{input.model}" \
      -d "{input.dbsnp}" -- "{output.vcf}"
    """

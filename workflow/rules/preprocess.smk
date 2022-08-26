#!/usr/bin/env snakemake

"""
Pre-process samples through alignemnt, duplicate marking, and BQSR
"""

import multiprocessing

#### Output ####
MINIMAP2 = "ramdisk/minimap2/{sample}/{read_idx}/alignments.bam"
MERGE_PB = "ramdisk/merge_pb/{sample}/merged.bam"

N_SPLITS = "n_splits/n_splits.txt"
BWA = "ramdisk/bwa_mem/{sample}/{read_idx}/{n_subsets}/{subset}/alignments.bam"
LC = "ramdisk/dedup/lc/{sample}/score.txt.gz"
DEDUP = "ramdisk/dedup/dedup/{sample}/dedup.bam"
BQSR = "qualcal/{sample}/recal.table"

#### Functions ####
# Get the downsampled fastq file
def get_downsample_fq(pair):
  def _get_downsample_fq(wldc):
    return DOWNSAMPLE.format(pair=pair, **wldc.__dict__)
  return _get_downsample_fq

# Collect all of the BAMs output by BWA
def collect_bwa_bams(suffix=""):
  def _collect_bwa_bams(wldc):
    expected = []
    sample_list = config["input"]["samples"][wldc.sample]

    n_splits = 1
    with checkpoints.numa_splits.get().output[0].open() as fh:
      n_splits = int(fh.read().rstrip())

    for read_idx in range(len(sample_list)):
      for subset in range(n_splits):
        expected.append(BWA.format(subset=subset, n_subsets=n_splits, read_idx=read_idx, sample=wldc.sample) + suffix)
    return expected
  return _collect_bwa_bams

def get_qualcal_aln(wldc):
    return DEDUP.format(sample=wldc.sample)

def get_qualcal_ref(suffix=""):
  def _get_qualcal_ref(wldc):
    if "Ultima" in wldc.sample:
      return CRAM_REF + suffix
    else:
      return REFERENCE + suffix
  return _get_qualcal_ref

def get_known_vcfs(wldc):
  expected = []
  for vcf_name in config["input"]["known_sites"].keys():
    expected.append(SITE_VCF.format(vcf_name=vcf_name))
  return expected

def get_mm2_bam(suffix=""):
  def _get_mm2_bam(wldc):
    sample_fq = config["input"]["samples"][wldc.sample]
    n_fq = len(sample_fq)

    expected = []
    for i in range(len(sample_fq)):
      expected.append(MINIMAP2.format(sample=wldc.sample, read_idx=str(i)) + suffix)
    return expected
  return _get_mm2_bam

#### Rules ####
rule minimap2:
  input:
    fq = get_downsample_fq(0),
    fa = rules.download_ref.output.fa,
    sentieon = config["tools"]["sentieon"],
  output:
    bam = temp(MINIMAP2),
    bai = temp(MINIMAP2 + ".bai"),
  benchmark:
    MINIMAP2 + ".benchmark.txt"
  log:
    stdout = MINIMAP2 + ".stdout",
    stderr = MINIMAP2 + ".stderr",
  threads:
    192
  params:
    rg = lambda wldc: config["input"]["samples"][wldc.sample][int(wldc.read_idx)]["rg"],
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.bam}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    "{input.sentieon}" minimap2 -a -R '{params.rg}' -t {threads}  -x map-hifi \
      {input.fa} "{input.fq}" | \
      {input.sentieon} util sort -t {threads} -i - --sam2bam -o "{output.bam}"
    """

rule merge_pb:
  input:
    bam = get_mm2_bam(),
    bai = get_mm2_bam(".bai"),
    fa = rules.download_ref.output.fa,
    fai = rules.download_ref.output.fai,
    sentieon = config["tools"]["sentieon"],
  output:
    bam = temp(MERGE_PB),
    bai = temp(MERGE_PB + ".bai"),
  benchmark:
    MERGE_PB + ".benchmark.txt"
  log:
    stdout = MERGE_PB + ".stdout",
    stderr = MERGE_PB + ".stderr",
  threads:
    96
  params:
    input = lambda wldc, input: " -i '" + "' -i '".join(input.bam) + "'",
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.bam}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    {input.sentieon} driver {params.input} -t {threads} -r {input.fa} \
      --algo ReadWriter "{output.bam}"
    """

checkpoint numa_splits:
  output:
    N_SPLITS
  benchmark:
    N_SPLITS + ".benchmark.txt"
  log:
    N_SPLITS + ".log"
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output}")
    mkdir -p "$outdir"
    exec &>"{log}"

    lscpu | grep "NUMA node(s):" | sed 's/^NUMA node.* //' > "{output}"
    """

rule bwa_mem_subset:
  input:
    fq1 = get_downsample_fq(0),
    fq2 = get_downsample_fq(1),
    fa = rules.download_ref.output.fa,
    fai = rules.download_ref.output.fai,
    amb = REFERENCE + ".amb",
    ann = REFERENCE + ".ann",
    bwt = REFERENCE + ".bwt",
    pac = REFERENCE + ".pac",
    sa = REFERENCE + ".sa",
    sentieon = config["tools"]["sentieon"],
    sentieon_old = config["tools"]["sentieon_old"],
    igzip = config["tools"]["igzip"],
  output:
    bam = temp(BWA),
    bai = temp(BWA + ".bai"),
  benchmark:
    BWA + ".benchmark.txt"
  log:
    stdout = BWA + ".stdout",
    stderr = BWA + ".stderr",
    bwa = BWA + ".bwa.log",
    sort = BWA + ".sort.log",
  params:
    bwa_k = 100000000,
    rg = lambda wldc: config["input"]["samples"][wldc.sample][int(wldc.read_idx)]["rg"],
  threads:
    lambda wldc: int(multiprocessing.cpu_count() / int(wldc.n_subsets))
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.bam}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}"

    numa_cpulist=()
    n_splits=$(echo {wildcards.n_subsets})
    subset=$(echo {wildcards.subset})
    for i in $(seq 1 "$n_splits"); do
        i=$((i - 1))
        numa_cpulist+=($(lscpu | grep "NUMA node$i" | sed 's/^NUMA.* //'))
    done
    cpulist="${{numa_cpulist[$subset]}}"
    
    mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{{print $2}}')
    bwt_mem=$((mem_kb / 1024 / 1024 / {wildcards.n_subsets} - 6))
    export bwt_max_mem="$bwt_mem"G

    perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)';
    taskset -c "$cpulist" {input.sentieon} bwa mem -R "{params.rg}" \
      -K {params.bwa_k} -t {threads} -p "{input.fa}" \
      <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
        {input.sentieon} fqidx extract -F "${{subset}}"/"$n_splits" -K {params.bwa_k} \
        <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
          {input.igzip} -dc "{input.fq1}") \
        <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
          {input.igzip} -dc "{input.fq2}")) \
      2>"{log.bwa}" | \
      taskset -c  "$cpulist" {input.sentieon_old} util sort -t {threads} --sam2bam \
      -o "{output.bam}" -i - --bam_compression 1 2>"{log.sort}"
    """

rule locuscollector:
  input:
    bam = collect_bwa_bams(),
    bai = collect_bwa_bams(".bai"),
    fa = rules.download_ref.output.fa,
    fai = rules.download_ref.output.fai,
    sentieon = config["tools"]["sentieon"],
  output:
    LC,
  benchmark:
    LC + ".benchmark.txt"
  log:
    stdout = LC + ".stdout",
    stderr = LC + ".stderr",
  threads:
    192
  params:
    input_str = lambda wldc, input: '-i "' + '" -i "'.join([x for x in input.bam]) + '"',
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}"

    {input.sentieon} driver --traverse_param 500000/10000 -r {input.fa} \
      -t {threads} {params.input_str} --algo LocusCollector --fun score_info \
      "{output}"
    """

rule dedup:
  input:
    bam = collect_bwa_bams(),
    bai = collect_bwa_bams(".bai"),
    fa = rules.download_ref.output.fa,
    fai = rules.download_ref.output.fai,
    lc = rules.locuscollector.output,
    sentieon = config["tools"]["sentieon"],
  output:
    bam = temp(DEDUP),
    bai = temp(DEDUP + ".bai"),
  benchmark:
    DEDUP + ".benchmark.txt"
  log:
    stdout = DEDUP + ".stdout",
    stderr = DEDUP + ".stderr",
  threads:
    192
  params:
    input_str = lambda wldc, input: '-i "' + '" -i "'.join([x for x in input.bam]) + '"',
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.bam}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}"

    {input.sentieon} driver --traverse_param 500000/10000 -r {input.fa} \
      -t {threads} {params.input_str} --algo Dedup --rmdup \
      --score_info "{input.lc}" --bam_compression 1 "{output.bam}"
    """

rule qualcal:
  input:
    aln = get_qualcal_aln,
    known_vcfs = get_known_vcfs,
    fa = get_qualcal_ref(),
    fai = get_qualcal_ref(".fai"),
    autosomes_bed = rules.autosomes_bed.output,
    sentieon = config["tools"]["sentieon"],
  output:
    BQSR,
  benchmark:
    BQSR + ".benchmark.txt"
  log:
    stdout = BQSR + ".stdout",
    stderr = BQSR + ".stderr",
  threads:
    192
  params:
    vcf_str = lambda wldc, input: '-k "' + '" -k "'.join([x for x in input.known_vcfs]) + '"',
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}"

    {input.sentieon} driver --traverse_param 500000/10000 -r {input.fa} \
      --interval {input.autosomes_bed} -t {threads} -i "{input.aln}" \
      --algo QualCal {params.vcf_str} "{output}"
    """

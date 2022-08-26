#!/usr/bin/env snakemake

"""
Downsample input datasets to a standard coverage of bases
"""

#### Output ####
COUNT_BASES = "count_bases/{sample}/count.txt"
DOWNSAMPLE_SEED = "downsample_seed/{sample}/{read_idx}/seed.txt"
DOWNSAMPLE = "downsample/{sample}/{read_idx}/{pair}/reads.fq.gz"

#### Functions ####
def get_sample_fq(wldc):
  expected = []

  sample_list = config["input"]["samples"][wldc.sample]
  for read_idx in range(len(sample_list)):
    fq_list = sample_list[read_idx]["fq"]
    for pair in range(len(fq_list)):
      expected.append(FASTQ.format(sample=wldc.sample, **locals()))

  return expected

#### Rules ####
rule count_bases:
  input:
    fq = get_sample_fq,
    pigz = ancient(config["tools"]["pigz"]),
  output:
    COUNT_BASES
  benchmark:
    COUNT_BASES + ".benchmark.txt"
  log:
    stdout = COUNT_BASES + ".stdout",
    stderr = COUNT_BASES + ".stderr",
  params:
    fq_str = lambda wldc, input: '"' + '" "'.join([x for x in input.fq]) + '"'
  threads:
    8
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    {input.pigz} -p {threads} -dc {params.fq_str} | paste - - - - | cut -f 2 | tr -d '\n' | wc -c > "{output}"
    """

rule downsample_seed:
  output:
    DOWNSAMPLE_SEED
  benchmark:
    DOWNSAMPLE_SEED + ".benchmark.txt"
  log:
    stdout = DOWNSAMPLE_SEED + ".stdout",
    stderr = DOWNSAMPLE_SEED + ".stderr",
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}"

    echo $RANDOM > "{output}"
    """

rule downsample:
  input:
    fq = rules.download_samples.output.fq,
    base_count = ancient(rules.count_bases.output),
    seed = ancient(rules.downsample_seed.output),
    seqtk = ancient(config["tools"]["seqtk"]),
    pigz = ancient(config["tools"]["pigz"]),
  output:
    DOWNSAMPLE,
  benchmark:
    DOWNSAMPLE + ".benchmark.txt"
  log:
    stdout = DOWNSAMPLE + ".stdout",
    stderr = DOWNSAMPLE + ".stderr",
  threads:
    8
  params:
    target_bases = config["params"]["n_target_bases"],
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    n_counted=$(cat "{input.base_count}")
    seed=$(cat "{input.seed}")
    n_target="{params.target_bases}"
    frac=$(bc -l <<< "$n_target / $n_counted")

    {input.seqtk} sample -s "$seed" "{input.fq}" 0"$frac" | \
      {input.pigz} -p {threads} -c > "{output}"
    """

#!/usr/bin/env snakemake

"""
Evaluation of VCF files against the GIAB truthset
"""

#### Output ####
HAPPY = "eval/{pipeline}/{sample}/{hcr_region}/sample"

#### Functions ####
def get_calls_vcf(suffix=""):
  def _get_calls_vcf(wldc):
    if wldc.pipeline == "DNAseq":
      return HAPLOTYPER.format(sample=wldc.sample) + suffix
    elif "HiFi" in wldc.sample:
      return DNASCOPE_LR.format(sample=wldc.sample) + suffix
    else:
      return DNASCOPE.format(sample=wldc.sample) + suffix
  return _get_calls_vcf

def get_lcr_bed(wldc):
  if wldc.hcr_region == "Ultima-HCR":
    return [TRUTHSET + "_ultima-lcr.bed"]
  else:
    return []

#### Rules ####
rule happy:
  input:
    calls_vcf = get_calls_vcf(),
    calls_tbi = get_calls_vcf(".tbi"),
    truth_vcf = rules.download_truthset.output.vcf,
    truth_tbi = rules.download_truthset.output.vcf_idx,
    truth_bed = rules.download_truthset.output.bed,
    lcr_bed = get_lcr_bed,
    ref = rules.download_ref.output.fa,
    sdf = rules.download_ref.output.sdf,
    rtgtools = config["tools"]["rtgtools"],
    python2 = config["tools"]["python2"],
    happy = config["tools"]["hap.py"],
    bedtools = config["tools"]["bedtools"],
  output:
    summary = HAPPY + ".summary.csv",
    extended = HAPPY + ".extended.csv",
    vcf = HAPPY + ".vcf.gz",
    tbi = HAPPY + ".vcf.gz.tbi",
  benchmark:
    HAPPY + ".benchmark.txt"
  log:
    stdout = HAPPY + ".stdout",
    stderr = HAPPY + ".stderr",
  threads:
    8
  resources:
    mem_mb=39000,
  params:
    happy_xargs = "--no-decompose --no-leftshift",
    lcr_bed = lambda wldc, input: input.lcr_bed[0] if len(input.lcr_bed) > 0 else "",
  shell:
    """
    set -exvo pipefail
    outdir=$(dirname "{output.summary}")
    mkdir -p "$outdir"
    exec >"{log.stdout}" 2>"{log.stderr}"
    # Link the input file into the output directory and cd to outdir. Hap.py 
    # cannot handle complex path characters (', =, ' ', etc.)
    in_vcf=$(realpath "{input.calls_vcf}")
    rm "$outdir"/calls.vcf.gz "$outdir"/calls.vcf.gz.tbi || $(exit 0)
    ln -s "$in_vcf" "$outdir"/calls.vcf.gz
    ln -s "$in_vcf".tbi "$outdir"/calls.vcf.gz.tbi
    outpre=$(basename "{output.summary}")
    outpre="${{outpre%%.summary.csv}}"

    eval_bed=$(realpath "{input.truth_bed}")
    ref_path=$(realpath "{input.ref}")
    sdf_path=$(realpath "{input.sdf}")
    truth_path=$(realpath "{input.truth_vcf}")
    lcr_path=$(pwd)/"{params.lcr_bed}"
    cd "$outdir"
    rtg_dir=$(dirname "{input.rtgtools}")
    export PATH=$rtg_dir:$PATH

    # Generate a HCR bed
    lcr_bed="{params.lcr_bed}"
    if [[ -n "$lcr_bed" ]]; then
      {input.bedtools} subtract -a "$eval_bed" -b "$lcr_path" > HG002_hcr.bed
      eval_bed="HG002_hcr.bed"
    fi

    # Evaluate
    {input.python2} {input.happy} "$truth_path" \
      calls.vcf.gz --verbose -o "$outpre" -r "$ref_path" \
      -f "$eval_bed" --threads \
      {threads} --engine=vcfeval --engine-vcfeval-template "$sdf_path" \
      {params.happy_xargs}
    """

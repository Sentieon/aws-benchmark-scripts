#!/usr/bin/env snakemake

"""
Perform alignment, preprocessing, variant calling and evaluation of different pipelines
"""

include: "rules/download.smk"
include: "rules/downsample.smk"
include: "rules/preprocess.smk"
include: "rules/call_variants.smk"
include: "rules/evaluate.smk"

#### Output ####
DOWNLOAD_PREPROCESS_ALL = "all_downloaded.txt"
ALL_PIPELINES = "run_all.txt"
HC_PIPELINE = "run_hc.txt"

#### Functions ####
def get_fq_all_samples(wldc):
  expected = []
  for sample in config["input"]["samples"].keys():
    if "Ultima" not in sample:
      expected.append(DOWNSAMPLE.format(sample=sample))
  return expected

def get_ultima_cram(wldc):
  expected = []
  for sample in config["input"]["samples"].keys():
    if "Ultima" in sample:
      expected.append(DOWNSAMPLE_CRAM.format(sample=sample))
  return expected

def get_all_vcfs(wldc):
  expected = []
  for vcf_name in config["input"]["known_sites"].keys():
    expected.append(SITE_VCF.format(vcf_name=vcf_name))
  return expected

def get_all_models(wldc):
  expected = []
  for platform in config["input"]["dnascope_model"].keys():
    expected.append(DNASCOPE_MODEL.format(platform=platform))
  return expected

def get_dnaseq_eval(wldc):
  expected = []
  sample_dict = config["input"]["samples"]
  for sample in sample_dict.keys():
    if "Illumina" in sample:
      expected.append(HAPPY.format(pipeline="DNAseq", sample=sample))
  return expected

def get_dnascope_eval(wldc):
  expected = []
  sample_dict = config["input"]["samples"]
  for sample in sample_dict.keys():
    expected.append(HAPPY.format(pipeline="DNAscope", sample=sample))
  return expected

#### Rules ####
rule download_all:
  input:
    get_fq_all_samples,
    get_ultima_cram,
    TRUTHSET + ".vcf.gz",
    REFERENCE,
    AUTOSOMES,
    get_all_vcfs,
    get_all_models,
    DNASCOPE_LR_MODEL,
  output:
    DOWNLOAD_PREPROCESS_ALL,
  shell:
    """ touch "{output}" """

rule run_hc:
  input:
    get_dnaseq_eval,
  output:
    HC_PIPELINE,
  shell:
    """ touch "{output}" """

rule run_all:
  input:
    HC_PIPELINE,
    get_dnascope_eval,
  output:
    ALL_PIPELINES,
  shell:
    """ touch "{output}" """

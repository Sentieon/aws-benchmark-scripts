#!/usr/bin/env python3

"""
Collect the DNAseq runtimes from one or more benchmark runs
"""

import argparse
import collections
import glob
import os.path
import sys

pipelines = {
    "DNAseq": collections.OrderedDict([
        (
            "Alignment",
            {
                "path": "ramdisk/bwa_mem/{sample}/0/*/*/alignments.bam.benchmark.txt",
                "reduce_func": max,
            },
        ),
        (
            "LocusCollector",
            {
                "path": "ramdisk/dedup/lc/{sample}/score.txt.gz.benchmark.txt",
            },
        ),
        (
            "Dedup",
            {
                "path": "ramdisk/dedup/dedup/{sample}/dedup.bam.benchmark.txt",
            },
        ),
        (
            "QualCal",
            {
                "path": "qualcal/{sample}/recal.table.benchmark.txt",
            },
        ),
        (
            "Haplotyper",
            {
                "path": "haplotyper/{sample}/calls.vcf.gz.benchmark.txt",
            }
        ),
    ]),
    "DNAscope CRAM": collections.OrderedDict([
        (
            "DNAscope_tmp",
            {
                "path": "dnascope_tmp/{sample}/calls.vcf.gz.benchmark.txt",
            },
        ),
        (
            "DNAscope",
            {
                "path": "dnascope/{sample}/calls.vcf.gz.benchmark.txt",
            },
        ),
    ]),
    "DNAscope": collections.OrderedDict([
        (
            "Alignment",
            {
                "path": "ramdisk/bwa_mem/{sample}/0/*/*/alignments.bam.benchmark.txt",
                "reduce_func": max,
            },
        ),
        (
            "LocusCollector",
            {
                "path": "ramdisk/dedup/lc/{sample}/score.txt.gz.benchmark.txt",
            },
        ),
        (
            "Dedup",
            {
                "path": "ramdisk/dedup/dedup/{sample}/dedup.bam.benchmark.txt",
            },
        ),
        (
            "DNAscope_tmp",
            {
                "path": "dnascope_tmp/{sample}/calls.vcf.gz.benchmark.txt",
            },
        ),
        (
            "DNAscope",
            {
                "path": "dnascope/{sample}/calls.vcf.gz.benchmark.txt",
            },
        ),
    ]),
    "DNAscope LongRead": collections.OrderedDict([
        (
            "Alignment",
            {
                "path": "ramdisk/minimap2/{sample}/*/alignments.bam.benchmark.txt",
            },
        ),
        (
            "Merge",
            {
                "path": "ramdisk/merge_pb/{sample}/merged.bam.benchmark.txt",
            },
        ),
        (
            "DNAscope-LR",
            {
                "path": "dnascope_lr/{sample}/calls.vcf.gz.benchmark.txt",
            },
        ),
    ])
}

pipeline_samples = {
    "DNAseq": [
        "Illumina_HG002_HiSeqX_40x",
        "Illumina_HG002_NovaSeq_40x",
    ],
    "DNAscope": [
        "Illumina_HG002_HiSeqX_40x",
        "Illumina_HG002_NovaSeq_40x",
        "Element_HG002_100x",
    ],
    "DNAscope CRAM": [
        "Ultima_HG002_cram",
    ],
    "DNAscope LongRead": [
        "PacBio_HG002_HiFi_Chem2",
    ],
}

row_order = [
    "Alignment", "Merge", "LocusCollector", "Dedup", "QualCal", "Haplotyper",
    "DNAscope_tmp", "DNAscope", "DNAscope-LR", "Total (s)", "Total (min)",
    "$/hr", "Total compute cost"
]

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input_directory", nargs="+", help="A directory containing benchmark results")
    parser.add_argument("--instance_price", nargs="+", help="The instance price in $/hr")
    parser.add_argument("--instance_type", nargs="+", help="The instance type")
    parser.add_argument("--outfile", default=sys.stdout, type=argparse.FileType('w'))
    return parser.parse_args(argv)

def bm_file_to_sec(bm_path):
    """ Extract the runtime in seconds from the snakemake benchmark file """
    fh = open(bm_path)
    _hdr = fh.readline()
    values = fh.readline().rstrip().split('\t')
    return float(values[0])

def main(args):
    assert len(args.input_directory) == len(args.instance_price)
    assert len(args.input_directory) == len(args.instance_type)

    # Iterate over each instance type
    overall_results = collections.OrderedDict()
    for in_dir, price, instance_type in zip(args.input_directory, args.instance_price, args.instance_type):
        price = float(price)
        for pipeline, samples in pipeline_samples.items():
            for sample in samples:
                runtime_in_sec = collections.OrderedDict()
                missing_stage = False
                for stage, d in pipelines[pipeline].items():
                    path = d["path"].format(sample=sample)
                    fullpath = os.path.join(in_dir, path)
                    stage_results = glob.glob(fullpath)
                    if not stage_results:
                        missing_stage = True
                        print(f"Missing files for instance '{instance_type}', pipeline '{pipeline}', and sample '{sample}' at: {fullpath}", file=sys.stderr, flush=True)
                        break

                    reduce_func = d.get("reduce_func", sum)
                    runtime_in_sec[stage] = reduce_func([bm_file_to_sec(x) for x in stage_results])
                if missing_stage:
                    print(f"Skipping pipeline '{pipeline}' for instance '{instance_type}' and sample '{sample}'")
                    continue

                # summary results
                runtime_in_sec["Total (s)"] = sum(runtime_in_sec.values())
                runtime_in_sec["Total (min)"] = runtime_in_sec["Total (s)"] / 60
                runtime_in_sec["$/hr"] = price
                runtime_in_sec["Total compute cost"] = runtime_in_sec["Total (s)"] / 60 / 60 * price

                key = " -- ".join((pipeline, sample))
                if instance_type not in overall_results:
                    overall_results[instance_type] = collections.OrderedDict()
                overall_results[instance_type][key] = runtime_in_sec.copy()

    # Output a nice table
    row = [[k] * len(v) for k,v in overall_results.items()]
    row = ["Instance type"] + [item for sublist in row for item in sublist]  # flatten
    print('\t'.join(row), file=args.outfile)

    row = ["Pipeline"] + [key for v in overall_results.values() for key in v.keys()]
    print('\t'.join(row), file=args.outfile)

    for row_name in row_order:
        row = [row_name]
        for _, v in overall_results.items():
            for _, d in v.items():
                row.append(d.get(row_name, 0.0))
        print('\t'.join([str(x) for x in row]), file=args.outfile)

    return 0

if __name__ == "__main__":
    args = parse_args()
    main(args)

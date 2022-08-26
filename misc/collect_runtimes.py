#!/usr/bin/env python3

"""
Collect the DNAseq runtimes from one or more benchmark runs
"""

import argparse
import collections
import glob
import os.path
import sys

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
        runtime_in_sec = collections.OrderedDict()
        price = float(price)

        # Alignment runtimes
        alignment_results = glob.glob(os.path.join(in_dir, "ramdisk/bwa_mem/Illumina_HG002_NovaSeq_40x/0/*/*/alignments.bam.benchmark.txt"))
        runtime_in_sec["Alignment"] = max([bm_file_to_sec(x) for x in alignment_results])  # One or more alignemnt jobs run in parallel. All jobs need to finish before the next step, so total runtime is the max of all runtimes

        # Other stages
        runtime_in_sec["LocusCollector"] = bm_file_to_sec(os.path.join(in_dir, "ramdisk/dedup/lc/Illumina_HG002_NovaSeq_40x/score.txt.gz.benchmark.txt"))
        runtime_in_sec["Dedup"] = bm_file_to_sec(os.path.join(in_dir, "ramdisk/dedup/dedup/Illumina_HG002_NovaSeq_40x/dedup.bam.benchmark.txt"))
        runtime_in_sec["QualCal"] = bm_file_to_sec(os.path.join(in_dir, "qualcal/Illumina_HG002_NovaSeq_40x/recal.table.benchmark.txt"))
        runtime_in_sec["Haplotyper"] = bm_file_to_sec(os.path.join(in_dir, "haplotyper/Illumina_HG002_NovaSeq_40x/calls.vcf.gz.benchmark.txt"))

        # summary results
        results = collections.OrderedDict()
        results["Total (s)"] = sum(runtime_in_sec.values())
        results["Total (min)"] = sum(runtime_in_sec.values()) / 60
        results["$/hr"] = price
        results["Total compute cost"] = sum(runtime_in_sec.values()) / 60 / 60 * price

        # Collect everything together
        runtime_in_sec.update(results)
        overall_results[instance_type] = runtime_in_sec.copy()

    # Output a nice table
    row = ["Instance type"] + args.instance_type
    print('\t'.join(row), file=args.outfile)
    row_names = list(overall_results[args.instance_type[0]].keys())
    for row_name in row_names:
        row = [row_name] + [str(overall_results[x][row_name]) for x in args.instance_type]
        print('\t'.join(row), file=args.outfile)

    return 0

if __name__ == "__main__":
    args = parse_args()
    main(args)

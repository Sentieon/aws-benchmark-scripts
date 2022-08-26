#!/usr/bin/env bash

set -exvuo pipefail

DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
WORK_DIR="/home/ec2-user/work"
USE_RAMDISK=false  # Switched to `true` on the hpc6a.48xlarge instance

if [[ "$USE_RAMDISK" == true ]]; then
    # Increase the size of /dev/shm
    sudo mount -o remount,size=350G /dev/shm
    ln -s /dev/shm "$WORK_DIR"/ramdisk
fi

mkdir -p "$WORK_DIR/ramdisk/tmp"

# An initial run generated these files. If they are present, we can copy them
# from the FSX mount into the working directory to save time downloading and
# pre-processing the input files

## Link files in FsX into the work directory
#FSX_PATH=(
#    /fsx/benchmark_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
#    /fsx/benchmark_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.amb
#    /fsx/benchmark_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.ann
#    /fsx/benchmark_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.bwt
#    /fsx/benchmark_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai
#    /fsx/benchmark_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.pac
#    /fsx/benchmark_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.sa
#    /fsx/benchmark_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.sdf
#    /fsx/benchmark_data/input_data/fastq/ElementBio_wgs_HG002_pcr-free_100x_R1.fastq.gz
#    /fsx/benchmark_data/input_data/fastq/ElementBio_wgs_HG002_pcr-free_100x_R2.fastq.gz
#    /fsx/benchmark_data/input_data/fastq/HG002.hiseqx.pcr-free.40x.R1.fastq.gz
#    /fsx/benchmark_data/input_data/fastq/HG002.hiseqx.pcr-free.40x.R2.fastq.gz
#    /fsx/benchmark_data/input_data/fastq/HG002.novaseq.pcr-free.40x.R1.fastq.gz
#    /fsx/benchmark_data/input_data/fastq/HG002.novaseq.pcr-free.40x.R2.fastq.gz
#    /fsx/benchmark_data/dnascope_models/SentieonDNAscopeModelElementBio0.3.model
#    /fsx/benchmark_data/dnascope_models/DNAscopeHiFiBeta0.4.pipeline/DNAscopeHiFiBeta0.4.model
#    /fsx/benchmark_data/dnascope_models/DNAscopeHiFiBeta0.4.pipeline/dnascope_HiFi.sh
#    /fsx/benchmark_data/dnascope_models/DNAscopeHiFiBeta0.4.pipeline/hapgt.py
#    /fsx/benchmark_data/dnascope_models/DNAscopeHiFiBeta0.4.pipeline/vcf_mod.py
#    /fsx/benchmark_data/dnascope_models/SentieonDNAscopeModelUltima0.2.model
#    /fsx/benchmark_data/dnascope_models/SentieonDNAscopeModel1.1.model
#    /fsx/benchmark_data/truthset/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
#    /fsx/benchmark_data/truthset/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
#    /fsx/benchmark_data/truthset/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
#    /fsx/benchmark_data/truthset/strat
#    /fsx/benchmark_data/vcfs/Homo_sapiens_assembly38.dbsnp138.vcf.gz
#    /fsx/benchmark_data/vcfs/Homo_sapiens_assembly38.known_indels.vcf.gz
#    /fsx/benchmark_data/vcfs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#    /fsx/benchmark_data/vcfs/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi
#    /fsx/benchmark_data/vcfs/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
#    /fsx/benchmark_data/vcfs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
#    /fsx/benchmark_data/input_data/fastq/m64011_190830_220126.fastq.gz
#    /fsx/benchmark_data/input_data/fastq/m64011_190901_095311.fastq.gz
#    /fsx/benchmark_data/input_data/fastq/m64012_190920_173625.fastq.gz
#    /fsx/benchmark_data/input_data/fastq/m64012_190921_234837.fastq.gz
#    /fsx/benchmark_data/input_data/fastq/m64015_190920_185703.fastq.gz
#    /fsx/benchmark_data/input_data/fastq/m64015_190922_010918.fastq.gz
#    /fsx/benchmark_data/cram/005401-UGAv3-1-CACATCCTGCATGTGAT.cram
#    /fsx/benchmark_data/cram/005401-UGAv3-1-CACATCCTGCATGTGAT.cram.crai
#    /fsx/benchmark_data/truthset/ultima-lcr.bed
#    /fsx/benchmark_data/count_bases
#    /fsx/benchmark_data/count_bases_cram
#    /fsx/benchmark_data/downsample_seed
#    /fsx/benchmark_data/downsample
#    /fsx/benchmark_data/downsample_cram
#)
#WORKDIR_PATH=(
#    "$WORK_DIR"/"reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
#    "$WORK_DIR"/"reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.amb"
#    "$WORK_DIR"/"reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.ann"
#    "$WORK_DIR"/"reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.bwt"
#    "$WORK_DIR"/"reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai"
#    "$WORK_DIR"/"reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.pac"
#    "$WORK_DIR"/"reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.sa"
#    "$WORK_DIR"/"reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.sdf"
#    "$WORK_DIR"/"input_fq/Element_HG002_100x/0/0/data.fq.gz"
#    "$WORK_DIR"/"input_fq/Element_HG002_100x/0/1/data.fq.gz"
#    "$WORK_DIR"/"input_fq/Illumina_HG002_HiSeqX_40x/0/0/data.fq.gz"
#    "$WORK_DIR"/"input_fq/Illumina_HG002_HiSeqX_40x/0/1/data.fq.gz"
#    "$WORK_DIR"/"input_fq/Illumina_HG002_NovaSeq_40x/0/0/data.fq.gz"
#    "$WORK_DIR"/"input_fq/Illumina_HG002_NovaSeq_40x/0/1/data.fq.gz"
#    "$WORK_DIR"/"dnascope_model/Element/dnascope.model"
#    "$WORK_DIR"/"dnascope_lr_script/DNAscopeHiFiBeta0.4.pipeline/DNAscopeHiFiBeta0.4.model"
#    "$WORK_DIR"/"dnascope_lr_script/DNAscopeHiFiBeta0.4.pipeline/dnascope_HiFi.sh"
#    "$WORK_DIR"/"dnascope_lr_script/DNAscopeHiFiBeta0.4.pipeline/hapgt.py"
#    "$WORK_DIR"/"dnascope_lr_script/DNAscopeHiFiBeta0.4.pipeline/vcf_mod.py"
#    "$WORK_DIR"/"dnascope_model/Ultima/dnascope.model"
#    "$WORK_DIR"/"dnascope_model/Illumina/dnascope.model"
#    "$WORK_DIR"/"truthset/v4.2.1/HG002/truthset.bed"
#    "$WORK_DIR"/"truthset/v4.2.1/HG002/truthset.vcf.gz"
#    "$WORK_DIR"/"truthset/v4.2.1/HG002/truthset.vcf.gz.tbi"
#    "$WORK_DIR"/"truthset/v4.2.1/HG002/truthset/strat"
#    "$WORK_DIR"/"sites_vcfs/dbsnp/calls.vcf.gz"
#    "$WORK_DIR"/"sites_vcfs/known_indels/calls.vcf.gz"
#    "$WORK_DIR"/"sites_vcfs/mills/calls.vcf.gz"
#    "$WORK_DIR"/"sites_vcfs/dbsnp/calls.vcf.gz.tbi"
#    "$WORK_DIR"/"sites_vcfs/known_indels/calls.vcf.gz.tbi"
#    "$WORK_DIR"/"sites_vcfs/mills/calls.vcf.gz.tbi"
#    "$WORK_DIR"/"input_fq/PacBio_HG002_HiFi_Chem2/0/0/data.fq.gz"
#    "$WORK_DIR"/"input_fq/PacBio_HG002_HiFi_Chem2/1/0/data.fq.gz"
#    "$WORK_DIR"/"input_fq/PacBio_HG002_HiFi_Chem2/2/0/data.fq.gz"
#    "$WORK_DIR"/"input_fq/PacBio_HG002_HiFi_Chem2/3/0/data.fq.gz"
#    "$WORK_DIR"/"input_fq/PacBio_HG002_HiFi_Chem2/4/0/data.fq.gz"
#    "$WORK_DIR"/"input_fq/PacBio_HG002_HiFi_Chem2/5/0/data.fq.gz"
#    "$WORK_DIR"/"input_cram/Ultima_HG002_cram/ultima.cram"
#    "$WORK_DIR"/"input_cram/Ultima_HG002_cram/ultima.cram.crai"
#    "$WORK_DIR"/"truthset/v4.2.1/HG002/truthset_ultima-lcr.bed"
#    "$WORK_DIR"/
#    "$WORK_DIR"/
#    "$WORK_DIR"/
#    "$WORK_DIR"/
#    "$WORK_DIR"/
#)
#
#for i in $(seq 1 ${#FSX_PATH[@]}); do
#    i=$((i - 1))
#    fsx_path="${FSX_PATH[$i]}"
#    workdir_path="${WORKDIR_PATH[$i]}"
#    workdir_base=$(dirname "$workdir_path")
#    mkdir -p "$workdir_base"
#    ln -s "$fsx_path" "$workdir_path"
#done

snakefile="$DIR"/../workflow/Snakefile.smk

mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')
mem_mb=$((mem_kb / 1024))
export LD_PRELOAD=/home/ec2-user/programs/jemalloc/jemalloc-5.2.1/lib/libjemalloc.so.2
export MALLOC_CONF=metadata_thp:auto,background_thread:true,dirty_decay_ms:30000,muzzy_decay_ms:30000
export SENTIEON_TMPDIR="$WORK_DIR/ramdisk/tmp"
export SENTIEON_LICENSE=  # FQDN:port or localhost file for the Sentieon license

architechure=$(uname -m)
configfile="$DIR"/config/config.yaml
if [[ "$architechure" == "x86_64" ]]; then
    configfile="$DIR"/config/config.yaml
elif [[ "$architechure" == "arm" || "$architechure" == "aarch64" ]]; then
    configfile="$DIR"/config/config_arm.yaml
else
    echo "UNKNOWN ARCHITECHURE"
    exit 1
fi

targets=(
  "eval/DNAseq/Illumina_HG002_HiSeqX_40x"
  "eval/DNAseq/Illumina_HG002_NovaSeq_40x"
  "eval/DNAscope/Illumina_HG002_HiSeqX_40x"
  "eval/DNAscope/Illumina_HG002_NovaSeq_40x"
  "eval/DNAscope/Element_HG002_100x"
  "eval/DNAscope/Ultima_HG002_cram"
  "eval/DNAscope_LongRead/PacBio_HG002_HiFi_Chem2"
)

for target in "${targets[@]}"; do
    expected=()
    for hcr_region in "none" "Ultima-HCR"; do
        expected+=("$target"/"$hcr_region"/sample.summary.csv)
    done
    $snakemake --reason --resources mem_mb="$mem_mb" --configfile "$configfile" -s "$snakefile" -j $(nproc) -d "$WORK_DIR" -p --verbose --keep-going ${expected[@]}
done

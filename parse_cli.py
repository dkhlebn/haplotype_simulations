#!/home/dkhlebnikov/miniconda3/envs/PofO_Sim/bin/python3

import argparse
import pandas as pd

from pathlib import Path

def chrom_lengths(fname):
    """Get katyotype chromosome lengths"""
    dt = dict()
    with open(fname, "r") as fin:
        for line in fin:
            chrom, leng = line.strip().split("\t")
            dt[chrom] = int(leng)
    return dt
    

def hap_orig_mapping(fname):
    """Get imprinting data for haplotypes"""
    dt = dict()
    with open(fname, "r") as fin:
        for line in fin:
            hap, mapping = line.strip().split("\t")
            dt[hap] = mapping
    dt = pd.DataFrame(dt.items(), columns=["Haplotype", "Parent"])
    dt.sort_values("Haplotype").reset_index().iloc[:, 1:]
    return dt


def get_popul_data(fname):
    popul_data = pd.read_excel(
        "support_data/20130606_sample_info.xlsx", sheet_name="HD Genotypes"
    ).iloc[:, [0, 1]]
    return popul_data

def parse_arguments():
    parser = argparse.ArgumentParser(description="Simulation.")

    parser.add_argument('-s', '--step', type=str, required=True, help="Step number, technical")
    parser.add_argument('-a', '--agg_dir', type=Path, required=True, help="Directory for results aggregation")
    parser.add_argument('-vcf', '--vcf_path', type=Path, default="joint_vcfs", help="Path to merged 1000 Genomes and haplotypes VCFs")
    parser.add_argument('-p', '--pops', type=Path, required=True, help="1000Genome Population information location")
    parser.add_argument('-ip', '--imprinting_labels', type=Path, required=True, help="Path to file with imprinting labels")
    parser.add_argument('-chrom', '--chrom_sizes', type=Path, required=True, help="Path to file with chromosome sizes")
    parser.add_argument('-wd', '--workdir', type=Path, default="temp_dir", help="Working directory base")
    parser.add_argument('-m', '--mode', type=str, required=True, help="Modality of the analysis")

    args = parser.parse_args()
    args.perm_flag, args.hide_flag, args.perm_type, args.dist = args.mode.split('_')
    args.hide_flag = args.hide_flag == "Hide"
    args.perm_flag = args.perm_flag == "Permute"
    args.chr_dt = chrom_lengths(args.chrom_sizes)
    args.wd = args.workdir / args.mode / args.step
    args.imprinting_data = hap_orig_mapping(args.imprinting_labels)
    args.pops = get_popul_data(args.pops)
    return args

ARGS = parse_arguments()
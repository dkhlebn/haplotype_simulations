#!/home/dkhlebnikov/miniconda3/envs/PofO_Sim/bin/python3


import pandas as pd
from parse_cli import parse_arguments
from utils import create_directories, remove_subdirs
from permutation_funcs import run_step

ARGS = parse_arguments()
result = None

while result is None:
  create_directories(ARGS.wd, ARGS.chr_dt)
  result = run_step(ARGS.chr_dt)
  remove_subdirs(ARGS.wd)
result.loc[:, ["Haplotype", "Cluster"]].to_csv(f"{ARGS.agg_dir}/step_{ARGS.step_n}.tsv", sep="\t")

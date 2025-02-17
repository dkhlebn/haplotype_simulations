import pandas as pd
from parse_cli import parse_arguments
from utils import create_directories, remove_subdirs
from permutation_funcs import run_step

ARGS = parse_arguments()
result = None

while result is None:
  create_directories(ARGS.wd, ARGS.chr_dt)
  result = run_step(chr_dt, wd_name)
  remove_subdirs(ARGS.wd)
result.loc[:, ["Haplotype", "Cluster"]].to_csv(f"{ARGS.agg_dir}/step_{args.step_n}.tsv", sep="\t")
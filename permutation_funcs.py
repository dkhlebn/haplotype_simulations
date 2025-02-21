#!/home/dkhlebnikov/miniconda3/envs/PofO_Sim/bin/python3

import subprocess as sp
import pickle
import shlex as sh
from collections import ChainMap
import numpy as np
import pandas as pd
import pysam as ps
from pqdm.threads import pqdm as pqdm_t
from pqdm.processes import pqdm as pqdm_p
from sklearn.cluster import AgglomerativeClustering
from .parse_cli import ARGS
from .utils import process_snp, create_tabix_index, get_distances, classical_mds
from .randomization_utils import generate_chromosome_sectors, permute_imprinting


def create_mock_haps(chrom, chrom_sectors):
    """Create mock merged 1000G + Haps file in
    corresponding directory. Returns filename
    for the mock file"""
    merged = ps.VariantFile(f"{ARGS.vcf_path}/{chrom}/1kg.haps.{chrom}.merged.vcf.gz")
    mock_f = ps.VariantFile(
        f"{ARGS.wd}/merged/{chrom}/1kg_haps.merged.vcf.gz", "wz", header=merged.header
    )
    for i, (start, end) in enumerate(chrom_sectors):
        flag = i % 2 #np.random.random() > 0.5
        for record in merged.fetch(chrom, start, end):
            if flag:
                record = process_snp(record)
            mock_f.write(record)
    merged.close(), mock_f.close()
    create_tabix_index(f"{ARGS.wd}/merged/{chrom}/1kg_haps.merged.vcf.gz")
    return f"{ARGS.wd}/merged/{chrom}/1kg_haps.merged.vcf.gz"


def run_plink_cmds(chrom, merged_path):
    """PLINK wrapper to generate distamce matrices for 1000G
    specimen and haplotype pseudo-specimen"""
    pdata_dir=f"{ARGS.wd}/plink_data"
    pdist_dir=f"{ARGS.wd}/plink_dist"
    cmd = f"sh haplotype_simulations/plink_script.sh {chrom} {merged_path} {pdata_dir} {pdist_dir}"
    proc = sp.Popen(sh.split(cmd))
    _, _ = proc.communicate()
    return f"{pdist_dir}/{chrom}/donor2.{chrom}.distances"


def execute_data_svd(plink_fs_prefix, n_cmp=5):
    """Run PCA on the distance matrix"""
    distances = pd.read_csv(f"{plink_fs_prefix}.mdist", header=None, sep=" ").iloc[:, :-1]
    sample_ids = pd.read_csv(f"{plink_fs_prefix}.mdist.id", sep='\t', header=None, names=["famids", "IID"])
    distances -= distances.values.mean()
    _, _, vh = np.linalg.svd(distances, full_matrices=False) # tweak False if result is bad
    principal_components = vh.T[:, :n_cmp]

    pca_df = pd.DataFrame(principal_components, columns=[f"PC{i+1}" for i in range(n_cmp)])
    pca_df.index = sample_ids["famids"]
    res = pca_df.merge(ARGS.pops, left_index=True, right_on="Sample").set_index("Sample")
    return res


def perform_MDS(plink_fs_prefix, n_cmp=5):
    """Run classical MDS on obtained distance matrix"""
    distances = pd.read_csv(f"{plink_fs_prefix}.mdist", header=None, sep=" ").iloc[:, :-1]
    sample_ids = pd.read_csv(f"{plink_fs_prefix}.mdist.id", sep='\t',
                             header=None, names=["famids", "IID"])
    pca_coords = classical_mds(distances, k=n_cmp)
    pca_df = pd.concat([sample_ids, pca_coords], axis=1)
    pca_df = pca_df.merge(ARGS.pops, left_on="famids", right_on="Sample").iloc[:, 2:]
    pca_df = pca_df.set_index('Sample')
    return pca_df


def distances_on_chrom(chrom, pca_df, n_cmp=5):
    """Perform distance calculation on PCA data and cluster haps into two groups"""
    pc_list = [f"PC{i}" for i in range(1, n_cmp + 1)]
    hap_df = (pca_df.query("Population == 'Our_Haplotype'")
              .reset_index()[['Sample', *pc_list]])
    haps = hap_df.set_index('Sample').apply(tuple, axis=1).to_dict()
    return get_distances(haps, pca_df, chrom, ARGS.dist.lower(), n_cmp=n_cmp)


def run_mocking(chrom, sectors):
    """Wrapper func to run simulation in parallel
    Runs on processes"""
    mock_path = create_mock_haps(chrom, sectors)
    plink_pref = run_plink_cmds(chrom, mock_path)
    return (chrom, plink_pref)


def run_dimreduction(chrom, plink_pref):
    """Wrapper func to run simulation in parallel
    Runs on threads"""
    df_chrom = perform_MDS(plink_pref)
    dt = distances_on_chrom(chrom, df_chrom)
    return dt


def assign_parents(step_df):
    """Infer mock haplotypes PofO"""
    clusts = AgglomerativeClustering().fit(step_df)
    df = pd.DataFrame({"Haplotype": step_df.index, "Cluster": clusts.labels_}).merge(
        ARGS.imprinting_data, how="left"
    )
    df = permute_imprinting(df)
    mapping = df.groupby("Cluster")["Parent"].agg(lambda x: x.mode().iloc[0])
    df["Cluster"] = df["Cluster"].map(mapping)
    if df["Cluster"].nunique() == 1:
        return None
    return df


def run_sim_step(chr_dt):
    """Run simulation step"""
    sectors = generate_chromosome_sectors(chr_dt)
    args = [{"chrom": chrom, "sectors": sectors}
                         for chrom, sectors in sectors.items()]
    chr2pref = pqdm_p(args, run_mocking, n_jobs=22, argument_type="kwargs")
    args = [{"chrom": t[0], "plink_pref": t[1]} for t in chr2pref]
    step_results = pqdm_t(
        args, run_dimreduction, n_jobs=22, argument_type="kwargs", disable=True
    )
    merged_dict = dict(ChainMap(*reversed(step_results)))
    step_result = pd.DataFrame(merged_dict).transpose()
    return assign_parents(step_result)


def compute_real_data(path_to_mvcfs):
    """Run computations on real data"""
    chr_list = [f'chr{i}' for i in [*range(1, 23)]]
    create_directories(path_to_mvcfs, "./", chr_list)
    args = [{"chrom": chrom, 
             "merged_path": f"{path_to_mvcfs}/{chrom}/1kg.haps.{chrom}.merged.vcf.gz"} for chrom in chr_list]
    chr2pref = pqdm_p(args, run_plink_cmds, n_jobs=22, argument_type="kwargs")
    ddts = []
    for chrom, ppref in tqdm(zip(chr_list, chr2pref)):
        df_chrom = perform_MDS(ppref)
        dt = distances_on_chrom(chrom, df_chrom)
        ddts.append(dt)
    merged_dict = dict(ChainMap(*ddts))
    step_result = pd.DataFrame(merged_dict).transpose()
    remove_subdirs(f"{path_to_mvcfs}/plink_data")
    remove_subdirs(f"{path_to_mvcfs}/plink_dist")
    return step_result

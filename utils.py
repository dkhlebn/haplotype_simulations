#!/home/dkhlebnikov/miniconda3/envs/PofO_Sim/bin/python3

import os
import shutil
import shlex as sh
import subprocess as sp
import pandas as pd
from scipy.spatial import distance
from scipy.linalg import inv
from numpy.linalg import eigh

# File and Directories manipulation
def create_directories(base_dir, res_dir, chromosomes):
    """Initiate dirs in working directory"""
    subdirs = ["merged", "plink_data", "plink_dist"]
    os.makedirs(base_dir, exist_ok=True)
    os.makedirs(res_dir, exist_ok=True)
    for subdir in subdirs:
        path = os.path.join(base_dir, subdir)
        os.makedirs(path, exist_ok=True)
        for chrom in chromosomes:
            nested_path = os.path.join(path, chrom)
            os.makedirs(nested_path, exist_ok=True)


def remove_subdirs(base_dir):
    """Remove directories after success"""
    if not os.path.exists(base_dir):
        raise FileNotFoundError(f"{base_dir} does not exist.")
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path):
            shutil.rmtree(item_path)


def chrom_lengths(fname):
    """Get katyotype chromosome lengths"""
    dt = {}
    with open(fname, "r") as fin:
        for line in fin:
            chrom, leng = line.strip().split("\t")
            dt[chrom] = int(leng)
    return dt


def hap_orig_mapping(fname):
    """Get imprinting data for haplotypes"""
    dt = {}
    with open(fname, "r") as fin:
        for line in fin:
            hap, mapping = line.strip().split("\t")
            dt[hap] = mapping
    dt = pd.DataFrame(dt.items(), columns=["Haplotype", "Parent"])
    return dt


def create_tabix_index(fname):
    """Wrapper to run tabix to index vcf files"""
    tabix = "/home/dkhlebnikov/tools/bin/tabix"
    cmd = f"{tabix} -f -p vcf {fname}"
    proc = sp.Popen(sh.split(cmd))
    proc.communicate()
    return 0

# processing utils
def process_snp(r):
    """Helper to flip value of SNP GT field"""
    if r.samples["HAP1"]["GT"] != (0, 0) or r.samples["HAP2"]["GT"] != (0, 0):
        r.samples["HAP1"]["GT"], r.samples["HAP2"]["GT"] = (
            r.samples["HAP2"]["GT"][::-1],
            r.samples["HAP1"]["GT"][::-1],
        )
    return r


def get_centroids(popul_df, n_cmp=5):
    """Calculate centroids"""
    df = (popul_df[popul_df['Population'] != "Our_Haplotype"]
                         .groupby('Population')
                         .mean()
                         .iloc[:, :n_cmp]
                         .apply(tuple, axis=1)
                         .to_dict())
    return df


def classical_mds(d, k=5):
    n, x = d.shape[0], d ** 2
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ x @ H
    eigvals, eigvecs = eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]
    pos_eigvals = eigvals[:k]
    pos_eigvals = pos_eigvals[pos_eigvals > 0]
    evecs = eigvecs[:, :len(pos_eigvals)]

    points = evecs * np.sqrt(pos_eigvals)
    return pd.DataFrame(points, columns=[f"PC{i}" for i in range(1, k+1)])


def calc_dist(crds1, crds2, cov_inv, mode="euclidean"):
    """Helper func to calculate vector distances"""
    if mode == "mahalanobis":
        if cov_inv is None:
            raise ValueError("Mahalanobis distance requires a covariance matrix.")
        return distance.mahalanobis(crds1, crds2, cov_inv)
    #mode: euclidean, cityblock, cosine
    return getattr(distance, mode)(crds1, crds2)


def get_distances(haps, pca_df, chrom, mode, n_cmp=5):
    """Get distances of haps to centroids"""
    centroids = get_centroids(pca_df, n_cmp)
    haps = {h: coords[:n_cmp] for h, coords in haps.items()}
    result_dt = {f"{chrom}_{hap}": {} for hap in haps}
    for popul in centroids:
        cov_matrix = (pca_df[pca_df['Population'] == popul]
                      .iloc[:, :n_cmp]
                      .cov()
                      .values
                     )
        cov_inv = inv(cov_matrix)
        for hap in haps:
            result_dt[f"{chrom}_{hap}"][popul] = calc_dist(haps[hap],
                                                           centroids[popul],
                                                           cov_inv,
                                                           mode=mode)
    return result_dt

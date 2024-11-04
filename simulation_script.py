import pysam as ps
import numpy as np
import pandas as pd
import subprocess as sp
import shlex as sh
import os
import glob

from pqdm.processes import pqdm
from itertools import repeat
from sklearn.decomposition import PCA


HEAD_DIR = "temp_sim_dir"


def init_dirs(head_dir_name):
    """Create directories for simulations"""
    if not os.path.exists(head_dir_name):
        os.mkdir(head_dir_name)
    if len(os.listdir(head_dir_name)) < 3:
        os.mkdir(f"{head_dir_name}/merged")
        os.mkdir(f"{head_dir_name}/plink_data")
        os.mkdir(f"{head_dir_name}/plink_dist")


def generate_chromosome_sectors(chromosome_lengths, min_size=100_000):
    """Generate chromosome sectors for exchanging haplotype SNPs"""
    chromosome_sectors = {}
    for chrom, length in chromosome_lengths.items():
        sectors = []
        n_dots = length // min_size
        borders = np.sort(np.random.randint(1, length, size=n_dots))
        start = 0
        for b in borders:
            sectors.append((start, b))
            start = b + 1
        chromosome_sectors[chrom] = sectors
    return chromosome_sectors


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
    return dt


def merge_hap_1kg(chrom, hap1_path, hap2_path, outdir=f"{HEAD_DIR}/merged"):
    """BCFTools merge wrapper to merge 1000G and donor/mock haps data"""
    if not os.path.exists(f"{outdir}/{chrom}"):
        os.mkdir(f"{outdir}/{chrom}")

    outfile = f"{outdir}/{chrom}/1kg_haps.merged.vcf"
    bcftools = "/home/dkhlebn/tools/bcftools"
    cmd = f"{bcftools} merge 1kg_vcfs/1kg_{chrom}.vcf.gz {hap1_path} {hap2_path} --missing-to-ref -Oz -o {outfile}"
    proc = sp.Popen(sh.split(cmd))
    _, stderr = proc.communicate()
    return outfile


def create_tabix_index(fname):
    """Wrapper to run tabix to index vcf files"""
    tabix = "/home/dkhlebn/tools/tabix"
    cmd = f"{tabix} -f -p vcf {fname}"
    proc = sp.Popen(sh.split(cmd))
    proc.communicate()
    return


def process_snp(r):
    """Helper to flip value of SNP GT field"""
    if r.samples["HAP1"]["GT"] != (0, 0) or r.samples["HAP2"]["GT"] != (0, 0):
        # worst python code line ever written
        r.samples["HAP1"]["GT"], r.samples["HAP2"]["GT"] = (
            r.samples["HAP2"]["GT"][::-1],
            r.samples["HAP1"]["GT"][::-1],
        )
    return r


def create_mock_haps(chrom, chrom_sectors, outdir=f"{HEAD_DIR}/merged"):
    """Create mock merged 1000G + Haps file in
    corresponding directory. Returns filename
    for the mock file"""
    if not os.path.exists(f"{outdir}/{chrom}"):
        os.mkdir(f"{outdir}/{chrom}")

    merged = ps.VariantFile(f"joint_vcfs/{chrom}/1kg.haps.{chrom}.merged.vcf.gz")
    mock_f = ps.VariantFile(
        f"{outdir}/{chrom}/1kg_haps.merged.vcf.gz", "w", header=merged.header
    )

    for i, (start, end) in enumerate(chrom_sectors):
        flag = i % 2 == 0
        for record in merged.fetch(chrom, start, end):
            if flag:
                record = process_snp(record)
            mock_f.write(record)
    merged.close(), mock_f.close()

    create_tabix_index(f"{outdir}/{chrom}/1kg_haps.merged.vcf.gz")
    return f"{outdir}/{chrom}/1kg_haps.merged.vcf.gz"


def run_plink_cmds(
    chrom,
    merged_path,
    pdata_dir=f"{HEAD_DIR}/plink_data",
    pdist_dir=f"{HEAD_DIR}/plink_dist",
):
    """PLINK wrapper to generate distamce matrices for 1000G
    specimen and haplotype pseudo-specimen"""
    if not os.path.exists(f"{pdata_dir}/{chrom}") and not os.path.exists(
        f"{pdist_dir}/{chrom}"
    ):
        os.mkdir(f"{pdata_dir}/{chrom}")
        os.mkdir(f"{pdist_dir}/{chrom}")

    cmd = f"sh ./plink_script.sh {chrom} {merged_path} {pdata_dir} {pdist_dir}"
    proc = sp.Popen(sh.split(cmd))
    _, _ = proc.communicate()
    return f"{pdist_dir}/{chrom}/donor2.{chrom}.distances"


def execute_data_pca(plink_fs_prefix, n_cmp=5):
    """Perform a dimensionality reduction on PLINK distance matrices"""
    distances = pd.read_table(f"{plink_fs_prefix}.mdist", header=None, sep=" ").iloc[
        :, :-1
    ]
    sample_ids = pd.read_table(f"{plink_fs_prefix}.mdist.id", header=None).rename(
        columns={0: "famids", 1: "IID"}
    )
    pca_data = pd.DataFrame(PCA(n_components=n_cmp).fit_transform(distances)).rename(
        columns={i: f"pc{i + 1}" for i in range(n_cmp)}
    )
    pca_data.index = sample_ids.famids
    return pca_data.merge(POPS_DATA, left_index=True, right_on="Sample").set_index(
        "Sample"
    )


def calc_dist(crds1, crds2, mode="euclid"):
    """Helper func to calculate vector distances"""

    def _calc_euclid(crds1, crds2):
        s = 0
        for i, j in zip(crds1, crds2):
            s += (i - j) ** 2
        s = s ** (1 / 2)
        return s

    def _calc_manh(crds1, crds2):
        s = 0
        for i, j in zip(crds1, crds2):
            s += np.abs(i - j)
        return s

    if mode == "euclid":
        return _calc_euclid(crds1, crds2)
    if mode == "manhattan":
        return _calc_mahn(crds1, crds2)
    raise ValueError("Wrong option for distance!")


def cluster_on_chrom(chrom, pca_df):
    """Perform distance calculation on PCA data and cluster haps into two groups"""
    haps = pca_df[pca_df.Population == "Our_Haplotype"].iloc[:, :5]
    haps = {x[0]: x[1:3] for x in haps.to_records().tolist()}
    result_dt = dict()
    centroids = (
        pca_df[pca_df.Population != "Our_Haplotype"]
        .reset_index()
        .iloc[:, 1:]
        .groupby("Population")
        .mean()
    )
    centroids = {x[0]: x[1:3] for x in centroids.to_records().tolist()}
    for hap in haps:
        result_dt[chrom + "_" + hap] = {
            popul: calc_dist(haps[hap], centroids[popul]) for popul in centroids
        }
    return result_dt


def run_clustering(chrom, sectors):
    """Wrapper func to run simulation in parallel"""
    mock_path = create_mock_haps(chrom, sectors)
    plink_pref = run_plink_cmds(chrom, mock_path)
    df_chrom = execute_data_pca(plink_pref)
    dt = cluster_on_chrom(chrom, df_chrom)
    return dt


HEAD_DIR = "temp_sim_dir"
init_dirs(HEAD_DIR)


POPS_DATA = pd.read_excel(
    "support_data/20130606_sample_info.xlsx", sheet_name="HD Genotypes"
).iloc[:, [0, 1]]
IMPRINT_DATA = hap_orig_mapping("support_data/data_imprinting.tsv")
chr_dt = chrom_lengths("support_data/hg38.chr.genome.tab")
chromosomes = list(chr_dt.keys())


# simulation part
sim_dt = dict()
for i in range(1):
    chrom_dt = dict()
    sectors = generate_chromosome_sectors(chr_dt)
    args = [x for x in zip(chromosomes, sectors.values(), repeat(chrom_dt))]
    args = [
        {key: value for key, value in zip(["chrom", "sectors", "chrom_dt"], arg_l)}
        for arg_l in args
    ]
    step_results = pqdm(args, run_clustering, n_jobs=22, argument_type="kwargs")
    dt = step_results[0]
    for res_dict in step_results[1:]:
        dt.update(res_dict)
    step_result = pd.DataFrame(chrom_dt).transpose()
    clusts = AgglomerativeClustering().fit(step_result)

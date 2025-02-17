#!/home/dkhlebnikov/miniconda3/envs/PofO_Sim/bin/python3

import numpy as np
from parse_cli import ARGS

def poisson_disk_1d(range_min, range_max, min_distance, num_samples):
    """Generates a sample of intervals borders for chromosome segmentation"""
    cell_size = min_distance
    num_cells = int((range_max - range_min) / cell_size) + 1
    grid = [None] * num_cells

    points, attempts = [], 0
    while len(points) < num_samples and attempts < 10 * num_samples:
        candidate = int(np.random.uniform(range_min, range_max))
        cell_index = int((candidate - range_min) / cell_size)
        if cell_index < 0 or cell_index >= num_cells:
            continue
        valid = True
        for neighbor in range(max(0, cell_index - 1), min(num_cells, cell_index + 2)):
            if grid[neighbor] is not None and abs(candidate - grid[neighbor]) < min_distance:
                valid = False
                break
        if valid:
            points.append(candidate)
            grid[cell_index] = candidate
        attempts += 1
    return sorted(points)


def generate_chromosome_sectors(chromosome_lengths, min_size=30_000, sect_size=100_000):
    """Generate chromosome sectors for exchanging haplotype SNPs"""
    chromosome_sectors = {}
    for chrom, length in chromosome_lengths.items():
        borders = poisson_disk_1d(0, length, min_size, length // sect_size)
        sectors, start = [], 0
        for b in borders:
            sectors.append((start, b))
            start = b + 1
        chromosome_sectors[chrom] = sectors
    return chromosome_sectors


def hide_labels(df, ratio=0.3):
    """Hides a portion of imprinting labels if necessary"""
    mask = df['Parent'].isin(['F', 'M'])
    n = int(mask.sum() * ratio)
    random_indices = np.random.choice(df[mask].index, size=n, replace=False)
    df.loc[random_indices, 'Parent'] = np.nan
    return df


def permute_labels(df, mode):
    """Permutes imprinting labels"""
    def _pair_permute(df):
        should_permute = np.random.rand(len(df) // 2) < 0.5
        for i, permute in enumerate(should_permute):
            hap1_idx, hap2_idx = 2 * i, 2 * i + 1
            if permute:
                df.loc[[hap1_idx, hap2_idx], 'Parent'] = (df
                                                     .loc[[hap2_idx, hap1_idx], 'Parent'].values)
        return df

    def _all_permute(df):
        df_imp = df.copy()
        df_imp["Parent"] = np.random.permutation(df_imp["Parent"])
        return df_imp

    def _total_permute(df):
        mask = df['Parent'].isin(['F', 'M'])
        indices = df[mask].index
        shuffled_labels = df.loc[indices, 'Parent'].values
        np.random.shuffle(shuffled_labels)
        df.loc[indices, 'Parent'] = shuffled_labels
        return df

    if mode == "Pair":
        return _pair_permute(df)
    if mode == "Total":
        return _total_permute(df)
    if mode == "All":
        return _all_permute(df)
    return df


def permute_imprinting(df):
    """Randomizes imprinting labels if necessary"""
    if ARGS.perm_flag:
        df = permute_labels(df, ARGS.perm_type)
    if ARGS.hide_flag:
        df = hide_labels(df)
    return df

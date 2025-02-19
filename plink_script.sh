#!/bin/bash


chrom=$1;
merged_path=$2;
pdata_dir=$3;
pdist_dir=$4;

PLINK="/home/dkhlebnikov/tools/PLINK/plink";
$PLINK --vcf ${merged_path} --out ${pdata_dir}/${chrom}/ --silent;
$PLINK --bfile ${pdata_dir}/${chrom}/ --autosome --geno 0.1 --mind 0.1 --maf 0.05 --nonfounders --allow-no-sex --recode --out ${pdist_dir}/${chrom}/ --silent;
$PLINK --allow-no-sex --nonfounders --file ${pdist_dir}/${chrom}/ --distance-matrix --out ${pdist_dir}/${chrom}/donor2.${chrom}.distances --silent;

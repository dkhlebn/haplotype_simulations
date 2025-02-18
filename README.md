## PofO inference simulations model
This repo is for the code used for statistical analysis of PofO assignment in my project.
We utilize creation of mock haplotypes by mixing the haplotype data we have and mapping 
them onto the 1000 Genomes data.

There are a few modes we can get the results with: 
* Permute or not permute true imprinting labels (`${PERM} = Perm / NoPerm`);
* Hide or not hide some of the true imprinting labels (`${HIDE} = Hide / NoHide`);
* if permuted, how do we permute labels: in pairs (`${MODE} = Pair`), between known labels (`${MODE} = Total`) or between all labels including the unknown (`${MODE} = All`);
* Distance can also be tweaked, accepted options are `${DIST} = Euclidean / Cityblock / Cosine / Mahalanobis`.

### Run example:
To run on example data, run:
```
python3 haplotype_simulations/main.py -s example_step -a ./results -p example/20130606_sample_info.xlsx -ip example/data_imprinting.tsv -chrom example/hg38.chr.genome.tab -wd ./temp_dir -m Perm_NoHide_Total_Euclidean
```

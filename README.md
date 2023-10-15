# CisTrans-ECAS
This procedure integrates the association of cis- and trans-components of gene expression with the phenotype or expression traits (e-traits), to prioritize causal genes and construct regulatory networks.


## Dependencies
CisTrans-ECAS is mainly tested in R-3.5.1. It depends on the following software:
+ [PLINK](https://www.cog-genomics.org/plink/): Used to extract a subset of genotype data from binary genotype files.
+ [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview): Used to construct prediction models for gene expression.

It also depends on the following R packages:
+ [snpStats](https://bioconductor.org/packages/3.8/bioc/html/snpStats.html): Used to read genotype data from binary genotype files.
+ [cpgen](https://rdrr.io/cran/cpgen/man/cGWAS.emmax.html): Used to fit linear mixed models for the association analysis of traits or e-traits with the cis- or trans-components of gene expression.
+ [fdrtool](https://cran.r-project.org/web/packages/fdrtool/index.html): Estimates false discovery rates (FDR).

## Applications
###1. Decomposition of gene expression into cis- and trans-components
#### gcta_cis.R
This script constructs prediction models for gene expression based on SNPs in the vicinity of genes, depending on GCTA.
##### Usage: 
```
$ Rscript gcta_cis.R [options]
```
##### Options:
```
  --gffs_file=value   An RData file containing gene annotation information, including a data frame named 'gffs' with columns Gene, chr, start, end, strand, Qsymbols, and Note.
  --exp_file=value   An RData file containing a gene expression matrix with a matrix named 'exp_m', where row names represent gene IDs and column names represent individual IDs.
  --genodir=value   Directory containing genotype data in plink's binary format, including .bed, .bim, and .fam files, saved per chromosome.
  --gfile_prefix=value   Pattern for the file names of genotype data, where chromosome number is represented by %s, exclusive of suffix.
  --plinkdir=value   Path to the plink program.
  --gctadir=value   Path to the GCTA program.
  --out_dir=value   Output directory.
  --extend=value   The range of cis-regulatory regions around genes. SNPs within a distance less than 'extend' from genes are used to build gene expression prediction models. e.g., extend=100000.
  --ncor=value   The number of threads to use.
  --help   Display this help message.
```
##### Example:
```
$ cd ~/test/
$ Rscript  gcta_cis.R  --gffs_file=./test_data/gffs.RData  --exp_file=./test_data/exp_matrix.RData  --genodir=./test_data/  --gfile_prefix=%s_529_test  --plinkdir=plink  --gctadir=gcta  --out_dir=./test_res/gcta_cis  --extend=1e5  --ncor=2
```

#### merge_cis_res.R
This script is used to merge the results generated by gcta_cis.R for different chromosomes.
##### Usage: 
```
$ Rscript merge_cis_res.R [options]
```
##### Options:
  --file_dir=value   Directory of output results for gcta_cis.R
  --help         Display this help message.
##### Example:
```
$ Rscript merge_cis_res.R --file_dir=./test_res/gcta_cis
```

###2. Associating phenotypes with cis- and trans-components of gene expression to prioritize candidate genes
#### cistrans_twas.R
This script is used to separately associate phenotypes with the cis- and trans-components of gene expression to identify candidate genes that influence the phenotype.
##### Usage: 
```
$ Rscript cistrans_twas.R [options]
```
##### Options:
```
  --gffs_file=value   An RData file containing gene annotation information, including a data frame named 'gffs' with columns Gene, chr, start, end, strand, Qsymbols, and Note.
  --K_file=value   An RData file containing kinship matrix, including a matrix named 'K'.
  --pheno_file=value   An RData file containing phenotype data, including a data frame named 'pheno' with row names as individual IDs and column names as traits.
  --vc_file=value   An RData file generated by 'merge_cis_res.R' that includes the variance explained by cis-genetic variants for each gene and the corresponding p-value.
  --cis_file=value   An RData file generated by 'merge_cis_res.R' that includes a matrix of the cis-component of gene expression, named 'res_cis_m', with genes as rows and individuals as columns.
  --trans_file=value   An RData file generated by 'merge_cis_res.R' that includes a matrix of the trans-component of gene expression, named 'res_trans_m', with genes as rows and individuals as columns.
  --out_dir=value    Output directory.
  --p_threshold=value   P-value threshold for cis-genetic variance.
  --help         Display this help message.
```
##### Example:
```
$ Rscript cistrans_twas.R   --gffs_file=./test_data/gffs.RData  --K_file=./test_data/Kinship.RData  --pheno_file=./test_data/pheno.RData  --vc_file=./test_res/gcta_cis/gcta_cis_vc.RData  --cis_file=./test_res/gcta_cis/gcta_cis_cis.RData  --trans_file=./test_res/gcta_cis/gcta_cis_trans.RData  --out_dir=./test_res/cistrans_twas  --p_threshold=1.62e-6 
```

###3. Associating e-traits with cis- and trans-components of gene expression to identify gene upstream regulatory factors
#### cistrans_etwas.R
This script is used to separately associate e-traits (gene expression levels) with the cis- and trans-components of gene expression to identify upstream regulatory factors that influence e-trait expression.
##### Usage: 
```
$ Rscript cistrans_etwas.R [options]
```
##### Options:
```
  --gffs_file=value   An RData file containing gene annotation information, including a data frame named 'gffs' with columns Gene, chr, start, end, strand, Qsymbols, and Note.
  --K_file=value   An RData file containing kinship matrix, including a matrix named 'K'.
  --exp_file=value An RData file containing a gene expression matrix with a matrix named 'exp_m', where row names represent gene IDs and column names represent individual IDs.
  --gene_file=value An RData file, which contains a vector named 'genes,' including e-traits used for cistrans_etwas.
  --vc_file=value   An RData file generated by 'merge_cis_res.R' that includes the variance explained by cis-genetic variants for each gene and the corresponding p-value.
  --cis_file=value   An RData file generated by 'merge_cis_res.R' that includes a matrix of the cis-component of gene expression, named 'res_cis_m', with genes as rows and individuals as columns.
  --trans_file=value   An RData file generated by 'merge_cis_res.R' that includes a matrix of the trans-component of gene expression, named 'res_trans_m', with genes as rows and individuals as columns.
  --out_dir=value    Output directory.
  --p_threshold=value   P-value threshold for cis-genetic variance.
  --help         Display this help message.
```
##### Example:
```
$  Rscript cistrans_etwas.R   --gffs_file=./test_data/gffs.RData  --K_file=./test_data/Kinship.RData  --exp_file=./test_data/exp_matrix.RData  --gene_file=./test_data/genes.RData  --vc_file=./test_res/gcta_cis/gcta_cis_vc.RData  --cis_file=./test_res/gcta_cis/gcta_cis_cis.RData  --trans_file=./test_res/gcta_cis/gcta_cis_trans.RData  --out_dir=./test_res/cistrans_etwas  --p_threshold=1.62e-6 
```

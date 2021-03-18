Example for trait-trajectory analysis
================

This is an example of the analysis performed for a dataset of pancreatic
islet development from [Byrnes et
al](https://doi.org/10.1038/s41467-018-06176-3).

The Monocle object is available to download from figshare.com:
[link](https://figshare.com/articles/dataset/Monocle_Objects_-_V2_Dataset/6783554?backTo=/collections/Lineage_dynamics_of_murine_pancreatic_development_at_single-cell_resolution/4158458).

Load the following required packages:

``` r
library(monocle)
library(Seurat)
library(data.table)
```

Source the functions from this github directory:

``` r
source('functions_GWAS_traj.R')
```

### Preprocessing

First, we extract the expression matrix and metadata, including the
pseudotime of the cells. We Load the Monocole object:

``` r
load('E14_fev_lineage_monocle_ob.Rdata')
```

We extract the count matrix and cell metadata from this object:

``` r
# count matrix
exp <- HSMM_seur_var@assayData$exprs

# metadata
meta <- HSMM_seur_var@phenoData@data
```

Might be a good idea to remove the monocle object and save the new
objects:

``` r
rm('HSMM_seur_var')
saveRDS(object = exp, file = "panc.exp.RDS")
saveRDS(object = meta, file = "panc.meta.RDS")
```

#### Processing the count matrix

We keep genes expressed in at least 10 cells

``` r
genes.keep <- rowSums(as.matrix(exp) > 0) >= 10
exp <- exp[genes.keep,]
```

We normalize the data using Seurat’s log-normalization:

``` r
exp <- Seurat::NormalizeData(exp)
```

The following converts mouse genes to human orthologues:

``` r
exp <- conv_hs(exp)
```

#### Calculating cell-trait association scores

This part involves running MAGMA’s gene property analysis for each cell
separately. This might take long, so consider using multiple cores.
First, we generate ‘covariate files’ for MAGMA for each cell. Each
output file is a table with the columns: gene, normalized expression in
the cell, and average normalized expression in the dataset. The function
create a folder “cov\_files”, and write the files into it.

``` r
generate_covs_cells(exp = exp, cors = 20, wd = ".")
```

The next step, it running the gene property analysis, for each cell.
This analysis firs the following the regression model to each cell:

![](https://github.com/eldadshulman/scGWAS/blob/master/data/pic/eq.PNG)

where Z is the vector of the gene’s Z-score converted from the p-values
obtained from MAGMA’s gene analysis for the trait. B is a matrix of
technical confounders, including gene length and SNPs’ LDs, calculated
by MAGMA. C is the vector of normalized expression of the genes in the
cell, and A is a vector of the average normalized expression for the
gene in the dataset. The t-statistic of \(\beta_c\) is taken as the
score for the association between the cell and trait.

if(\!(dir.exists(“gene\_prop”))) dir.create(“gene\_prop”) This function
will Note that this step takes even longer.

``` r
gene_prop_cells(cor = 30)
```

Last, extract the cell scores from MAGMA’s output, namely the
t-statistics.

``` r
tss <- get_scores(cors = 20)
```

#### Calculating trajectory-trait association

We now examine the correlation between cell-trait association scores and
the pseudotime assigned to the cells. First, alanlyzing the trajectory
without considering branch identity, we fit the following linear
regression model:

Cell Score \~ Pseudotime

and perform a one-sided test (Pseudotime coefficient \> 0). Significance
indicates that cells’ association with the trait increases along the
trajectory’s pseudotime.

For trajectories with more than one branch (pancreas, neurons, B cells,
kidney PTA, and heart datasets), to examine if the trait association was
trait trait-dependent, we used a procedure similar to Monocle’s method
for identifying branch-dependent genes (Census) \[14\]. Briefly, we
compared the goodness of fit of the following two linear regression
models: Cell Score \~ Pseudotime Cell Score \~
Pseudotime+Branch:Pseudotime where Where the first model assumes that
the trait is not branch dependent, and the second model, which contains
the interaction term between branch and pseudotime (Branch:Pseudotime),
assumes branch dependence (Fig. S1). Likelihood ratio test was carried
out using the R package VGAM \[46\] (lrtest function). As this approach
requires the assignment of branch identity to each cell, including
progenitor the cells before the branch point, we followed Monocle’s
approach: We divided the unbranched progenitors into branches by
ordering them according to pseudotime, and assigning odd and even ranked
cells to the first and second branches, respectively. We assigned the
first progenitor to both branches. (Similarly, in the case of three
branches, progenitors were split according to ranking, assigning cells
with rankings that follow the arithmetic sequences 21+23n,32+23n,43+23n
(n=0,1,2,3, … )…) to the first, second, and third branches,
respectively). Furthermore, to examine if there is a positive
correlation between trait-association scores and pseudotime within a
branch, we used the following model (Fig. S1C): Cell Score \~
Branch:Pseudotime Note that this model differs from the interaction
model used for testing branch dependency, which includes both pseudotime
and the interaction term, and therefore, does not have a coefficient
that directly indicates for the direction of the correlation.

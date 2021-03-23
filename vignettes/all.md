Genetic mapping of developmental trajectories for complex traits and
diseases - vignette
================
Eldad David Shulman & Ran Elkon

-   [Step 1: Identification of connections between developmental
    trajectories and
    traits](#step-1-identification-of-connections-between-developmental-trajectories-and-traits)
    -   [1a. Converting GWAS variant scores into gene-trait association
        scores](#1a-converting-gwas-variant-scores-into-gene-trait-association-scores)
    -   [1b. Calculating cell-trait association
        scores](#1b-calculating-cell-trait-association-scores)
    -   [1c. Trajectory inference.](#1c-trajectory-inference)
    -   [1d. Examining the association between trait and
        trajectory](#1d-examining-the-association-between-trait-and-trajectory)
-   [Step 2: Elucidate molecular pathways that underlie the link between
    the trajectory and
    trait](#step-2-elucidate-molecular-pathways-that-underlie-the-link-between-the-trajectory-and-trait)
    -   [2a. Finding pathways enriched in the
        trajectory](#2a-finding-pathways-enriched-in-the-trajectory)
    -   [2b. Examining the trajectory-enriched pathwaysfor trait
        association](#2b-examining-the-trajectory-enriched-pathwaysfor-trait-association)
-   [Step 3: Prioritize genes that carry the link between the pathway,
    trait and
    trajectory](#step-3-prioritize-genes-that-carry-the-link-between-the-pathway-trait-and-trajectory)

This vignette contains scripts, explanations, and examples for our
pipeline for genetic mapping of developmental trajectories for complex
traits and diseases. The pipeline is based on integrative analysis of
Genome-Wide Association Studies (GWAS) and single-cell RNA-seq
(scRNA-seq). The analysis performs the following three main tasks:

1.  Identification of connections between developmental trajectories and
    traits.
2.  Elucidate molecular pathways that underlie the link between the
    trajectory and trait
3.  Prioritize genes that carry the link between the pathway, trait and
    trajectory

For the example here we use scRNA-seq dataset of pancreatic islet
development [(Byrnes et
al.)](https://doi.org/10.1038/s41467-018-06176-3), and a GWAS dataset of
type 2 diabetes [(Mahajan et
al.)](https://dx.doi.org/10.1038%2Fs41588-018-0084-1).

The following flowchart summarizes the analysis steps:
![](https://github.com/ElkonLab/scGWAS/blob/master/data/pic/flow.PNG)

# Step 1: Identification of connections between developmental trajectories and traits

## 1a. Converting GWAS variant scores into gene-trait association scores

This is performed using MAGMA gene analysis. We provide example output
[files](https://github.com/ElkonLab/scGWAS/tree/master/data/magma_outputs)
for the type 2 diabetes dataset. Refer to MAGMA’s
[website](https://ctg.cncr.nl/software/magma) and
[manual](https://ctg.cncr.nl/software/MAGMA/doc/manual_v1.07.pdf) for a
detailed explanation. The command used is:

    magma --bfile g1000.eur --pval summary.stat.file  use=rs_id,pval ncol=sample_size  --gene-annot gene.loc.g1000.eur.genes.annot  --out T2D_European.BMIadjusted

Where the *–pval* is the path to the GWAS summary statistics. The
*–gene-annot* is the gene annotation file. See
[manual](https://ctg.cncr.nl/software/MAGMA/doc/manual_v1.07.pdf) for
instructions. We use an annotation window of 10-kbp around the gene. The
*–bfile* gets a path to the file that specifies reference data used to
estimate LD between SNPs. We processed this file using the European
population VCF files from the [1000 Genome
project](https://www.internationalgenome.org/).

## 1b. Calculating cell-trait association scores

Load the following required packages:

``` r
library(monocle)
library(Seurat)
library(data.table)
```

### Preprocessing

We will use extract the count matrix from the monocle object. Fir, Load
the object:

``` r
load('E14_fev_lineage_monocle_ob.Rdata')
```

The count matrix is obtained by:

``` r
exp <- HSMM_seur_var@assayData$exprs
```

To save working memory, it might be a good idea to remove the monocle
object from the environment.

``` r
rm('HSMM_seur_var')
```

We keep genes expressed in at least 10 cells:

``` r
genes.keep <- rowSums(as.matrix(exp) > 0) >= 10
exp <- exp[genes.keep,]
```

Normalize the data using Seurat’s log-normalization:

``` r
exp <- Seurat::NormalizeData(exp)
```

Last, we convert the mouse genes to human orthologues. We provid a
conversion table for ‘one to one’ orthologs
[here](https://github.com/ElkonLab/scGWAS/blob/master/data/rds/ortholog_one2one.RDS).
Other tables can be obtained from
[BioMart](https://www.ensembl.org/biomart/martview/68453bf127464cacc5a5d064c92158e0).

``` r
ortholog_one2one <- readRDS("ortholog_one2one.RDS")
head(ortholog_one2one, 3)
# convert the expression matrix into a data.frame
exp <- cbind.data.frame(data.frame(GENE = rownames(exp)), as.data.frame(exp))
# merge
exp <- merge(ortholog_one2one, exp, by = "GENE")
exp <- exp[,-2]
exp <- exp[!duplicated(exp$GENE),]
```

#### Calculating cell-trait association scores

This part involves running MAGMA’s gene property analysis for each cell
separately. This might take long, so consider using multiple cores.
First, we generate ‘covariate files’ for MAGMA for each cell. Each
output file is a table with the columns: gene, normalized expression in
the cell, and average normalized expression in the dataset. The function
takes as input the data.frame of the normalized expression matrix. The
user specifies the name of the output directory (out.dir) and the number
of cores (cor).

``` r
generate_covs_cells(exp = exp, out.dir = "cov_panc",cors = 20)
```

The next step, it running the gene property analysis, for each cell.
This analysis firs the following the regression model to each cell:

![](https://github.com/eldadshulman/scGWAS/blob/master/data/pic/eq.PNG)

where Z is the vector of the gene’s Z-score converted from the p-values
obtained from MAGMA’s gene analysis for the trait. B is a matrix of
technical confounders, including gene length and SNPs’ LDs, calculated
by MAGMA. C is the vector of normalized expression of the genes in the
cell, and A is a vector of the average normalized expression for the
gene in the dataset. The t-statistic of *β*<sub>*c*</sub> is taken as
the score for the association between the cell and trait.

This function will run the regression for each cell. The input is:

-   magma.path - The path to magma directory
    (“/path/to/MAGMA\_directory/”)

-   covariate.files.dir.path - The path to the covariate files directory
    created by the previews step.

-   raw.file.path - The path to the output (raw) of MAGMA’s gene
    analysis for the trait. We provide the output
    [here](https://github.com/ElkonLab/scGWAS/blob/master/data/magma_outputs/T2D_European.BMIadjusted.genes.raw).

-   out.dir - The name of the output directory.

-   cor - the number of cores to use.

We will get MAGMA gene property output files for each cell. Note that
this step takes even longer.

``` r
gene_prop_cells(cor = 30, raw.file.path = "T2D_European.BMIadjusted.genes.raw", 
                covariate.files.dir.path = "cov_panc", out.dir = "T2D", magma.path = ".")
```

Last, extract the cell scores from MAGMA’s output, namely the
t-statistics.

``` r
cs.pan.t2d <- get_scores(gsa.files.path = "T2D", cors = 20)
```

``` r
head(cs.pan.t2d)
```

    ##              cells         cs
    ## 1 AAACCTGAGAGAGCTC  0.4663375
    ## 2 AAACCTGAGAGTGACC -2.0146530
    ## 3 AAACCTGCATGTAGTC  1.3549249
    ## 4 AAACCTGCATTCACTT -0.9757387
    ## 5 AAACCTGGTAGAAGGA -0.7142417
    ## 6 AAACCTGGTTACTGAC -0.4747647

The first column is cell barcode and the second is the cell score (T2D
risk association score).

``` r
saveRDS(object = cs.pan.t2d, file = "cs.pan.t2d.RDS")
```

## 1c. Trajectory inference.

This step is performed using tools such as Monocle 2, Monocle 3,
destiny. In principle, any tool that gives quantitative maturation
scores for cells, e.g., pseudotime, is suitable. See Monocle’s [version
2](http://cole-trapnell-lab.github.io/monocle-release/docs/), and
[version 3](https://cole-trapnell-lab.github.io/monocle3/) websites for
thorough explanations. For the pancreatic development dataset used here,
the analysis codes were published by the authors and are available for
download
[here](https://figshare.com/articles/software/Scripts_for_Analysis/6783569?backTo=/collections/Lineage_dynamics_of_murine_pancreatic_development_at_single-cell_resolution/4158458).
Also, the output, Monocle object, is available
[here](https://figshare.com/articles/dataset/Monocle_Objects_-_V2_Dataset/6783554?backTo=/collections/Lineage_dynamics_of_murine_pancreatic_development_at_single-cell_resolution/4158458).

## 1d. Examining the association between trait and trajectory

In this step, we will examine the association between trait scores and
pseudotimes of cells.

#### Load the data

Load the following required packages:

``` r
library(monocle)
library(ggplot2)
library(ggpubr)
library(dplyr)
```

From [**step
1b**](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/1b.md),
we got a data.frame with trait (T2D risk) association scores for each
cell. You can download the output
[here](https://github.com/ElkonLab/scGWAS/blob/master/data/rds/cs.pan.t2d.RDS).
We load it:

``` r
cs.pan.t2d <- readRDS("cs.pan.t2d.RDS")
head(cs.pan.t2d, 3)
```

    ##              cells         cs
    ## 1 AAACCTGAGAGAGCTC  0.4663375
    ## 2 AAACCTGAGAGTGACC -2.0146530
    ## 3 AAACCTGCATGTAGTC  1.3549249

From **step 1c**, we have the trajectory, pseudotime, and states. The
monocole object can be downloaded from
[here](https://figshare.com/articles/dataset/Monocle_Objects_-_V2_Dataset/6783554?backTo=/collections/Lineage_dynamics_of_murine_pancreatic_development_at_single-cell_resolution/4158458).
We load it:

``` r
load("E14_fev_lineage_monocle_ob.Rdata")
```

### Preprocessing

We will extract the cell metadata from the monocle object:

``` r
cell.meta <- HSMM_seur_var@phenoData@data
head(cell.meta,3)
```

    ##                  nGene  nUMI orig.ident percent.mito res.0.8 res.1.4
    ## AAACCTGAGAGAGCTC  2607  6448    Fev_E14   0.02946650       4       4
    ## AAACCTGCATTCACTT  3655 10313    Fev_E14   0.02036462       6       5
    ## AAACCTGGTAGAAGGA  2411  5718    Fev_E14   0.02868136       0       0
    ##                  Size_Factor Pseudotime State
    ## AAACCTGAGAGAGCTC   0.8565298   12.32431     3
    ## AAACCTGCATTCACTT   1.3698101   16.46710     3
    ## AAACCTGGTAGAAGGA   0.7595592   16.44184     3

And the coordinates of the cells in the dimensional reduction plot.

``` r
reducedDim <- t(HSMM_seur_var@reducedDimS)
colnames(reducedDim) <- c("dim1", "dim2")
reducedDim <- -reducedDim
head(reducedDim,3)
```

    ##                      dim1       dim2
    ## AAACCTGAGAGAGCTC 1.041375 -0.7681857
    ## AAACCTGCATTCACTT 4.560717 -2.9894874
    ## AAACCTGGTAGAAGGA 4.624069 -2.8805481

To get one data.frame with all the variables, we merge the three
objects:

``` r
pan.t2d <- merge(cs.pan.t2d, cell.meta, by.x = "cells", by.y = "row.names")
pan.t2d <- merge(pan.t2d, reducedDim, by.x = "cells", by.y = "row.names")
head(pan.t2d, 3)
```

    ##              cells         cs nGene nUMI orig.ident percent.mito res.0.8
    ## 1 AAACCTGAGAGAGCTC  0.4663375  2607 6448    Fev_E14   0.02946650       4
    ## 2 AAACCTGAGAGTGACC -2.0146530  2709 6486    Fev_E14   0.02343509       2
    ## 3 AAACCTGCATGTAGTC  1.3549249  2931 7928    Fev_E14   0.02245207       1
    ##   res.1.4 Size_Factor Pseudotime State       dim1       dim2
    ## 1       4   0.8565298   12.32431     3  1.0413749 -0.7681857
    ## 2       1   0.8615776   11.07033     4 -0.5003008 -0.8450683
    ## 3       3   1.0531278   13.33874     5 -0.7624532  5.3630297

We can plot the trajectory, coloring the cells according to pseudotime:

``` r
p <- ggplot(pan.t2d, ggplot2::aes(y = dim2 , x = dim1,
                                  color =  Pseudotime))  +
  geom_point(size=1) + xlab("Component 1") + ylab("Component 2") +
  ggpubr::theme_pubr()  +
  scale_color_continuous(high = "#132B43", low = "#56B1F7") 
print(p)
```

![](https://github.com/ElkonLab/scGWAS/blob/master/data/pic/pseudotime.png)

And, coloring the cells according to T2D risk scores.

``` r
p <- ggplot(pan.t2d, ggplot2::aes(y = dim2 , x = dim1,
                                  color =  cs))  +
  geom_point(size=1) + xlab("Component 1") + ylab("Component 2") +
  ggpubr::theme_pubr()  +
  scale_color_gradient2(low ="#4575B4", mid = "#91BFDB",
                        high = "#D73027", name = "T2D Risk Score")
print(p)
```

![](https://github.com/ElkonLab/scGWAS/blob/master/data/pic/risk.png)

#### Analysis without branches

Now we can perform the regression analysis: cell scores \~ pseudotime.

``` r
mod <- lm(data = pan.t2d, formula = cs ~ Pseudotime)
res <- summary(mod)
print(res$coefficients)
```

    ##                Estimate  Std. Error   t value     Pr(>|t|)
    ## (Intercept) -0.39567031 0.055645311 -7.110578 1.471984e-12
    ## Pseudotime   0.04059159 0.004270521  9.505067 4.243757e-21

To get the p-value, we will us the single-sided test pseudotime
coefficients &gt; 0.

``` r
pv <- pt(coef(res)[, 3], mod$df, lower = FALSE)[2]
print(paste("p-value = ", signif(pv, 2)))
```

    ## [1] "p-value =  2.1e-21"

-   We can visualize the result with a scatter plot:

``` r
p <- ggplot2::ggplot(pan.t2d,ggplot2::aes(x = Pseudotime, y = cs)) + ggplot2::geom_point() + 
  ggpubr::theme_pubr() + ggplot2::ylab("T2D Risk") + ggplot2::xlab("Pseudotime") + 
  ggplot2::coord_fixed(ratio=4) + ggplot2::geom_smooth(method='lm')
print(p)
```

![](https://github.com/ElkonLab/scGWAS/blob/master/data/pic/scatterplot.PNG)

#### Branch Analysis

First, we need to assign branch identity to each cell. For this, we
follow monocle and use the “State” variable.

``` r
p <- ggplot(pan.t2d, ggplot2::aes(y = dim2 , x = dim1,
                                  color =  State))  +
  geom_point(size=1) + xlab("Component 1") + ylab("Component 2") +
  ggpubr::theme_pubr()
print(p)
```

![](https://github.com/ElkonLab/scGWAS/blob/master/data/pic/states.png)

Here, we know that state 1 are the unbranched progenitors. The
beta-branch is state 5, and 2,3,4 are the aplha branch. We divide the
unbranched progenitors into branches by ordering them according to
pseudotime, and assigning odd and even ranked cells to the first and
second branches, respectively. We assigne the first progenitor to both
branches. This is achieved by the following:

``` r
pan.t2d <- branch_assign(pan.t2d, progenitors = 1, branch1 = c(2,3,4), 
                         branch2 = 5, b.names = c("Alpha", "Beta"))

head(pan.t2d)
```

    ##              cells          cs nGene  nUMI orig.ident percent.mito res.0.8
    ## 1 ACACTGAAGAATGTTG  0.55804175  3375 10687    Fev_E14   0.02358226      10
    ## 2 ACACTGAAGAATGTTG  0.55804175  3375 10687    Fev_E14   0.02358226      10
    ## 3 CCATTCGTCGTTGCCT -1.12164422  4215 15785    Fev_E14   0.03376623       8
    ## 4 CGAATGTCAGCAGTTT  0.36772213  3807 13108    Fev_E14   0.03135729      10
    ## 5 CATCAGAAGAGGACGG  0.36563265  3257  9750    Fev_E14   0.02810256      10
    ## 6 CGAACATAGCTTATCG  0.09042509  3726 12361    Fev_E14   0.01860691      10
    ##   res.1.4 Size_Factor Pseudotime State      dim1      dim2 rank  cell.type
    ## 1      12    1.419491 0.00000000     1 -9.781157 -3.676372    0 Progenitor
    ## 2      12    1.419491 0.00000000     1 -9.781157 -3.676372    1 Progenitor
    ## 3       7    2.096824 0.02374406     1 -9.879779 -3.464093    2 Progenitor
    ## 4      12    1.741088 0.12551330     1 -9.723341 -3.528040    3 Progenitor
    ## 5      12    1.295156 0.16231976     1 -9.631442 -3.610558    4 Progenitor
    ## 6      12    1.641992 0.33112614     1 -9.527802 -3.454580    5 Progenitor
    ##   Branch
    ## 1  Alpha
    ## 2   Beta
    ## 3  Alpha
    ## 4   Beta
    ## 5  Alpha
    ## 6   Beta

To examine if the D2T risk association is branch-dependent, we compare
the goodness following models:

``` r
# Null model, assumes that there is no branch-dependency:   
mod0 <- lm(data = pan.t2d, formula = cs ~ Pseudotime)
# model that assumes branch-dependency: 
mod1 <- lm(data = pan.t2d, formula = cs ~ Pseudotime:Branch + Pseudotime)
```

We use Likelihood ratio test to compare the models:

``` r
lkrt <- lmtest::lrtest(mod1, mod0)
pv <- paste0("p-value = ", signif(lkrt$`Pr(>Chisq)`[2],2))
cat(pv)
```

    ## p-value = 5e-21

To examine if there is a positive correlation between T2D
risk-association scores and pseudotime within each branch, we use the
following model.

``` r
mod <- lm(data = pan.t2d, formula = cs ~ Pseudotime:Branch)
res <- summary(mod)
# For one-sided test p-value:
pvs <- pt(coef(res)[, 3], mod$df, lower = FALSE)


pv <- paste0("Alpha: p-value = ", signif(pvs[2],2),
             "\nBeta: p-value = ", signif(pvs[3],2))
cat(pv)
```

    ## Alpha: p-value = 1.2e-19
    ## Beta: p-value = 7.1e-40

-   We can visualize the result with a scatter plot:

``` r
# For reggression lines:
pre <- data.frame(Pseudotime = pan.t2d$Pseudotime)
pre$Branch <- levels(pan.t2d$Branch)[1]

pre2 <- data.frame(Pseudotime = pan.t2d$Pseudotime)
pre2$Branch = levels(pan.t2d$Branch)[2]

pan.t2d$line1 <- predict(mod1, pre)
pan.t2d$line2 <- predict(mod1, pre2)
p <- ggplot2::ggplot(pan.t2d,ggplot2::aes(x = Pseudotime, y = cs, 
                                          color = cell.type)) + 
  ggplot2::geom_point() + ggpubr::theme_pubr() + 
  ggplot2::ylab("T2D Risk") + ggplot2::xlab("Pseudotime") + 
  scale_color_manual(values = c( "grey", "#72ECFA", "#FA8072")) +
  geom_line(aes(Pseudotime ,line1), color = "#08cbe2",
            size =2, show.legend = F) +
  geom_line(aes(Pseudotime,line2), color = "#FA8072", 
            size =2, show.legend = F) +
  ggplot2::coord_fixed(ratio=4)
print(p)
```

![](https://github.com/ElkonLab/scGWAS/blob/master/data/pic/scatterplot_branch.png)
We save the object pan.t2d for the next step.

``` r
saveRDS(object = pan.t2d, 
        file = "meta.pancreas.RDS")
```

# Step 2: Elucidate molecular pathways that underlie the link between the trajectory and trait

## 2a. Finding pathways enriched in the trajectory

In this step we will identify gene sets that are enriched in both the
trajectory of pancreas development.

### Load packeges

Load the following required packages:

``` r
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(speedglm)
```

### Load data

From [**step
1b**](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/1d.md),
we got a data.frame with the expression matrix normalized, filtered, and
with genes converted to human orthologs. You can download the output
[here](https://github.com/ElkonLab/scGWAS/blob/master/data/rds/cs.pan.t2d.RDS).
We load it:

``` r
exp <- readRDS("cm.pancreas.RDS")
head(exp[,1:5], 3)
```

    ##          AAACCTGAGAGAGCTC AAACCTGCATTCACTT AAACCTGGTAGAAGGA AAACGGGTCTTGTACT
    ## C14orf1         0.9364339        1.0782354         1.011187                0
    ## KIAA0141        0.0000000        0.0000000         0.000000                0
    ## KIAA1841        0.0000000        0.6779036         1.011187                0
    ##          AAAGATGGTGCGGTAA
    ## C14orf1         0.7458719
    ## KIAA0141        0.0000000
    ## KIAA1841        0.0000000

From **step 1d**, we go the cell metadata. We load it.

``` r
meta <- readRDS("meta.pancreas.RDS")
head(meta, 3)
```

    ##              cells         cs nGene  nUMI orig.ident percent.mito res.0.8
    ## 1 ACACTGAAGAATGTTG  0.5580418  3375 10687    Fev_E14   0.02358226      10
    ## 2 ACACTGAAGAATGTTG  0.5580418  3375 10687    Fev_E14   0.02358226      10
    ## 3 CCATTCGTCGTTGCCT -1.1216442  4215 15785    Fev_E14   0.03376623       8
    ##   res.1.4 Size_Factor Pseudotime State      dim1      dim2 rank  cell.type
    ## 1      12    1.419491 0.00000000     1 -9.781157 -3.676372    0 Progenitor
    ## 2      12    1.419491 0.00000000     1 -9.781157 -3.676372    1 Progenitor
    ## 3       7    2.096824 0.02374406     1 -9.879779 -3.464093    2 Progenitor
    ##   Branch
    ## 1  Alpha
    ## 2   Beta
    ## 3  Alpha

We load the C5 gene set (GO terms) from msigdbr R packege.

``` r
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g,3)
```

    ## Registered S3 method overwritten by 'cli':
    ##   method     from    
    ##   print.boxx spatstat

    ## # A tibble: 3 x 2
    ##   gs_name                                         entrez_gene
    ##   <chr>                                                 <int>
    ## 1 GO_1_4_ALPHA_OLIGOGLUCAN_PHOSPHORYLASE_ACTIVITY      390637
    ## 2 GO_1_4_ALPHA_OLIGOGLUCAN_PHOSPHORYLASE_ACTIVITY        4507
    ## 3 GO_1_4_ALPHA_OLIGOGLUCAN_PHOSPHORYLASE_ACTIVITY        5834

### Ordering the genes according to pseudotime effect

First, we match the order of the cells in the rows of the metadata, with
that of the expression matrix.

``` r
meta <- meta[match(colnames(exp), meta$cells),]
```

To use GSEA, we first created a ranked list for the genes in the
dataset, according to the incremental change in their expression along
the trajectory. We estimated these incremental changes using GLM to
model the effect of pseudotime on genes’ expression. We use
log-normalized gene expression values and a Gaussian error distribution
GLM. The number of genes detected in each cell was added as a covariate
(nGene).

``` r
p.effect <- df.genes(exp = exp, meta = meta, model = '~Pseudotime:Branch + nGene', brach = 'beta', cors = 10)
```

``` r
head(p.effect, 3)
```

    ##       gene   Estimate Std. Error       t Pr(>|t|)
    ## 1  C14orf1  1.404e-04  2.531e-03  0.0555 9.56e-01
    ## 2 KIAA0141  2.569e-03  1.949e-03  1.3180   0.1876
    ## 3 KIAA1841 -1.836e-02  2.173e-03 -8.4490 4.72e-17

### Running GSEA

To use GSEA from clusterProfiler, we covert gene symbols to ENTREZ IDs.

``` r
  gene.df <- bitr(p.effect$SYMBOL, fromType = "SYMBOL" ,
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  p.effect <- merge(p.effect, gene.df, by = "SYMBOL")
```

We use the t-statistics of the pseudotime:branch coefficients for
ranking the genes.

``` r
  p.effect <- dplyr::arrange(p.effect, dplyr::desc(t))
  geneList <- as.numeric(p.effect$t)
  names(geneList) <- p.effect$gene
```

``` r
 head(geneList, 3)
```

    ##    5126    3375    7276 
    ## 39.7417 39.3662 36.5741

``` r
tail(geneList,3)
```

    ##     2877     6142    50674 
    ## -35.0638 -38.3016 -53.2655

We run GSEA using the ontology gene set (GO, C5, 50 &lt; set size &lt;
500) from MSigDB.

``` r
  edo <-  clusterProfiler::GSEA(geneList = geneList, TERM2GENE =m_t2g, maxGSSize = 500, 
              minGSSize = 50, pvalueCutoff = 1, nPerm = 100000)
  edo <- setReadable(x = edo, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
```

``` r
 head(edo,2)
```

    ##                                                                                                        ID
    ## GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION
    ## GO_CELL_DIVISION                                                                         GO_CELL_DIVISION
    ##                                                                                               Description
    ## GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION
    ## GO_CELL_DIVISION                                                                         GO_CELL_DIVISION
    ##                                                      setSize enrichmentScore
    ## GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION     495      -0.4172673
    ## GO_CELL_DIVISION                                         492      -0.4068760
    ##                                                            NES       pvalue
    ## GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION -1.659186 1.286389e-05
    ## GO_CELL_DIVISION                                     -1.617105 1.288012e-05
    ##                                                          p.adjust      qvalues
    ## GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION 0.0008616726 0.0006815505
    ## GO_CELL_DIVISION                                     0.0008616726 0.0006815505
    ##                                                      rank
    ## GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION 2245
    ## GO_CELL_DIVISION                                     2896
    ##                                                                        leading_edge
    ## GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION tags=26%, list=18%, signal=22%
    ## GO_CELL_DIVISION                                     tags=42%, list=23%, signal=33%
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    core_enrichment
    ## GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION                                                                                                                                                                                                                                                                                                                                                                                                                         MEX3D/NCOR2/UNK/RBMS3/C1QBP/EIF4B/PSMD8/HIST1H3F/RBM24/ZFP36L1/ILF3/ROCK1/UPF3A/HNRNPLL/PSMD12/PSMA7/DKC1/PSMA4/NEURL1/HABP4/PSMD5/POLR2C/NUP85/NEMF/DDX1/ATF4/AGO1/EIF4H/SECISBP2/ELAVL1/ENC1/ZC3H14/PIWIL4/CNOT6/AAAS/PIWIL1/MYCN/TYMS/RBM10/CDC37/EIF3K/NANOS1/LSM14A/PSMC3/MAP2K1/IGF2BP2/SMAD1/PINK1/TAF15/GRB7/FMR1/CNOT2/SYNCRIP/MIF4GD/MTPN/PSMC4/TNRC6B/EIF2S1/ERBB2/GTPBP2/TGFB1/H3F3B/ATG14/HIST2H3A/EIF2B2/PSMB1/APP/CNOT7/ZFP36/PAIP2/PSMD7/APEX1/CNBP/RPL30/PAIP1/DNAJC1/KHDRBS1/EIF4EBP2/YTHDC1/PSMD11/PSMA2/FXR1/MAPK3/DDX39B/E2F1/RNASEL/HIST2H4B/EIF4E3/PSMD10/YTHDF2/SET/TSN/EIF1/EIF2A/NUDT21/HNRNPD/HSPA8/SERBP1/EIF5B/TP53/PSMC2/CNOT6L/MAGOH/PSMD4/TOB1/RPS27A/HNRNPA2B1/PSMC6/NDC1/HNRNPM/ELAVL4/VIM/PPP3CA/EIF3H/DDX5/EEF2/RPL5/UBA52/RPL26/EIF3E/GSPT1/CDK4/RPS9/RPS14/SOX4/RPS4X/BTG2/RPS3
    ## GO_CELL_DIVISION                                     STAMBP/WWTR1/TRIOBP/MAU2/HAUS1/CCNE2/VPS4B/PTCH1/RBL1/E2F8/NEDD1/MASTL/CNTROB/SKA3/CDK7/MAP9/CCNO/SPAG5/BORA/PLK4/CEP63/SIRT2/KATNA1/ERCC6L/NSMCE2/MLLT3/ZW10/GPSM2/FBXO5/RAB35/CCNF/DSN1/RACGAP1/KIF14/ANAPC7/WEE1/USP39/ANAPC4/NCAPH/WASL/CIB1/PLK2/PSRC1/CDC16/NUP62/CKAP2/DIXDC1/NCAPD3/BUB3/CDC6/CHFR/ARPP19/RCC2/ROCK1/SPDL1/FGF9/ECT2/ARHGEF2/MAEA/FGFR2/HAUS8/BRCC3/E2F7/BUB1B/FAM83D/SPICE1/SPAST/ASPM/SMC3/PIK3C3/SNX9/NCAPG2/ZNF830/NCAPD2/FIGN/PKN2/CDK1/SMC4/ZNF207/IL1B/THOC2/SAC3D1/CCND1/BUB1/INCENP/LIG1/CENPW/MAD2L1/PARD3/CDC25C/PLK1/RBBP8/TEAD3/MCMBP/EVI5/CEP55/OIP5/CDCA8/CDCA5/BCL2L1/KIF2A/ZWILCH/CENPJ/HAUS7/SKA1/PTN/NEK2/AKNA/NCAPG/AURKB/CDCA2/LZTS2/ETV5/PPP2R2D/VANGL2/MIS18BP1/HAUS3/TGFB1/SAPCD2/CDCA3/NUF2/TSG101/AURKA/MAD2L2/NEK6/CHMP4B/CCNA2/BIRC5/TOP2A/TERF1/PPP1CB/KIF20B/FGFR1/UBE2C/KIF2C/UBE2S/TIMELESS/PPP1CC/TXNL4A/TACC3/CCND3/CDC42/SPC24/KIF23/TIPIN/SMC1A/NDE1/LBH/CCNB2/HMGA2/KIF20A/REEP3/CCNB1/DLL1/KNSTRN/CDC26/SEPT3/TTC28/RALA/NDC80/KIF11/CDC20/PRC1/TPX2/NUSAP1/SMC2/RHOB/CECR2/CENPE/CALM1/CENPF/CDK2/SH3GLB1/HELLS/RTKN/CAT/CKS2/REEP4/SYCE2/IGF2/STOX1/SDE2/RHOC/PGF/TOP1/PRPF40A/CENPA/SEPT7/SFN/CFL1/SON/PLK5/KIT/NUMBL/CLTA/CDC14B/TXNIP/CCNG2/UBE2I/PLK3/CDK4/STMN1/MDK/RPS3

We can plot, e.g. the “REGULATION\_OF\_INSULIN\_SECRETION” gene-set:

``` r
enrichplot::gseaplot2(edo, geneSetID = "GO_REGULATION_OF_INSULIN_SECRETION", pvalue_table = T)
```

![](https://github.com/ElkonLab/scGWAS/blob/master/data/pic/GSEA.png)

Notice that in our manuscript, the plots are inverted. There we choose
rank genes in an increasing rather than decreasing order. We did it so
that the pseudotime direction is kept from left to right. Therefore,
while in the manuscript we sought for enrichment score (ES) &lt; 0, here
we get the same results from ES &gt; 0. We keep gene-sets with q-value
&lt; 0.05 and NES &gt; 0.

``` r
edo <- edo@result
edo <- edo[which(edo$qvalues < 0.05 & edo$NES > 0),]
saveRDS(object = edo, file = "pan.gene.sets.RDS")
```

## 2b. Examining the trajectory-enriched pathwaysfor trait association

After identifying the gene sets induced during beta-cell development in
[step
2a](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/2a.md), we
now examine their leading-edge set for trait, namely T2D risk,
enrichment. For this, we are using MAGMA gene set analysis.

### Load packeges

Load the following required packages:

``` r
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
```

### Load data

From [step
2a](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/2a.md), we
got the beta=cells enriched gene sets.  
We load it:

``` r
gsea.res <- readRDS("pan.gene.sets.RDS")
head(gsea.res, 3)
```

    ##                                                                  ID
    ## GO_HORMONE_ACTIVITY                             GO_HORMONE_ACTIVITY
    ## GO_RESPIRASOME                                       GO_RESPIRASOME
    ## GO_PROTEIN_N_LINKED_GLYCOSYLATION GO_PROTEIN_N_LINKED_GLYCOSYLATION
    ##                                                         Description setSize
    ## GO_HORMONE_ACTIVITY                             GO_HORMONE_ACTIVITY      60
    ## GO_RESPIRASOME                                       GO_RESPIRASOME      60
    ## GO_PROTEIN_N_LINKED_GLYCOSYLATION GO_PROTEIN_N_LINKED_GLYCOSYLATION      64
    ##                                   enrichmentScore      NES       pvalue
    ## GO_HORMONE_ACTIVITY                     0.6142094 2.025791 2.560098e-05
    ## GO_RESPIRASOME                          0.6237416 2.057230 2.560098e-05
    ## GO_PROTEIN_N_LINKED_GLYCOSYLATION       0.6878882 2.295250 2.585048e-05
    ##                                      p.adjust      qvalues rank
    ## GO_HORMONE_ACTIVITY               0.001038636 0.0008215221  763
    ## GO_RESPIRASOME                    0.001038636 0.0008215221  991
    ## GO_PROTEIN_N_LINKED_GLYCOSYLATION 0.001038636 0.0008215221 1355
    ##                                                     leading_edge
    ## GO_HORMONE_ACTIVITY                tags=22%, list=6%, signal=20%
    ## GO_RESPIRASOME                     tags=42%, list=8%, signal=39%
    ## GO_PROTEIN_N_LINKED_GLYCOSYLATION tags=39%, list=11%, signal=35%
    ##                                                                                                                                                                                                   core_enrichment
    ## GO_HORMONE_ACTIVITY                                                                                                                                  IAPP/TTR/PYY/GIP/NPY/COPA/GAST/PPY/CARTPT/UCN3/GCG/CHGB/PNOC
    ## GO_RESPIRASOME                    NDUFV3/NDUFA13/UQCRQ/NDUFA1/NDUFB4/NDUFB7/NDUFA2/NDUFB10/NDUFA12/NDUFB2/NDUFC1/NDUFS5/NDUFA6/NDUFA3/COX6B1/UQCRB/COX7A1/NDUFB3/NDUFS6/NDUFA8/NDUFS8/NDUFS3/NDUFB6/COX7A2/NDUFS7
    ## GO_PROTEIN_N_LINKED_GLYCOSYLATION                            OSTC/DAD1/DDOST/RPN2/TMEM258/KRTCAP2/RPN1/GFPT1/DPM3/MAGT1/STT3A/LMAN1/DERL3/OST4/UGGT1/TUSC3/DPM2/FUT8/PGM3/ALG12/MOGS/STT3B/ST6GAL1/GORASP1/DPAGT1

### Input files for MAGMA

-   **–gene-results**: This is the output from MAGMA’s gene analysis.
    The T2D is available
    [here](https://github.com/ElkonLab/scGWAS/blob/master/data/magma_outputs/T2D_European.BMIadjusted.genes.raw).

-   **–set-annot**: This is the gene-set file. We prepare a gene-set
    file from the leading-edge subsets of the beta-cell enriched gene
    sets.

``` r
# Loop over gene-sets
  for(j in 1:nrow(gsea.res)){
    # Extract leading-edge genes, from the j's gene set.
    g.set <- data.frame(gene = strsplit(gsea.res$core_enrichment[j], "/")[[1]])
    # Assign gene-set name column
    g.set$set <- gsea.res$ID[j]
    g.set <- g.set[,c(2,1)]
    # Combine with other gene sets.
    if(!exists("set.annot")) set.annot <- g.set
    if(exists("set.annot")) set.annot <- rbind.data.frame(set.annot,g.set)
    
  } 
```

``` r
head(set.annot, 3)
```

    ##                    V1   V2
    ## 1 GO_HORMONE_ACTIVITY IAPP
    ## 2 GO_HORMONE_ACTIVITY  TTR
    ## 3 GO_HORMONE_ACTIVITY  PYY

``` r
  # write file:
  fwrite(x = set.annot, sep = "\t", col.names = F, row.names = F, quote = F, 
         file = "panc.leading-edge.txt")
```

### Using MAGMA’s gene-set analysis

We use MAGMA, providing the paths to the above inputs (–gene-results,
–set-annot).

    magma --settings gene-info --gene-results T2D_European.BMIadjusted.genes.raw --set-annot panc.leading-edge.txt col=2,1 --out t2d.panc

The output file, ‘t2d.panc.gsa.out’, gives us the association between
T2D and the leading-edge gene-sets.

``` r
gsa <- read.table("t2d.panc.gsa.out", header = T)

head(gsa, 3)
```

                             VARIABLE TYPE NGENES      BETA   BETA_STD      SE
    1 GO_ACTIVE_TRANSMEMBRANE_TRAN...  SET     33  0.144420  0.0064990 0.18230
    2 GO_ADENYLATE_CYCLASE_ACTIVAT...  SET      8  0.099064  0.0021967 0.35044
    3        GO_ATP_METABOLIC_PROCESS  SET     48 -0.056490 -0.0030646 0.13070
            P
    1 0.21413
    2 0.38871
    3 0.66720
                                                                         FULL_NAME
    1                                 GO_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY
    2 GO_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY
    3                                                     GO_ATP_METABOLIC_PROCESS

# Step 3: Prioritize genes that carry the link between the pathway, trait and trajectory

In this step we will identify gene sets that are enriched in both the
trajectory of pancreas development.

### Load packeges

Load the following required packages:

``` r
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
```

### Load data

From **step 1a**, we got an output file (“gene.out”) from MAGMA’s gene
analysis, giving gene-trait (T2D) association p-values. It can be
downloaded from
[here](https://github.com/ElkonLab/scGWAS/blob/master/data/magma_outputs/T2D_European.BMIadjusted.genes.out).
We load it:

``` r
gene.p <- read.table("T2D_European.BMIadjusted.genes.out", header = T)
head(gene.p, 3)
```

    ##     GENE CHR  START   STOP NSNPS NPARAM       N   ZSTAT       P
    ## 1 SAMD11   1 859993 879961    13      8 1740000 0.20813 0.41756
    ## 2  NOC2L   1 879583 894679    19     13 1740000 0.97740 0.16419
    ## 3 KLHL17   1 895967 901099    19     13 1740000 0.40072 0.34431

From [**step
2b**](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/2b.md),
we got the output from MAGMA’s gene-set analysis. This file has T2D risk
association p-values for the leading-edge gene sets enriched in the
pancreas trajectory. It can be downloaded from
[here](https://github.com/ElkonLab/scGWAS/blob/master/data/magma_outputs/t2d.panc.gsa.out).
We load it:

``` r
gsa <- read.table("t2d.panc.gsa.out", header = T)

head(gsa, 3)
```

                             VARIABLE TYPE NGENES      BETA   BETA_STD      SE
    1 GO_ACTIVE_TRANSMEMBRANE_TRAN...  SET     33  0.144420  0.0064990 0.18230
    2 GO_ADENYLATE_CYCLASE_ACTIVAT...  SET      8  0.099064  0.0021967 0.35044
    3        GO_ATP_METABOLIC_PROCESS  SET     48 -0.056490 -0.0030646 0.13070
            P
    1 0.21413
    2 0.38871
    3 0.66720
                                                                         FULL_NAME
    1                                 GO_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY
    2 GO_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY
    3                                                     GO_ATP_METABOLIC_PROCESS

From the same step, we also got a file with the trajectory significant
gene-sets and their leading-edge genes.

``` r
set.annot <- read.table("panc.leading-edge.txt", header = F)
colnames(set.annot) <- c("Gene-set", "GENE")
head(set.annot, 3)
```

                 Gene-set GENE
    1 GO_HORMONE_ACTIVITY IAPP
    2 GO_HORMONE_ACTIVITY  TTR
    3 GO_HORMONE_ACTIVITY  PYY

### Processing the data

First, we keep gene-sets significantly associated with T2D:

``` r
gsa <- gsa[which(gsa$P < 0.05),]

sig.gene.sets <- gsa$FULL_NAME
```

Now, we keep the genes from the leading-edge of the significant
(trajectory and trait) gene sets:

``` r
risk.genes <- set.annot[which(set.annot$`Gene-set` %in% sig.gene.sets),]
```

Many genes are members of more than one set and therefore appear
multiple times in this table. We use dplyr to collapse the rows of such
genes as follows:

``` r
risk.genes <- dplyr::group_by(risk.genes, GENE) %>% dplyr::summarise(`Gene-sets` = paste(`Gene-set`, collapse = ", "))
```

    `summarise()` ungrouping output (override with `.groups` argument)

``` r
head(risk.genes,3)
```

    # A tibble: 3 x 2
      GENE  `Gene-sets`                                                             
      <chr> <chr>                                                                   
    1 AACS  GO_HORMONE_TRANSPORT, GO_REGULATION_OF_PEPTIDE_SECRETION, GO_REGULATION~
    2 ABAT  GO_HORMONE_TRANSPORT, GO_REGULATION_OF_PEPTIDE_SECRETION, GO_REGULATION~
    3 ABCA2 GO_NEGATIVE_REGULATION_OF_TRANSPORT, GO_VACUOLAR_MEMBRANE               

Now we have a table for genes from the leading-edge gene sets associated
with both trajectory and trait. We merge this table with the gene-trait
p-values from MAGMA to prioritize the genes, keeping only genes with
p-value &lt; 0.05.

``` r
risk.genes <- merge(risk.genes, gene.p[,c(1,8,9)], by = "GENE") 
risk.genes <- dplyr::filter(risk.genes, P < 0.05 ) %>% dplyr::arrange(P)
head(risk.genes,3)
```

         GENE
    1 SLC30A8
    2   KCNQ1
    3    WFS1
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Gene-sets
    1                                                                                                                                                                                                                   GO_TRANSPORT_VESICLE_MEMBRANE, GO_HORMONE_TRANSPORT, GO_REGULATION_OF_PEPTIDE_SECRETION, GO_TRANSPORT_VESICLE, GO_REGULATION_OF_HORMONE_LEVELS, GO_PEPTIDE_SECRETION, GO_SIGNAL_RELEASE, GO_SECRETORY_GRANULE_MEMBRANE, HP_ADULT_ONSET, GO_METAL_ION_HOMEOSTASIS, HP_ABNORMAL_GLUCOSE_HOMEOSTASIS, GO_CELLULAR_ION_HOMEOSTASIS, GO_DIVALENT_INORGANIC_CATION_HOMEOSTASIS
    2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  GO_REGULATION_OF_CARDIAC_MUSCLE_CONTRACTION, GO_REGULATION_OF_STRIATED_MUSCLE_CONTRACTION
    3 GO_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE, GO_CELLULAR_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN, GO_TRANSPORT_VESICLE_MEMBRANE, GO_PROTEIN_MATURATION, GO_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS, GO_INTRINSIC_COMPONENT_OF_ORGANELLE_MEMBRANE, GO_TRANSPORT_VESICLE, GO_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN, GO_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE, HP_ADULT_ONSET, GO_INTRINSIC_COMPONENT_OF_ENDOPLASMIC_RETICULUM_MEMBRANE, GO_METAL_ION_HOMEOSTASIS, HP_ABNORMAL_GLUCOSE_HOMEOSTASIS, GO_CELLULAR_ION_HOMEOSTASIS, GO_DIVALENT_INORGANIC_CATION_HOMEOSTASIS
        ZSTAT          P
    1 11.5250 4.9257e-31
    2 10.6780 6.4236e-27
    3  8.4311 1.7116e-17

We color the cells in the trajectory plot according to the expression of
the genes as follows. We load the cell metadata from [**step
1d**](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/1d.md).

``` r
meta <- readRDS("meta.pancreas.RDS")
```

And the normalized, human ortholog converted matrix from [**step
1b**](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/1b.md).

``` r
exp <- readRDS("cm.pancreas.RDS")
```

``` r
# genes to plot
genes.plot <- risk.genes$GENE[1:3]
print(genes.plot)
```

    [1] "SLC30A8" "KCNQ1"   "WFS1"   

``` r
#Subset the expression matrix.
exp <- exp[which(rownames(exp) %in% genes.plot),]

# scale data using Seurat
exp.scaled <- Seurat::ScaleData(exp)
exp.scaled <- t(exp.scaled)
# Floor using Seurat
exp.scaled <- apply(exp.scaled, 2, Seurat::MinMax, 
                    min = -1.5, max = 1.5)
exp.scaled <- cbind.data.frame(data.frame(cells = rownames(exp.scaled)), 
                               as.data.frame(exp.scaled))

# Merge gene expression and dimensional reduction coordinates.
exp.scaled <- merge(exp.scaled, meta, by = "cells")

for(i in genes.plot){
p <- ggplot(exp.scaled, ggplot2::aes(y = dim2 , x = dim1,
                                  color =  get(i)))  +
  geom_point(size=1) + xlab("Component 1") + ylab("Component 2") +
  ggpubr::theme_pubr()  +
  scale_color_gradient2(low ="#FF00FF", mid = "#000000", midpoint = 0,
                        high ="#FFFF00", name = i)
print(p)
}
```

![](https://github.com/ElkonLab/scGWAS/blob/master/data/pic/slc.png)
![](https://github.com/ElkonLab/scGWAS/blob/master/data/pic/kcnq1.png)
![](https://github.com/ElkonLab/scGWAS/blob/master/data/pic/wfs1.png)

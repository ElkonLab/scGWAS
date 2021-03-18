Genetic mapping of developmental trajectories for complex traits and
diseases
================

This repository contains scripts, explanations, and examples for our
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

## Prerequisites

  - [MAGMA](https://ctg.cncr.nl/software/magma)
  - Download our R functions
    [here](https://github.com/ElkonLab/scGWAS/blob/master/R/functions_scGWAS.R)
  - Make sure the following packages are Installed: Monocole (v2), plyr,
    ggplot2, data.table, parallel, Seurat, clusterProfiler (a
    Bioconductor package).
  - Some input files are also required. For details and sample files,
    see individual vignettes bellow.

## Examples for usage

We provide vignettes for each of these steps, using scRNA-seq dataset of
pancreatic islet development [(Byrnes et
al.)](https://doi.org/10.1038/s41467-018-06176-3), and a GWAS dataset of
type 2 diabetes [(Mahajan et
al.)](https://dx.doi.org/10.1038%2Fs41588-018-0084-1). The following are
links to each vignette:

  - **Task 1**: Identification of connections between developmental
    trajectories and traits vignette:
    [link](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/first_step.md)
  - **Task 2**: Elucidate molecular pathways that underlie the link
    between the trajectory and trait vignette:
    [link](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/first_step.md)
  - **Task 3**: Prioritize genes that carry the link between the
    pathway, trait and trajectory vignette:
    [link](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/first_step.md)

## Authors

  - Eldad David Shulman

  - Prof. Ran Elkon

## License

This project is licensed under the BSD 3 License - see the
[LICENSE.md](https://github.com/ElkonLab/scGWAS/blob/master/LICENSE.md)
file for details.

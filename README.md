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

-   [MAGMA](https://ctg.cncr.nl/software/magma)
-   Download our R functions
    [here](https://github.com/ElkonLab/scGWAS/blob/master/R/functions_scGWAS.R)
-   Make sure the following packages are Installed: Monocole (v2), plyr,
    ggplot2, data.table, parallel, Seurat, clusterProfiler (a
    Bioconductor package).
-   Some input files are also required. For details and sample files,
    see individual vignettes bellow.

## Examples for usage

We provide vignettes for the analysis using scRNA-seq dataset of
pancreatic islet development [(Byrnes et
al.)](https://doi.org/10.1038/s41467-018-06176-3), and a GWAS dataset of
type 2 diabetes [(Mahajan et
al.)](https://dx.doi.org/10.1038%2Fs41588-018-0084-1). The following
details the analysis steps and provides links to the vignette: \* **Step
1: Identification of connections between developmental trajectories and
traits.** This involves the following:

-   **1a. Converting GWAS variant scores into gene-trait association
    scores.** This is performed using MAGMA gene analysis. We provide
    example output
    [files](https://github.com/ElkonLab/scGWAS/tree/master/data/magma_outputs)
    for the type 2 diabetes dataset. Refer to MAGMA’s
    [website](https://ctg.cncr.nl/software/magma) and
    [manual](https://ctg.cncr.nl/software/MAGMA/doc/manual_v1.07.pdf)
    for a detailed explanation.

-   **1b. Calculating cell-trait association scores.** See our
    [vignette](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/1b.md)

-   **1c. Trajectory inference.** This step is performed using tools
    such as Monocle 2, Monocle 3, destiny. In principle, any tool that
    gives quantitative maturation scores for cells, e.g., pseudotime, is
    suitable. See Monocle’s [version
    2](http://cole-trapnell-lab.github.io/monocle-release/docs/), and
    [version 3](https://cole-trapnell-lab.github.io/monocle3/) websites
    for thorough explanations. For the pancreatic development dataset
    used here, the analysis codes were published by the authors and are
    available for download
    [here](https://figshare.com/articles/software/Scripts_for_Analysis/6783569?backTo=/collections/Lineage_dynamics_of_murine_pancreatic_development_at_single-cell_resolution/4158458).
    Also, the output, Monocle object, is available
    [here](https://figshare.com/articles/dataset/Monocle_Objects_-_V2_Dataset/6783554?backTo=/collections/Lineage_dynamics_of_murine_pancreatic_development_at_single-cell_resolution/4158458).

-   **1d. Examining the association between trait and trajectory.** See
    our
    [vignette](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/1d.md).

-   **Step 2: Elucidate molecular pathways that underlie the link
    between the trajectory and trait.** This includes first

    -   **2a.** finding pathways enriched in the trajectory. This is
        covered
        [here](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/2a.md).
    -   **2b.** examining if they are enriched for the trait. This is
        covered
        [here](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/2a.md).

-   **Step 3: Prioritize genes that carry the link between the pathway,
    trait and trajectory.** See our [vignette
    here](https://github.com/ElkonLab/scGWAS/blob/master/vignettes/2a.md).

The following flowchart summarizes the analysis steps:
![](https://github.com/ElkonLab/scGWAS/blob/master/data/pic/flow.PNG)

## Authors

-   Eldad David Shulman

-   Prof. Ran Elkon

## License

This project is licensed under the BSD 3 License - see the
[LICENSE.md](https://github.com/ElkonLab/scGWAS/blob/master/LICENSE.md)
file for details.

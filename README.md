# Experimental Design and Differential Inference for Comparative Single-cell RNA-sequencing Studies

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)

## Overview

Nowadays, single-cell RNA-sequencing (scRNA-seq) experiments that profile cells from multiple treatment or biological conditions have become prevalent. However, despite the urgent demands, guidelines on designing a valid comparative scRNA-seq study and statistical methods that can rigorously quantify the differences between scRNA-seq datasets are still lacking.

Although solid statistical methods for calling differentially expressed (DE) genes have been the cornerstone of genomic research in the bulk data era, with Limma, DESeq2 and edgeR each serving tens of thousands of studies, these classic methods are not directly applicable to scRNA-seq data. The main challenge is that although we know the condition label of each cell---whether the cell was collected from a healthy individual or a patient, the cell-type label of each cell is unknown. Therefore, for differential inference of scRNA-seq data, we have to simultaneously (a) infer the cell-type label of each individual cell, (b) identify the cell-type-specific DE genes---genes that are DE between conditions for each cell type and (c) quantify the changes in cell-type proportions between conditions, in another word, infer differential abundance (DA).

Here, we develop DIFferential inference for scRNA-seq data (DIFseq), a statistical method that is able to rigorously quantify the condition effects on both cell-type-specific gene expression levels and cellular compositions for scRNA-seq data by accounting for all the uncertainties in the analysis with an interpretable hierarchical model. DIFseq simultaneously adjusts batch effects, characterizes the cell-type-specific size factors, the count data nature and the dropout events of scRNA-seq data, clusters cells into cell types and quantifies the condition effects.

## Repo Contents

- [R](./R): `R` code.
- [data](./data): the example data for the demo.
- [inst/doc](./inst/doc): compiled user's guide and illustration of the applications of `BUSseq` package to the demo dataset.
- [man](./man): help manual.
- [src](./src): `C++` source code.
- [tests](./tests): sample code for the demo dataset.
- [vignettes](./vignettes): source code for the user's guide.

## System Requirements

#### Hardware Requirements

The `BUSseq` package works on a standard personal computer (PC). The runtimes reported below were generated on an Ubuntu 18.04 operating system by a PC desktop with 8 GB RAM and 8 cores of 2.6 GHz.

### Software Requirements

The package supports *Windows* operating systems. It has been tested on Windows 10 Enterprise. Before installing the `BUSseq` package, users should have installed `R` with version 4.1 or higher. The users should also install [`Rtools`](https://cran.r-project.org/bin/windows/Rtools/).

## Installation Guide

From an `R` session, please type:

```
install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BiocGenerics")
BiocManager::install("edgeR")
BiocManager::install("S4Vectors")
BiocManager::install("SingleCellExperiment")
BiocManager::install("SummarizedExperiment")

library(devtools)

# Install DIFseq without building the vignette
install_github("songfd2018/DIFseq-Rpackage") # install DIFseq

# Install DIFseq with building the vignette
# install_github("songfd2018/DIFseq-Rpackage", build_vignettes = TRUE)
```

It takes about two minutes to install the package without building the vignette, and the vignette in the html format can be viewed [here](https://htmlpreview.github.io/?https://github.com/songfd2018/DIFseq_tutorial/blob/main/DIFseq%20User%20Guide.html). Alternatively, the package can also be installed with building the vignette. However, it may take about half an hour to install and build the package.

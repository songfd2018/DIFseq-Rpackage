# Experimental Design and Differential Inference for Comparative Single-cell RNA-sequencing Studies

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Instructions for use](#instructions-for-use)
- [License](./LICENSE)
- [Citation](#citation)

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

## Installation Guide

From an `R` session, type:

```
require(devtools)
install_github("songfd2018/DIFseq-Rpackage") # install BUSseq
```

## Demo

Please check the user's guide for detailed instructions on how to use the package by running the following code in the `R` session:

```
vignette("DIFseq_user_guide",package="DIFseq")  # view the vignettes
```

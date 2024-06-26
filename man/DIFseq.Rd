% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Object.R
\docType{class}
\name{DIFseq}
\alias{DIFseq}
\alias{CreateDIFseqObject}
\alias{DIFseq-class}
\title{The DIFseq constructor}
\usage{
CreateDIFseqObject(
  ...,
  batch = NULL,
  condition = NULL,
  sample = NULL,
  n.celltype = 5,
  hyperparam = list(),
  iter = list()
)
}
\arguments{
\item{...}{Arguments passed to the DIFseq constructor to fill the slots of the
base class. This should be either a \code{\linkS4class{SingleCellExperiment}}
or a feature-by-cell matrix}

\item{batch}{A vector with the same length as the number of cells indicating
the batch label of each cell. By default, all the cells are assigned to
the same batch.}

\item{condition}{A vector with the same length as the number of cells indicating
the condition label of each cell. By default, all the cells are assigned to
the same condition.}

\item{sample}{A vector with the same length as the number of cells indicating
the sample label of each cell. By default, each batch-condition pair is regarded
as a sample.}

\item{n.celltype}{the number of cell type}

\item{hyperparam}{A list of hyper-parameter values}

\item{iter}{A list of parameter to control the EM algorithm}
}
\value{
A DIFseq object
}
\description{
The DIFseq class is designed to store the parameter estimation and
differential inference results by DIFseq mode. It inherits
the \linkS4class{SingleCellExperiment} class and can be used in the same manner.
In addition, the class stores the parameter estimation
}
\details{
A DIFseq object can be constructed by inputting a gene-by-cell gene
expression matrix. In this case it simply constructs a SingleCellExperiment
and fills the relevant slots.
}
\examples{

library(SingleCellExperiment)
ux <- matrix(rpois(12000, 5), ncol=200)
vx <- log2(ux + 1)
pca <- prcomp(t(vx))

sce <- SingleCellExperiment(assays=list(counts=ux, logcounts=vx),
                            reducedDims=SimpleList(PCA=pca$x))

DIFseq_obj <- CreateDIFseqObject(sce)
DIFseq_obj

}
\author{
Fangda Song, Kevin Y. Yip and Yingying Wei
}

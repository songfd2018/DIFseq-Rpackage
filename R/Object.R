#' The DIFseq constructor
#'
#' The DIFseq class is designed to store the parameter estimation and
#' differential inference results by DIFseq mode. It inherits
#' the \linkS4class{SingleCellExperiment} class and can be used in the same manner.
#' In addition, the class stores the parameter estimation
#'
#' @param ... Arguments passed to the DIFseq constructor to fill the slots of the
#' base class. This should be either a \code{\linkS4class{SingleCellExperiment}}
#' or a feature-by-cell matrix
#' @param Batch A vector with the same length as the number of cells indicating
#' the batch label of each cell. By default, all the cells are assigned to
#' the same batch.
#' @param Condition A vector with the same length as the number of cells indicating
#' the condition label of each cell. By default, all the cells are assigned to
#' the same condition.
#' @param n.celltype the number of cell type
#' @param hyperparam A list of hyper-parameter values
#' @param iter A list of parameter to control the EM algorithm
#
# @details
# In this class the underlying structure is the gene/feature-by-cell expression
# data. The additional slots provide a link between these single cells and
# the neighbourhood representation. This can be further extended by the use
# of an abstracted graph for visualisation that preserves the structure of the
# single-cell KNN-graph
#
#' @details
#' A DIFseq object can be constructed by inputting a gene-by-cell gene
#' expression matrix. In this case it simply constructs a SingleCellExperiment
#' and fills the relevant slots.
#'
#' @return A DIFseq object
#'
#' @author Fangda Song, Kevin Y. Yip and Yingying Wei
#'
#' @examples
#'
#' library(SingleCellExperiment)
#' ux <- matrix(rpois(12000, 5), ncol=200)
#' vx <- log2(ux + 1)
#' pca <- prcomp(t(vx))
#'
#' sce <- SingleCellExperiment(assays=list(counts=ux, logcounts=vx),
#'                             reducedDims=SimpleList(PCA=pca$x))
#'
#' DIFseq_obj <- CreateDIFseqObject(sce)
#' DIFseq_obj
#'
#' @docType class
#' @name DIFseq
#'
#'
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Matrix Matrix
CreateDIFseqObject <- function(...,
                               batch = NULL,
                               condition = NULL,
                               sample = NULL,
                               n.celltype = 5,
                               hyperparam = list(),
                               iter = list()){

  sce <- SingleCellExperiment(...)

  .sce_to_DIFseq(sce,
                 batch = batch,
                 condition = condition,
                 sample = sample,
                 n.celltype = n.celltype,
                 hyperparam = hyperparam,
                 iter = iter)

}

#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
#' @importFrom BiocGenerics ncol
.sce_to_DIFseq <- function(sce,
                           batch = NULL,
                           condition = NULL,
                           sample = NULL,
                           n.celltype = 5,
                           hyperparam = list(),
                           iter = list()){

  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }

  n.cell <- ncol(sce)

  # Set default value for batch and condition informaiton
  if(is.null(batch)){
    batch = factor(rep(1,n.cell))
  }

  if(is.null(condition)){
    condition = factor(rep(1, n.cell))
  }

  if(is.null(sample)){
    sample = factor(paste0("B",batch,"_C",condition))
  }

  sce$Batch <- batch
  sce$Condition <- condition
  sce$Pair <- factor(paste0("B",batch,"_C",condition))
  sce$Sample <- sample

  # Set default value for hyper-parameter values
  if(is.null(hyperparam$hyper_pi)){
    hyperparam$hyper_pi = 2
  }

  if(is.null(hyperparam$hyper_gamma0)){
    hyperparam$hyper_gamma0 = 3
  }

  if(is.null(hyperparam$hyper_gamma1)){
    hyperparam$hyper_gamma1 = c(0.001, 0.01)
  }

  if(is.null(hyperparam$hyper_alpha)){
    hyperparam$hyper_alpha = 5
  }

  if(is.null(hyperparam$tau0sq)){
    hyperparam$tau0sq = 0.01
  }

  if(is.null(hyperparam$hyper_p)){
    hyperparam$hyper_p = c(1,3)
  }

  if(is.null(hyperparam$hyper_tau1sq)){
    hyperparam$hyper_tau1sq = c(2,500)
  }

  if(is.null(hyperparam$hyper_nu)){
    hyperparam$hyper_nu = 5
  }

  if(is.null(hyperparam$hyper_delta)){
    hyperparam$hyper_delta = 5
  }

  if(is.null(hyperparam$hyper_phi)){
    hyperparam$hyper_phi = c(1,0.1)
  }

  if(is.null(hyperparam$hyper_pbeta)){
    hyperparam$hyper_pbeta = c(1,4)
  }

  if(is.null(hyperparam$hyper_peta)){
    hyperparam$hyper_peta = c(1,9)
  }

  # Set default value for iteration control
  if(is.null(iter$learn_rate)){
    iter$learn_rate = c(1,1) # learning rate and decay rate
  }

  if(is.null(iter$mini_batch)){
    iter$mini_batch = min(1000, n.cell)
  }

  if(is.null(iter$early_stop)){
    iter$early_stop = 10
  }

  if(is.null(iter$iter_steps_limit)){
    iter$iter_steps_limit = c(1000, 100) # Maximum and minimum iteration number
  }

  if(is.null(iter$check_per_iter)){
    iter$check_per_iter = 5
  }

  if(is.null(iter$MC_rep)){
    iter$MC_rep = 10
  }

  if(is.null(iter$Valid_size)){
    iter$Valid_size = 2
  }

  if(is.null(iter$n.cores)){
    iter$n.cores = 8
  }

  if(is.null(iter$seed)){
    iter$seed = 1234
  }

  out <- new("DIFseq", sce,
             n.celltype = n.celltype,
             hyperparam = hyperparam,
             iter = iter,
             estimation = list())

  return(out)
}

# The DIFseq class
#
# @slot n.celltype Number of cell types
# @slot hyperparam A list of hyperparameters
# @slot iter A list of parameters to control the EM algorithm
# @slot estimation A list of estimated parameter values
#' @aliases DIFseq
#' @rdname DIFseq
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors SimpleList
setClass("DIFseq",
         contains = "SingleCellExperiment",
         slots=c(
           n.celltype = "numeric",
           hyperparam = "list",
           iter = "list",
           estimation = "list"
         ),
         prototype = list(
           n.celltype = 5,
           hyperparam = list(),
           iter = list(),
           estimation = list()
         )
)

# All Generics
#' @export
setGeneric("NumCond", function(x) standardGeneric("NumCond"))

#' @export
setGeneric("NumBatch", function(x) standardGeneric("NumBatch"))

#' @export
setGeneric("NumPair", function(x) standardGeneric("NumPair"))

#' @export
setGeneric("NumSample", function(x) standardGeneric("NumSample"))

#' @export
setGeneric("BatchInd", function(x) standardGeneric("BatchInd"))

#' @export
setGeneric("CondInd", function(x) standardGeneric("CondInd"))

#' @export
setGeneric("PairInd", function(x) standardGeneric("PairInd"))

#' @export
setGeneric("SampleInd", function(x) standardGeneric("SampleInd"))

#' @export
setGeneric("MiniBatch", function(x) standardGeneric("MiniBatch"))

#' @export
setGeneric("MiniBatch<-", function(x, value) standardGeneric("MiniBatch<-"))

#' @export
setGeneric("EarlyStop", function(x) standardGeneric("EarlyStop"))

#' @export
setGeneric("EarlyStop<-", function(x, value) standardGeneric("EarlyStop<-"))

#' @export
setGeneric("CheckLike", function(x) standardGeneric("CheckLike"))

#' @export
setGeneric("CheckLike<-", function(x, value) standardGeneric("CheckLike<-"))

#' @export
setGeneric("MCRep", function(x) standardGeneric("MCRep"))

#' @export
setGeneric("MCRep<-", function(x, value) standardGeneric("MCRep<-"))

#' @export
setGeneric("Ncores", function(x) standardGeneric("Ncores"))

#' @export
setGeneric("Ncores<-", function(x, value) standardGeneric("Ncores<-"))

#' @export
setGeneric("Seed", function(x) standardGeneric("Seed"))

#' @export
setGeneric("Seed<-", function(x, value) standardGeneric("Seed<-"))

#' @export
setGeneric("BIC_DIFseq", function(x) standardGeneric("BIC_DIFseq"))

#' @export
setGeneric("Prop", function(x) standardGeneric("Prop"))

#' @export
setGeneric("CellTypes", function(x) standardGeneric("CellTypes"))

#' @export
setGeneric("Baseline", function(x) standardGeneric("Baseline"))

#' @export
setGeneric("TypeEffects", function(x) standardGeneric("TypeEffects"))

#' @export
setGeneric("CondEffects", function(x) standardGeneric("CondEffects"))

#' @export
setGeneric("BatchEffects", function(x) standardGeneric("BatchEffects"))

#' @export
setGeneric("OverDisp", function(x) standardGeneric("OverDisp"))

#' @export
setGeneric("CellEffects", function(x) standardGeneric("CellEffects"))

#' @export
setGeneric("DropoutCoef", function(x) standardGeneric("DropoutCoef"))

#' @export
setGeneric("ImputedCounts", function(x) standardGeneric("ImputedCounts"))

#' @export
setGeneric("CorrectedCounts", function(x) standardGeneric("CorrectedCounts"))

#' @export
setGeneric("Intrinsic", function(x) standardGeneric("Intrinsic"))

#' @export
setGeneric("DiffExp", function(x) standardGeneric("DiffExp"))

#' @export
setGeneric("DiffAbund", function(x) standardGeneric("DiffAbund"))

# # @importFrom S4Vectors SimpleList
# # @importFrom Matrix Matrix
# # @import SingleCellExperiment
# .fromSCE <- function(sce){
#   # make the distance and adjacency matrices the correct size
#   out <- new("DIFseq", sce,
#              graph=list(),
#              nhoods=Matrix(0L, sparse=TRUE),
#              nhoodDistances=NULL,
#              nhoodCounts=Matrix(0L, sparse=TRUE),
#              nhoodIndex=list(),
#              nhoodExpression=Matrix(0L, sparse=TRUE),
#              .k=NULL)
#
#   reducedDims(out) <- reducedDims(sce)
#   altExps(out) <- list()
#
#   return(out)
# }
#
# # @importFrom Matrix Matrix
# # @importFrom S4Vectors DataFrame SimpleList
# # @importFrom SingleCellExperiment colData rowData altExps reducedDims colPairs rowPairs
# .fromMatrix <- function(mat){
#   # return an empty Milo object
#   out <- new("DIFseq",
#              SingleCellExperiment(mat),
#              graph=list(),
#              nhoods=Matrix(0L, sparse=TRUE),
#              nhoodDistances=NULL,
#              nhoodCounts=Matrix(0L, sparse=TRUE),
#              nhoodIndex=list(),
#              nhoodExpression=Matrix(0L, sparse=TRUE),
#              .k=NULL)
#
#   reducedDims(out) <- reducedDims(sce)
#   altExps(out) <- list()
#
#   if (objectVersion(out) >= "1.11.3"){
#     colPairs(out) <- SimpleList()
#     rowPairs(out) <- SimpleList()
#   }
#
#   return(out)
# }
#
# # @importFrom Matrix Matrix
# # @importFrom S4Vectors DataFrame SimpleList
# # @importFrom SingleCellExperiment colData rowData altExps reducedDims colPairs rowPairs
# .emptyMilo <- function(...){
#   # return an empty Milo object
#   out <- new("Milo",
#              graph=list(),
#              nhoods=Matrix(0L, sparse=TRUE),
#              nhoodDistances=NULL,
#              nhoodCounts=Matrix(0L, sparse=TRUE),
#              nhoodIndex=list(),
#              nhoodExpression=Matrix(0L, sparse=TRUE),
#              .k=NULL,
#              int_elementMetadata=DataFrame(),
#              int_colData=DataFrame())
#
#   altExps(out) <- SimpleList()
#   reducedDims(out) <- SimpleList()
#
#   if (objectVersion(out) >= "1.11.3"){
#     colPairs(out) <- SimpleList()
#     rowPairs(out) <- SimpleList()
#   }
#
#   out
# }



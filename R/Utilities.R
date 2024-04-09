#' @title
#' Methods for extracting statistical inference results from a DIFseq object
#'
#' @description
#' Get or set methods for Milo object slots. Generally speaking these methods
#' are used internally, but they allow the user to assign their own externally computed
#' values - should be used \emph{with caution}.
#'
#' @section Dimension information:
#' In the following descriptions \code{x} is always a \linkS4class{DIFseq} object.
#' \describe{
#' \item{\code{NumCond(x)}:}{Return the number of conditions \code{T}.}
#' \item{\code{NumBatch(x)}:}{Return the number of batches \code{B}.}
#' \item{\code{NumPair(x)}:}{Return the number of batch-condition pairs \code{P}.}
#' \item{\code{NumSample(x)}:}{Return the number of samples \code{S}.}
#' \item{\code{BatchInd(x)}:}{Return the batch indices for all the cells.}
#' \item{\code{CondInd(x)}:}{Return the condition indices for all the cells.}
#' \item{\code{PairInd(x)}:}{Return the pair indices for all the cells.}
#' \item{\code{SampleInd(x)}:}{Return the sample indices for all the cells.}
#'  }
#'
#' @section Adjust the running parameters of MCEM algorithm:
#' \describe{
#' \item{\code{MiniBatch(x)}:}{Return the mini-batch size of the MCEM algorithm.}
#' \item{\code{MiniBatch(x)<-value}:}{Set the mini-batch size of the MCEM algorithm.}
#' \item{\code{EarlyStop(x)}:}{Return the number of iterations allowing early stopping in the MCEM algorithm.}
#' \item{\code{EarlyStop(x)<-value}:}{Set the number of iterations allowing early stopping in the MCEM algorithm.}
#' \item{\code{CheckLike(x)}:}{Return the number of iterations checking the increase of logliklihood function in the MCEM algorithm.}
#' \item{\code{CheckLike(x)<-value}:}{Set the number of iterations checking the increase of logliklihood function in the MCEM algorithm.}
#' \item{\code{MCRep(x)}:}{Return the number of replicates in the Monte Carlo E step.}
#' \item{\code{MCRep(x)<-value}:}{Set the number of replicates in the Monte Carlo E step.}
#' \item{\code{Ncores(x)}:}{Return the number of cores for parallel computing in the MCEM algorithm.}
#' \item{\code{Ncores(x)<-value}:}{Set the number of cores for parallel computing in the MCEM algorithm.}
#' \item{\code{Seed(x)}:}{Return the seed of the MCEM algorithm.}
#' \item{\code{Seed(x)<-value}:}{Set the seed of the MCEM algorithm.}
#' }
#'
#' @section Estimated parameter values or results of posterior inference:
#' \describe{
#' \item{\code{BIC_DIFseq(x)}:}{Return the Bayesian Informtion Criterion (BIC) value.}
#' \item{\code{Prop(x)}:}{Return a \code{S} by \code{K} matrix of the estimated sample-specific cell type proportions.}
#' \item{\code{CellTypes(x)}:}{Return a vector of length \code{N} indicating the estimated cell type labels for all cells.}
#' \item{\code{Baseline(x)}:}{Return a vector of length \code{G} giving the estimated baseline log-scale gene expression levels.}
#' \item{\code{TypeEffects(x)}:}{Return a \code{G} by \code{K} matrix of the estimated cell type effects compared with the first cell type.
#' Note that the first column is zero as the first cell type is taken as the reference cell type.}
#' \item{\code{CondEffects(x)}:}{Return a \code{G} by (\code{K * T}) matrix of the estimated condition effects compared with the first condition.
#' Note that the K columns are zero as the first condition is taken as the reference condition.}
#' \item{\code{BatchEffects(x)}:}{Return a \code{G} by \code{B} matrix of the estimated batch effects compared with the first batch.
#' Note that the first column is zero as the first batch is taken as the reference batch.}
#' \item{\code{OverDisp(x)}:}{Return a \code{G} by \code{B} matrix of the estimated over-dispersion.}
#' \item{\code{CellEffects(x)}:}{Return a vector of length \code{G} of the estimated cell size factors.}
#' \item{\code{DropoutCoef(x)}:}{Return a \code{B} by \code{2} matrix of the estimated dropout logistic coefficients.}
#' \item{\code{ImputedCounts(x)}:}{Return a \code{G} by \code{N} matrix of the imputed read counts.}
#' \item{\code{CorrectedCounts(x)}:}{Return a \code{G} by \code{N} matrix of the corrected read counts after removing batch effects and cell-size factor.}
#' }
#'
#' @section Differential inference results:
#' \describe{
#' \item{\code{Intrinsic(x)}:}{Return a boolen vector of length \code{G} indicating the identified intrinsic genes across different cell types.}
#' \item{\code{DiffExp(x)}:}{Return a \code{G} by \code{K} matrix of the identified cell-type-specific DE genes across different conditions.}
#' \item{\code{DiffAbund(x)}:}{Return the statistics and log10 p-values in the differential abundance analysis.}
#' }
#'
#'
#' \describe{
#' \item{\code{show(x)}:}{Print information to the console regarding the \code{\linkS4class{DIFseq}} object.}
#' }
#'
#' @author Fangda Song, Kevin Y. Yip and Yingying Wei
#'
#' @name DIFseq-methods
#' @rdname methods
#' @docType methods
#' @aliases
#' NumCond
#' NumCond,DIFseq-method
#' NumBatch
#' NumBatch,DIFseq-method
#' NumPair
#' NumPair,DIFseq-method
#' NumSample
#' NumSample,DIFseq-method
#' CondInd
#' CondInd,DIFseq-method
#' BatchInd
#' BatchInd,DIFseq-method
#' PairInd
#' PairInd,DIFseq-method
#' SampleInd
#' SampleInd,DIFseq-method
#' MiniBatch
#' MiniBatch<-
#' MiniBatch,DIFseq-method
#' MiniBatch<-,DIFseq-method
#' EarlyStop
#' EarlyStop<-
#' EarlyStop,DIFseq-method
#' EarlyStop<-,DIFseq-method
#' CheckLike
#' CheckLike<-
#' CheckLike,DIFseq-method
#' CheckLike<-,DIFseq-method
#' MCRep
#' MCRep<-
#' MCRep,DIFseq-method
#' MCRep<-,DIFseq-method
#' Ncores
#' Ncores<-
#' Ncores,DIFseq-method
#' Ncores<-,DIFseq-method
#' Seed
#' Seed<-
#' Seed,DIFseq-method
#' Seed<-,DIFseq-method
#' BIC_DIFseq
#' BIC_DIFseq,DIFseq-method
#' Prop
#' Prop,DIFseq-method
#' CellTypes
#' CellTypes,DIFseq-method
#' Baseline
#' Baseline,DIFseq-method
#' TypeEffects
#' TypeEffects,DIFseq-method
#' CondEffects
#' CondEffects,DIFseq-method
#' BatchEffects
#' BatchEffects,DIFseq-method
#' OverDisp
#' OverDisp,DIFseq-method
#' CellEffects
#' CellEffects,DIFseq-method
#' DropoutCoef
#' DropoutCoef,DIFseq-method
#' ImputedCounts
#' ImputedCounts,DIFseq-method
#' CorrectedCounts
#' CorrectedCounts,DIFseq-method
#' Intrinsic
#' Intrinsic,DIFseq-method
#' DiffExp
#' DiffExp,DIFseq-method
#' DiffAbund
#' DiffAbund,DIFseq-method
#' show
#' show,DIFseq-method
#'
#' @examples
#' show(DIFseq_obj)
NULL

#' @export
setMethod("NumCond", "DIFseq", function(x) {
  if (!is.null(x$Condition)) {
    return(length(unique(x$Condition)))
  } else{
    warning("Condition indices do not exist!")
  }
})

#' @export
setMethod("NumBatch", "DIFseq", function(x) {
  if (!is.null(x$Batch)) {
    return(length(unique(x$Batch)))
  } else{
    warning("Batch indices do not exist!")
  }
})

#' @export
setMethod("NumPair", "DIFseq", function(x) {
  if (!is.null(x$Batch) & !is.null(x$Condition)) {
    return(length(unique(paste0(x$Batch, "_", x$Condition))))
  } else{
    warning("Batch or condition indices do not exist!")
  }
})

#' @export
setMethod("NumSample", "DIFseq", function(x) {
  if (!is.null(x$Sample)) {
    return(length(unique(x$Sample)))
  } else{
    warning("Sample indices do not exist!")
  }
})

#' @export
setMethod("CondInd", "DIFseq", function(x) {
  if (!is.null(x$Condition)) {
    return(factor(x$Condition))
  } else{
    warning("Condition indices do not exist!")
  }
})

#' @export
setMethod("BatchInd", "DIFseq", function(x) {
  if (!is.null(x$Batch)) {
    return(factor(x$Batch))
  } else{
    warning("Batch indices do not exist!")
  }
})

#' @export
setMethod("PairInd", "DIFseq", function(x) {
  if (!is.null(x$Pair)){
    return(factor(x$Pair))
  }else if(!is.null(x$Batch) & !is.null(x$Condition)) {
    return(factor(paste0(x$Batch, "_", x$Condition)))
  } else{
    warning("Batch, condition or pair indices do not exist!")
  }
})

#' @export
setMethod("SampleInd", "DIFseq", function(x) {
  if (!is.null(x$Sample)) {
    return(factor(x$Sample))
  } else{
    warning("Sample indices do not exist!")
  }
})

#' @export
setMethod("MiniBatch", "DIFseq", function(x)
  x@iter$mini_batch)

#' @export
setMethod("MiniBatch<-", "DIFseq", function(x, value) {
  if (is.numeric(value) && length(value) == 1) {
    x@iter$mini_batch <- round(value)
  } else{
    warning("Please input an integer as the mini-batch size.")
  }
  # validObject(x)
  x
})

#' @export
setMethod("EarlyStop", "DIFseq", function(x)
  x@iter$early_stop)

#' @export
setMethod("EarlyStop<-", "DIFseq", function(x, value) {
  if (is.numeric(value) && length(value) == 1) {
    x@iter$early_stop <- round(value)
  } else{
    message("Please input an integer as the times of non-increasing
            log-likelihood to stop the algorithm.")
  }
  x
})

#' @export
setMethod("CheckLike", "DIFseq", function(x)
  x@iter$check_per_iter)

#' @export
setMethod("CheckLike<-", "DIFseq", function(x, value) {
  if (is.numeric(value) && length(value) == 1) {
    x@iter$check_per_iter <- round(value)
  } else{
    message("Please input an integer as the iteration number for a loglikelihood checking.")
  }
  # validObject(x)
  x
})

#' @export
setMethod("MCRep", "DIFseq", function(x)
  x@iter$MC_rep)

#' @export
setMethod("MCRep<-", "DIFseq", function(x, value) {
  if (is.numeric(value) && length(value) == 1) {
    x@iter$MC_rep <- round(value)
  } else{
    warning("Please input an integer as the number of replications in Monte Carlo E step.")
  }
  # validObject(x)
  x
})

#' @export
setMethod("Ncores", "DIFseq", function(x)
  x@iter$n.cores)

#' @export
setMethod("Ncores<-", "DIFseq", function(x, value) {
  if (is.numeric(value) && length(value) == 1) {
    x@iter$n.cores <- round(value)
  } else{
    warning("Please input an integer as the number of cores used in the parallel computing.")
  }
  # validObject(x)
  x
})

#' @export
setMethod("Seed", "DIFseq", function(x)
  x@iter$seed)

#' @export
setMethod("Seed<-", "DIFseq", function(x, value) {
  if (is.numeric(value) && length(value) == 1) {
    x@iter$seed <- round(value)
  } else{
    warning("Please input an integer as the seed of RNG.")
  }
  # validObject(x)
  x
})

#' @export
setMethod("BIC_DIFseq", "DIFseq", function(x) {
  if (!is.null(x@estimation$BIC)) {
    return(x@estimation$BIC)
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to compute BIC.")
  }
})

#' @export
setMethod("Prop", "DIFseq", function(x) {
  .pi <- x@estimation$pi
  if (!is.null(.pi)) {
    if (is.factor(x$Sample)) {
      rownames(.pi) <- levels(x$Sample)
    } else{
      rownames(.pi) <- levels(factor(x$Sample))
    }
    return(.pi)
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to estimate cell type proportions.")
  }
})

#' @export
setMethod("CellTypes", "DIFseq", function(x) {
  .w <- x@estimation$W
  if (!is.null(.w)) {
    return(.w)
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to estimate cell type labels of all the cells")
  }
})

#' @export
setMethod("Baseline", "DIFseq", function(x) {
  .alpha <- x@estimation$alpha
  if (!is.null(.alpha)) {
    return(.alpha)
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to estimate baseline gene expression levels.")
  }
})

#' @export
setMethod("TypeEffects", "DIFseq", function(x) {
  .beta <- x@estimation$beta
  if (!is.null(.beta)) {
    return(.beta)
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to estimate cell type effects.")
  }
})

#' @export
setMethod("CondEffects", "DIFseq", function(x) {
  .eta <- x@estimation$eta
  if (!is.null(.eta)) {
    return(.eta)
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to estimate condition effects.")
  }
})

#' @export
setMethod("BatchEffects", "DIFseq", function(x) {
  .nu <- x@estimation$nu
  if (!is.null(.nu)) {
    return(.nu)
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to estimate batch effects.")
  }
})

#' @export
setMethod("OverDisp", "DIFseq", function(x) {
  .phi <- x@estimation$phi
  if (!is.null(.phi)) {
    return(.phi)
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to estimate overdisperison parameters.")
  }
})

#' @export
setMethod("CellEffects", "DIFseq", function(x) {
  .delta <- x@estimation$delta
  if (!is.null(.delta)) {
    return(.delta)
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to estimate cell size factors.")
  }
})

#' @export
setMethod("DropoutCoef", "DIFseq", function(x) {
  .gamma <- x@estimation$gamma
  if (!is.null(.gamma)) {
    return(.gamma)
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to estimate dropout coefficients.")
  }
})

#' @importFrom SummarizedExperiment assay assays
#' @export
setMethod("ImputedCounts", "DIFseq", function(x) {
  if ("imputed_count" %in% names(assays(x))) {
    return(assay(x, i = "imputed_count"))
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to obtain imputed read counts.")
  }
})

#' @importFrom SummarizedExperiment assay assays
#' @export
setMethod("CorrectedCounts", "DIFseq", function(x) {
  if ("corrected_count" %in% names(assays(x))) {
    return(assay(x, i = "corrected_count"))
  } else{
    warning("Please run the \"DIFseq_MCEM\" function to obtain corrected read counts.")
  }
})

#' @export
setMethod("Intrinsic", "DIFseq", function(x) {
  .D <- x@estimation$intri
  if (!is.null(.D)) {
    return(.D)
  } else{
    warning("Please run the \"DiffExpression\" function to conduct differetial expression analysis.")
  }
})

#' @export
setMethod("DiffExp", "DIFseq", function(x) {
  .t <- NumCond(x)
  if(.t > 1){
    .E <- x@estimation$DE
    if (!is.null(.E)) {
      return(.E)
    } else{
      warning("Please run the \"DiffExpression\" function to conduct differetial expression analysis.")
    }
  }else{
    warning("More than one condition are required to identify cell-type-specific DE genes.")
  }
})

# DiffAbund
#' @export
setMethod("DiffAbund", "DIFseq", function(x) {
  .DA <- x@estimation$DA
  if (!is.null(.DA)) {
    return(.DA)
  } else{
    warning("Please run the \"DiffAbundance\" function to conduct differential abundance analysis.")
  }
})

#' @importFrom S4Vectors coolcat
#' @importFrom methods callNextMethod
.DIFseq_show <- function(object) {
  callNextMethod()
  # coolcat("batch: %d\n", NumBatch(x))
  cat(sprintf("batch, condition and sample number: %d, %d, %d\n", NumBatch(object), NumCond(object), NumSample(object)))
  cat(sprintf("posterior inference: %s\n", length(object@estimation) > 0))
  cat(sprintf("differential expression: %s\n", !is.null(object@estimation$DE)))
  cat(sprintf("differential abundance: %s\n", !is.null(object@estimation$DA)))

}

#' @export
#' @import methods
setMethod("show", "DIFseq", .DIFseq_show)


#' Implement MCEM algorithm for the posterior inference of the DIFseq model
#'
#' The function \code{DIFseq_MCEM} implements an MCEM algorithm to fit the model of DIFferential inference for scRNA-seq data.
#' DIFseq is an interpretable Bayesian hierarchical model that closely follows the data-generating mechanism of
#' scRNA-seq experiments from multiple conditions. DIFseq can simultaneously correct batch effects, cluster cell types,
#' impute missing data caused by dropout events and detect differentially expressed genes without requiring a preliminary normalization step.
#' We develop an MCEM algorithm to conduct posterior inference for the BUSseq model. Here, we denote the batch number as B,
#' the gene number as G and the user-specific cell-type number as K.
#'
#' @param DIFseqObj A matrix, \code{\linkS4class{SingleCellExperiment}} or DIFseq object
#' containing feature-by-cell gene expression data.
#' @param count A string indicating the name of count data matrix in DIFseq object
#' @param working_dir directory to store the parameter values at all of the iterations
#' @param correction A Boolean value indiacting whether to compute the corrected count data
#'
#' @details
#' This function implements our proposed MCEM algorithm to conduct
#' the posterior inference for the DIFseq model.
#'
#'
#' @return A \code{\linkS4class{DIFseq}} object with the parameter estimation.
#'
#' @author Fangda Song, Kevin Y. Yip and Yingying Wei
#' @examples
#'
#' DIFseq_obj <- DIFseq_MCEM(DIFseq_obj)
#'
#' @importFrom SummarizedExperiment assay assay<-
#' @useDynLib DIFseq, .registration = TRUE
#' @name DIFseq_MCEM
#' @export
DIFseq_MCEM <- function(DIFseqObj, count = "counts",
                      working_dir = getwd(), correction = TRUE){

  if(is(DIFseqObj, "DIFseq")){
    # The input data should be a DIFseq Object created by the CreateDIFseqObject function
    read <- assay(DIFseqObj, i = count)
    N <- ncol(DIFseqObj)
    G <- nrow(DIFseqObj)

    b_infor <- factor(DIFseqObj$Batch)
    B <- length(levels(b_infor))

    t_infor <- factor(DIFseqObj$Condition)
    NumCond <- length(levels(t_infor))

    s_infor <- factor(DIFseqObj$Sample)
    S <- length(levels(s_infor))

    p_infor <- paste0("B",as.numeric(b_infor),"_C",as.numeric(t_infor))
    P <- length(unique(p_infor))

    BT_pair <- matrix(NA, P, 3)
    pair_cell <- table(p_infor)
    pair_name <- names(pair_cell)
    pair_name <- strsplit(pair_name, split = "_")
    pair_name <- sapply(pair_name, substring, first = 2)
    BT_pair[,1:2] <- as.numeric(t(pair_name))
    BT_pair[,3] <- as.numeric(pair_cell)


  }else{
    stop(paste0("Please create a DIFseqObj object by \"CreateDIFseqObject\" function!\n"))
  }

  K <- DIFseqObj@n.celltype

  # Directory to record the parameter values at each iteration on the hard disk
  if(!dir.exists(working_dir)){
    dir.create(working_dir)
  }
  iter_dir <- paste0(working_dir,"/MCEM_K",K)
  dir.create(iter_dir, showWarnings=FALSE)

  t.start <- Sys.time()
  message("   conducting MCEM algorithm...\n")

  # Initial cell type label W
  # w_inital <- .initial_w(count = read, n.celltype = K, batch = b_infor)
  w_inital <- sample(1:K, size = N, replace = TRUE)

  # prepare the input to C++ program
  dim <- c(N, S, G, B, NumCond, P, K, t(BT_pair))

  control_genes <- rep(0, G)

  cell_ind <- cbind(b_infor, t_infor, factor(p_infor), s_infor)

  iter_infor <- unlist(DIFseqObj@iter)
  # names(iter_infor) <- NULL

  hyper <- unlist(DIFseqObj@hyperparam)
  # names(hyper) <- NULL

  # Allocate memory for parameters
  alpha <- rep(0, G)
  beta <- rep(0, G * K)
  eta <- rep(0, G * K * NumCond)
  nu <- rep(0, G * B)
  delta <- rep(0, N)
  gamma <- rep(0, B * 2)
  phi <- rep(0, G * B)
  pi <- rep(0, S * K)
  ptau1 <- rep(0, 3)

  R <- DIFseqObj@iter$MC_rep
  w_MC <- rep(0, N * R)

  PrL <- rep(0, G * K)
  PrJ <- rep(0, G * K * NumCond)
  iter_stage <- rep(0, 3)

  MCEM_results<-.C("DIFseq_MCEM",
                  # count data
                  y_obs = as.integer(t(read)),
                  # dimension information
                  dim = as.integer(dim),
                  # cell index,
                  cell_ind = as.integer(cell_ind),
                  # Index for control genes
                  control_genes = as.integer(control_genes),
                  # iteration setting
                  iter_infor = as.double(iter_infor),
                  # output directory
                  dir_output = as.character(iter_dir),
                  # hyperparameter
                  hyper = as.double(hyper),
                  # initial cell type labels,
                  W = as.integer(w_inital),
                  # collect outputs
                  alpha = as.double(alpha),
                  beta_vec = as.double(beta),
                  eta_vec = as.double(eta),
                  nu_vec = as.double(nu),
                  delta = as.double(delta),
                  gamma_vec = as.double(gamma),
                  phi_vec = as.double(phi),
                  pi_vec = as.double(pi),
                  ptau1 = as.double(ptau1),
                  w_MC_vec = as.integer(w_MC),
                  PrL_vec = as.double(PrL),
                  PrJ_vec = as.double(PrJ),
                  iter_stage = as.integer(iter_stage),
                  # x_imputed
                  x_imputed_vec = as.integer(rep(0,N * G)),
                  # kappa = as.double(0),
                  loglike = as.double(rep(0,2))
  )

  t.end <- Sys.time()
  message(paste0("   The MCEM algorithm takes: ",
                 round(difftime(t.end, t.start,units="mins"), 3), " mins", "\n"))

  # Build the model estimation slot
  estimation <- list()
  estimation$W <- MCEM_results$W
  estimation$alpha <- MCEM_results$alpha
  estimation$beta <- matrix(MCEM_results$beta_vec, nrow = G, byrow = TRUE)
  estimation$eta <- matrix(MCEM_results$eta_vec, nrow = G, byrow = TRUE)
  estimation$nu <- matrix(MCEM_results$nu_vec, nrow = G, byrow = TRUE)
  estimation$delta <- MCEM_results$delta
  estimation$gamma <- matrix(MCEM_results$gamma_vec, nrow = B, byrow = TRUE)
  estimation$phi <- matrix(MCEM_results$phi_vec, nrow = G, byrow = TRUE)
  estimation$pi <- matrix(MCEM_results$pi_vec, nrow = S, byrow = TRUE)
  estimation$p_beta <- MCEM_results$ptau1[1]
  estimation$p_eta <- MCEM_results$ptau1[2]
  estimation$tau1 <- MCEM_results$ptau1[3]

  # MC sampling of cell type labels and the posterior probability of indicators
  estimation$w_MC <- matrix(MCEM_results$w_MC_vec, nrow = N, byrow = TRUE)
  estimation$logPrL0 <- matrix(MCEM_results$PrL_vec, nrow = G, byrow = TRUE)
  estimation$logPrJ0 <- matrix(MCEM_results$PrJ_vec, nrow = G, byrow = TRUE)

  # Store no. of iterations
  estimation$iter <- MCEM_results$iter_stage

  # Store imputed data in assay slot as "imputed_count"
  x_imputed <- matrix(MCEM_results$x_imputed, G, N, byrow = TRUE)
  assay(DIFseqObj, i = "imputed_count") <- x_imputed

  # Store the threshold for DE and BIC
  estimation$kappa = MCEM_results$kappa
  estimation$BIC = MCEM_results$loglike[2]

  DIFseqObj@estimation <- estimation

  # Correct read count data and store at "corrected_count" assay
  if(correction){
    message(paste0("   Correcting batch effects and cell size factor to obtain",
                   " the corrected read count data...\n"))
    assay(DIFseqObj, i = "corrected_count") <-
      .adjusted_values(x_imputed, cell_ind, K, estimation)
  }

  return(DIFseqObj)
}

# Quantile normalization for batch effects correction
#' @importFrom stats pnbinom runif qnbinom
.adjusted_values <- function(ReadCount, Indicators, .K,
                            Estimation){

  .alpha = Estimation$alpha
  .beta = Estimation$beta
  .eta = Estimation$eta
  .nu = Estimation$nu
  .delta = Estimation$delta
  .phi = Estimation$phi
  .w = Estimation$W

  CorrectedCount <- ReadCount
  .N <- ncol(ReadCount)
  .G <- nrow(ReadCount)

  for(i in 1:.N){

    b <- Indicators[i,1]
    t <- Indicators[i,2]
    w <- .w[i]

    # percentile
    px <- pnbinom(ReadCount[,i], size = .phi[,b], mu = exp(.alpha + .beta[,w] + .eta[,(t-1) * .K + w] + .nu[,b] + .delta[i]))
    pxminus1 <- pnbinom(ReadCount[,i] - 1, size = .phi[,b], mu = exp(.alpha + .beta[,w] + .eta[,(t-1) * .K + w] + .nu[,b] + .delta[i]))

    # get the aligned percentile
    local_u <- runif(.G) * (px - pxminus1) + pxminus1
    local_u <- ifelse(local_u > 0.9999, 0.9999, local_u)

    # obtain the quantile
    CorrectedCount[,i] <- qnbinom(local_u, size = .phi[,1], mu = exp(.alpha + .beta[,w] + .eta[,(t-1) * .K + w]))

    if(i %% 1000 == 0){
      print(paste("Finish the correction of", i, "cells..."))

    }
  }
  return(CorrectedCount)
}


# Differential expression
#' Identify intrinsic genes and cell-type-specific DE genes
#'
#' The function \code{DiffExpression} not only controls the Bayesian false discovery rate (FDR) to adjust for multiple hypothesis tests,
#' but also controls the minimum proportions of non-zero counts in a group of cells to ensure the biological meaning.
#'
#'
#' @param DIFseqObj A matrix, \code{\linkS4class{SingleCellExperiment}} or DIFseq object
#' containing feature-by-cell gene expression data.
#' @param fdr Bayesian false discovery rate of all intrinsic genes and cell-type-specific DE genes.
#' @param min.prop Minimum proportions of cells with non-zero read counts in a given group.
#' For intrinsic genes, it requires the proportions of non-zero counts in a cell type for a given gene to be greater than the threshold,
#'
#' @details
#' In differetianl expression analysis, there are `G * (K-1)` hypothesis tests to identify intrinsic genes, where `G` and `K` denote the number of genes and cell types.
#' Meanwhile, there are `G * (T-1) * K` hypothesis tests for cell-type-specific DE genes, where `T` represents the number of conditions.
#' Thus, we control the Bayesian FDR of all the `G * (T * K - 1)` tests to adjust for multiple hypotheses.
#'
#' @return A DIFseq object storing the identified intrinsic genes and cell-type-specific DE genes in the `estimation` slot
#'
#' @author
#' Fangda Song, Kevin Y. Yip and Yingying Wei
#'
#' @examples
#'
#' DIFseq_obj <- DiffExpression(DIFseq_obj)
#' DE_genes <- DiffExp(DIFseq_obj)
#'
#' @importFrom SummarizedExperiment assay
#' @name DiffExpression
#' @export
DiffExpression <- function(DIFseqObj, fdr = 0.05, min.prop = 0.1){

  if(length(DIFseqObj@estimation) > 0){

    # estiamted cell type labels
    .W <- DIFseqObj@estimation$W

  }else{
    stop(paste0("Please run \"DIFseq_MCEM\" function to conduct posterior inference!\n"))
  }

  # Obtain dimension infor
  .K <- DIFseqObj@n.celltype

  t_infor <- factor(DIFseqObj$Condition)
  .NumCond <- length(levels(t_infor))

  read <- assay(DIFseqObj, i = "counts")
  .G <- nrow(read)



  if(.NumCond > 1){
    # Output intrinsic genes and cell-type-specific DE genes
    # Extract log(Pr(L = 0)) and log(Pr(J = 0))
    PrL_est <- DIFseqObj@estimation$logPrL0
    PrJ_est <- DIFseqObj@estimation$logPrJ0

    PrL_est <- exp(PrL_est)
    PrJ_est <- exp(PrJ_est)

    # colnames(PrL_est) <- paste0("K",1:K_opt)
    # colnames(PrJ_est) <- paste0("T",rep(1:Num_Treatment, each = K_opt),"K",rep(1:K_opt, Num_Treatment))

    pval_col <- cbind(PrL_est[,-1], PrJ_est[, .K + 1:((.NumCond - 1) * .K)])

    xi_thres <- .postprob_DE_thr_fun(pval_col, fdr_threshold = fdr)
    DE_est <- .estimate_IG_indicators(pval_col, xi_thres)

    D_est <- apply(DE_est[,1:(.K-1)], 1, sum) > 0
    num_intri <- sum(D_est)

    E_est <- apply(DE_est[,.K - 1 + 1:((.NumCond - 1) * .K)], 1, sum) > 0
    num_DE <- sum(E_est)

    ################################################
    # consider the expression proportions of genes #
    ################################################
    ## intrinsic genes
    L_est <- DE_est[,1:(.K - 1)]
    L_adj <- matrix(0, .G, .K - 1)

    expr_prop <- array(NA, dim = c(.G, .K, 2))

    for(k in 2:.K){
      celltype_k <- which(.W == k)
      for(g in which(L_est[,k-1] == 1)){

        yg <- read[g, ]

        ygk <- yg[celltype_k]
        ygnk <- yg[-celltype_k]

        expr_prop[g,k,1] <- sum(ygk > 0)/length(ygk)
        expr_prop[g,k,2] <- sum(ygnk > 0)/length(ygnk)

        if(expr_prop[g,k,1] > min.prop & expr_prop[g,k,2] > min.prop){
          L_adj[g, k - 1] <- 1
        }
      }

    }

    D_adj <- apply(L_adj, 1, sum) > 0
    num_intri_adj <- sum(D_adj)

    message(num_intri_adj, " intrinsic genes are found.\n")
    DIFseqObj@estimation$intri <- D_adj
    # message("The output format is a vector implying the intrinsic gene",
    #        " indices.\n")

    # cell-type-specific DE genes
    J_est <- DE_est[,.K - 1 + 1:((.NumCond - 1) * .K)]
    J_adj <- matrix(0, .G, (.NumCond - 1) * .K)

    expr_prop_cond <- array(NA, dim = c(.G, .NumCond,.K, 2))

    for(t in 2:.NumCond){
      for(k in 1:.K){
        celltype_t <- which(.W == k & t_infor == t)
        celltype_other <- which(.W == k & t_infor != t)
        for(g in which(J_est[,k + (t-2) * .K] == 1)){

          yg <- read[g, ]
          ygt <- yg[celltype_t]
          ygnt <- yg[celltype_other]

          expr_prop_cond[g,t,k,1] <- sum(ygt > 0)/length(ygt)
          expr_prop_cond[g,t,k,2] <- sum(ygnt > 0)/length(ygnt)
          # print(expr_prop_treatment[g,t,k,])

          if(expr_prop_cond[g,t,k,1] > min.prop & expr_prop_cond[g,t,k,2] > min.prop){
            J_adj[g, k + (t-2) * .K] <- 1
          }
        }
      }
    }

    E_adj <- apply(J_adj, 1, sum) > 0
    num_DE_adj <- sum(E_adj)

    E_combinded <- J_adj[,1:.K]
    if(.NumCond > 2){
      for(t in 3:.NumCond){
        E_combinded <- E_combinded | J_adj[,1:.K + .K * (t-2)]
      }
    }

    num_intri_adj_acorss_cond <- apply(E_combinded,2,sum)

    message(paste(num_intri_adj_acorss_cond, collapse = ", "),
            " cell-type-specific DE genes are found in ",.K," cell types, respectively.\n")

    DIFseqObj@estimation$DE <- ifelse(E_combinded == 1, TRUE, FALSE)

  }else{
    # Only output intrinsic genes
    # Extract log(Pr(L = 0))
    PrL_est <- DIFseqObj@estimation$logPrL0
    PrL_est <- exp(PrL_est)

    pval_col <- cbind(PrL_est[,-1])

    xi_thres <- .postprob_DE_thr_fun(pval_col, fdr_threshold = fdr)
    DE_est <- .estimate_IG_indicators(pval_col, xi_thres)

    D_est <- apply(DE_est[,1:(.K-1)], 1, sum) > 0
    num_intri <- sum(D_est)

    ################################################
    # consider the expression proportions of genes #
    ################################################
    ## intrinsic genes
    L_est <- DE_est[,1:(.K - 1)]
    L_adj <- matrix(0, .G, .K - 1)

    expr_prop <- array(NA, dim = c(.G, .K, 2))

    for(k in 2:.K){
      celltype_k <- which(.W == k)
      for(g in which(L_est[,k-1] == 1)){

        yg <- read[g, ]

        ygk <- yg[celltype_k]
        ygnk <- yg[-celltype_k]

        expr_prop[g,k,1] <- sum(ygk > 0)/length(ygk)
        expr_prop[g,k,2] <- sum(ygnk > 0)/length(ygnk)

        if(expr_prop[g,k,1] > min.prop & expr_prop[g,k,2] > min.prop){
          L_adj[g, k - 1] <- 1
        }
      }

    }

    D_adj <- apply(L_adj, 1, sum) > 0
    num_intri_adj <- sum(D_adj)

    message(num_intri_adj, " intrinsic genes are found.\n")
    DIFseqObj@estimation$intri <- D_adj
  }

  return(DIFseqObj)
}

.fdrDEindicator <- function(xi, kappa){

  ind_intr <- xi <= kappa
  fdr <- sum(xi[ind_intr])/sum(ind_intr)

  return(fdr)
}

# Calculate the DE posterior probability threshold
.postprob_DE_thr_fun <- function(xi, fdr_threshold=0.05){

  kappa_fdr_matr <- NULL
  kappa_set <- sort(unique(xi))

  kappa_ind <- which(kappa_set < 0.5 & kappa_set > fdr_threshold)

  if(length(kappa_ind) == 0){
    kappa <- 0.5
  }else{

    for(i in kappa_ind){

      kappa <- kappa_set[i]
      fdr <- .fdrDEindicator(xi, kappa=kappa)

      if(fdr > fdr_threshold){
        break
      }
    }

    kappa <- kappa_set[i-1]
  }

  return(kappa)
}

# Estimate intrinsic gene indicators
.estimate_IG_indicators <- function(xi, postprob_DE_threshold = 0.5){

  EstL <- xi
  EstL[xi >= postprob_DE_threshold] <- 0
  EstL[xi <= postprob_DE_threshold] <- 1
  # message("The output format is a matrix.\n")
  # message(paste0("Each row represents a gene, and each column",
  #               " corresponds to a cell type from 2 to K\n"))
  return(EstL)
}

# # Intrinsic gene index
# .IG_index <- function(EstIGindicators){
#   ind <- which(rowSums(EstIGindicators) > 0)
#   message(c(length(ind), " intrinsic genes are found.\n"))
#   message("The output format is a vector implying the intrinsic gene",
#           " indices.\n")
#   return(ind)
# }

# Differential abundance
#' Conduct differential abundance analysis
#'
#' The function \code{DiffAbundance} apply rigorous hypothesis testing to
#' identify differential abundance (DA) across batches or conditions
#'
#' @param DIFseqObj A matrix, \code{\linkS4class{SingleCellExperiment}} or DIFseq object
#' containing feature-by-cell gene expression data.
#' @param ref An integer indicates the reference cell type. By default, the last cell type is regarded as the reference cell type
#' @param dim Dimension to be test. It has three possible options: \code{Batch}, \code{Cond} and `Pair`.
#' @param subset A numeric vector indicating the batches/conditions to be included in the comparison. By default,
#' consider all of the batches, conditions or pairs in differential abundance analysis
#'
#' @return A DIFseq object storing the overall and cell-type-specific test statistics
#' and p-values in the \code{estimation} slot.
#'
#' @details  The `dim` argument specifies the group of cells to conduct the DA analysis.
#' In particular, `Batch` indicates to check DA between batches involved `subset` argument over all of the conditions;
#' `Cond` means to check DA between conditions specified by `subset` argument over all of the batches;
#' similarly, `Pair` means to check DA between different pairs involved in `subset` argument.
#'
#' @author
#' Fangda Song, Kevin Y. Yip and Yingying Wei
#'
#' @examples
#'
#' DIFseq_obj <- DiffAbundance(DIFseq_obj, dim = "Cond")
#' DiffAbund(DIFseq_obj)
#'
#' DIFseq_obj <- DiffAbundance(DIFseq_obj, dim = "Batch")
#' DiffAbund(DIFseq_obj)
#'
#' @importFrom stats cov pnorm pchisq
#' @name DiffAbundance
#' @export
DiffAbundance <- function(DIFseqObj, ref = NULL,
                          dim = c("Cond", "Batch", "Pair"), subset = NULL){

    .pi <- Prop(DIFseqObj)
  w <- DIFseqObj@estimation$w_MC
  meta <- DIFseqObj@colData

  if(is.null(ref)){
    ref = ncol(.pi)
  }

  .K <- ncol(.pi)
  .S <- nrow(.pi)
  s_infor <- meta$Sample

  if(!is.null(colnames(.pi))){
    .cname <- colnames(.pi)
  }else{
    .cname <- paste0("Cell Type ",1:.K)
  }

  .cname <- .cname[-ref]
  # browser()

  if(dim == "Cond"){
    # get the indices of each cell and the number of treatments
    d_infor <- factor(meta$Condition)
  }else if(dim == "Batch"){
    d_infor <- factor(meta$Batch)
  }else if(dim == "Pair"){
    d_infor <- factor(meta$Pair)
  }else{
    stop(paste0("dim should be one of 'Cond', 'Batch' and 'Pair'!\n"))
  }

  if(is.null(subset)){
    .dname <- levels(d_infor)
    .d <- length(.dname)
    .dsub <- 1:.d
  }else{

    .dsub <- subset
    .d <- length(subset)
    .dname <- levels(d_infor)[subset]

  }

  # Sample_set stores samples belonging to each treatment
  Sample_set <- list()
  for (j in 1:.d) {
    Sample_set[[j]] <- as.numeric(unique(s_infor[as.numeric(d_infor) == .dsub[j]]))
  }

  res <- matrix(NA, .d * (.d - 1) / 2, 2 * .K)
  row_names <- NULL

  # compute Louis' estimator
  lvar <- array(NA, dim = c(.S, .K - 1, .K - 1))

  # put the reference to the last column
  if(ref != .K){
    .pi <- cbind(.pi[,-ref] , .pi[,ref])
    pos_ref <- which(w == ref)
    pos_lat <- which(w > ref)
    w[pos_ref] <- .K
    w[pos_lat] <- w[pos_lat] - 1
  }

  for(s in 1:.S){
    cell_index <- which(s_infor==s)
    lvar[s,,] <- .Louis_variance_MCESM(.pi[s,],w[cell_index,])
  }

  index <- 1
  for (j1 in 1:(.d - 1)) {
    for (j2 in (j1+1):.d) {

      set1 <- Sample_set[[j1]]
      set2 <- Sample_set[[j2]]

      ns1 <- length(set1)
      ns2 <- length(set2)

      # Difference in proportion
      if(ns1 > 1){
        pi_mean_1 <- apply(.pi[set1, ], 2, mean)
      }else{
        pi_mean_1 <- .pi[set1, ]
      }

      if(ns2 > 1){
        pi_mean_2 <- apply(.pi[set2, ], 2, mean)
      }else{
        pi_mean_2 <- .pi[set2, ]
      }

      pi_dif <- pi_mean_1 - pi_mean_2

      # Compute within-group covariance
      Var_within <- matrix(0, .K - 1, .K - 1)
      for (s in set1) {
        Var_within <- Var_within + lvar[s, ,] / ns1 ^ 2
      }

      for (s in set2) {
        Var_within <- Var_within + lvar[s, ,] / ns2 ^ 2
      }

      # Compute between-group covariance
      Var_est <- Var_within

      if(ns1 > 1){
        Var_est <- Var_est + cov(.pi[set1, 1:(.K - 1)])/ns1
      }

      if(ns2 > 1){
        Var_est <- Var_est + cov(.pi[set2, 1:(.K - 1)])/ns2
      }

      # test statistics
      stat_celltype <- pi_dif[1:(.K - 1)] / sqrt(diag(Var_est))
      stat_overall <- t(pi_dif[1:(.K - 1)]) %*% solve(Var_est) %*% pi_dif[1:(.K - 1)]

      # pval
      log10p_celltype <- log10(2) + pnorm(abs(stat_celltype), lower.tail = FALSE, log.p = TRUE)/log(10)
      log10p_overall <- pchisq(stat_overall, df = .K - 1, lower.tail = FALSE, log.p = TRUE)/log(10)

      # # adjust by BH
      # p_adj <- p.adjust(p_celltype, method = "BH")

      res[index, 2 * 1:.K - 1] <- c(stat_celltype, stat_overall)
      # res[index, 2 * 1:.K] <-  c(p_adj, p_overall)
      res[index, 2 * 1:.K] <-  c(log10p_celltype, log10p_overall)

      row_names <- c(row_names, paste0(dim,"_",.dname[j1], " vs ",dim,"_",.dname[j2]))

      index <- index + 1
    }
  }


  rownames(res) <- row_names
  colnames(res) <- paste(rep(c(.cname,"Overall"),each = 2),rep(c("stat", "log10(pval)"),.K),sep = ":")
  DIFseqObj@estimation$DA <- res

  return(DIFseqObj)
}


.Louis_variance_MCESM <- function(cell_proportion, cell_labels_MC){

  if(is.list(cell_proportion)){
    cell_proportion <- unlist(cell_proportion)
  }

  n_rep <- ncol(cell_labels_MC)
  K <- length(cell_proportion)
  Km1 <- length(cell_proportion) - 1

  sum_w <- matrix(NA, K, n_rep)
  for(r in 1:n_rep){
    for(k in 1:K){
      sum_w[k,r] <- sum(cell_labels_MC[,r]==k)
    }
  }

  # Add 0.01 for singular values
  sum_w <- sum_w + 0.01

  #print(sum_w)

  first_der <- rep(0, Km1)
  first_der_cross <- matrix(0, Km1, Km1)
  second_der <- matrix(0, Km1, Km1)

  for(r in 1:n_rep){
    temp_first <- sum_w[1:(K - 1),r]/cell_proportion[1:(K-1)] - sum_w[K]/cell_proportion[K]
    first_der <- first_der + temp_first
    first_der_cross <- first_der_cross + temp_first %*% t(temp_first)
    second_der <- second_der -
      diag(sum_w[1:(K - 1),r]/(cell_proportion[1:(K-1)])^2) -
      sum_w[K,r]/(cell_proportion[K])^2
  }

  first_der <- first_der/n_rep
  first_der_cross <- first_der_cross/n_rep
  second_der <- second_der/n_rep

  obs_variance <- -second_der + first_der_cross - first_der %*% t(first_der)
  obs_variance <- solve(obs_variance)

  return(obs_variance)
}


#' @importFrom edgeR cpm
#' @importFrom cluster pam
#' @importFrom factoextra get_dist
.initial_w <- function(count, n.celltype, batch){

  G <- nrow(count)
  N <- ncol(count)
  rownames(count) <- paste0("Gene_",1:G)
  colnames(count) <- paste0("Cell_",1:N)

  # normalization
  y_norm <- cpm(count, log = TRUE)
  obs_df <- data.frame(t(y_norm))

  # scaling
  obs_df <- scale(obs_df)

  # clustering on the first batch
  cells_in_batch1 <- which(batch == 1)
  random_cells <- sample(cells_in_batch1, 1000)
  batch_df <- obs_df[random_cells, ]

  cluster_batch <- pam(batch_df, k = n.celltype)

    # medoids
    Medoids_batch <- cluster_batch$medoids
    out <- rep(NA, N)
    for (i in 1:N) {
      # L1_dist <- colSums(abs(t(Medoids_batch)- obs_df[i,]))
      Man_dist <-
        get_dist(rbind(obs_df[i, ], Medoids_batch), method = "pearson")
      out[i] <- which.min(Man_dist[1:n.celltype])
    }

    return(out)
}

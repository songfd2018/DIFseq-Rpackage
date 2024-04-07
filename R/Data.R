#' An example dataset of the DIFseq object
#'
#' The object is a DIFseq object obtained from the simulated data in the "Example" of \code{DIFseq-package}.
#'
#' @format A DIFseq object contains the following slots
#' \describe{
#'   \item{assays}{Contains gene expression matrices}
#'   \item{colData}{Contains the batch, condition and sample information}
#'
#' }
#' @details{
#' The simulated count data is generated in two batches and two conditions. Each batch-condition pair contains
#' 5 donors or samples, and each donor includes 100 cells. For all the cells, 40 genes are measured.
#' Moreover, all cells come from three cell types.
#' }
#'
#' @examples{
#' \dontrun{
#' set.seed(1234)
#' # The number of conditions
#' Num_Cond <- 2
#'
#' # The number of batches under each condition
#' B <- 2
#'
#' # The number of active pairs
#' P <- 4
#'
#' # The number of donors in each pair
#' m <- 5
#'
#' # The number of donors
#' S <- P * m
#'
#' # The number of cells per batch
#' ns <- rep(100, S)
#'
#' # The grid of observed cells
#' BT_pair <- data.frame(Batch = c(1,1,2,2),
#'                       Condition = c(1,2,1,2),
#'                       nbt = rep(ns[1] * m, P))
#'
#'
#' b_infor <- rep(BT_pair[,1], BT_pair[,3])
#' t_infor <- rep(BT_pair[,2], BT_pair[,3])
#' p_infor <- rep(1:P, BT_pair[,3])
#' s_infor <- rep(1:S, ns)
#'
#' #The total number of cells
#' N <- sum(ns)
#'
#' #The number of genes
#' G <- 40
#'
#' #The number of cell types
#' K <- 3
#'
#' # the project name
#' proj <- "demo"
#'
#' # The first column of gamma.syn denotes the intercept of
#' # the logistic regression for dropout events
#' # The second column of gamma.syn denotes the odds ratios
#' # of the logistic regression for dropout events
#' gamma.syn<-matrix(0,B,2)
#' gamma.syn[1,]<-c(-0.5,-0.5)
#' gamma.syn[2,]<-c(-0.8,-0.5)
#'
#' #the log-scale baseline expression levels
#' alpha.syn<-rep(NA,G)
#' alpha.syn[1:(G/4)]<-rep(2,G/4)
#' alpha.syn[(G/4+1):(G/4*2)]<-rep(1,G/4)
#' alpha.syn[(G/4*2+1):(G/4*3)]<-rep(0.5,G/4)
#' alpha.syn[(G/4*3+1):G]<-rep(0,G/4)
#'
#' alpha.syn[1:3] <- 4
#' alpha.syn[G/4 + 1:3] <- 3
#' alpha.syn[G/4 * 2 + 1:3] <- 2.5
#' alpha.syn[G/4 * 3 + 1:3] <- 2
#'
#' #the cell-type effects
#' beta.syn<-matrix(0,G,K)
#'
#' # The first cell type is regarded as the reference cell type
#' # without cell-type effects
#' beta.syn[,1] <- 0
#'
#' #the cell-type effects of the second cell type
#' beta.syn[1:3,2] <- -2
#' beta.syn[4:6, 2] <- 2
#' beta.syn[G/4 + 1:3,2] <- -2
#' beta.syn[G/4 + 4:6,2] <- 2
#' beta.syn[G/4 * 2 + 1:3,2] <- -2
#' beta.syn[G/4 * 2 + 4:6,2] <- 2
#' beta.syn[G/4 * 3 + 1:3,2] <- -2
#' beta.syn[G/4 * 3 + 4:6,2] <- 2
#'
#' #the cell-type effects of the third cell type
#' beta.syn[1:3,3] <- -2
#' beta.syn[7:9, 3] <- 2
#' beta.syn[G/4 + 1:3,3] <- -2
#' beta.syn[G/4 + 7:9,3] <- 2
#' beta.syn[G/4 * 2 + 1:3,3] <- -2
#' beta.syn[G/4 * 2 + 7:9,3] <- 2
#' beta.syn[G/4 * 3 + 1:3,3] <- -2
#' beta.syn[G/4 * 3 + 7:9,3] <- 2
#'
#' #####################
#' # Condition effects #
#' #####################
#' eta.syn<-matrix(NA, G, Num_Cond * K)
#' # the first condition as the reference batch
#' eta.syn[,1:K] <- 0
#'
#' # the second condition
#' # Up-regulated or Down-regulated highly expressed genes
#' # in the corresponding cell types
#' eta.syn[,K + 1:K] <- 0
#' # for the first cell type
#' eta.syn[1:2, K + 1] <- -1
#' eta.syn[G/4 + 1:2,K + 1] <- 1
#'
#' # for the second cell type
#' eta.syn[4:5, K + 2] <- -1
#' eta.syn[G/4 + 4:5,K + 2] <- 1
#'
#' #################
#' # batch effects #
#' #################
#' nu.syn<-matrix(NA,G,B)
#'
#' #the first batch is taken as the reference batch
#' #without batch effects
#' nu.syn[,1] <- 0
#'
#' #the batch effect of the second batch
#' nu.syn[,2] <- rep(c(1, -1, 2, -2),each = G/4)
#'
#' ##################################
#' # cell-specific sequencing depth #
#' ##################################
#' delta.syn <- rep(NA, N)
#' nbt <- table(p_infor)
#'
#' # the first cell in each batch-condition pair is regarded as the reference cell
#' # with the cell-specific size factors being 0
#' delta.syn[1:(nbt[1]/2)] <- 0
#' delta.syn[(nbt[1]/2 + 1):nbt[1]] <- 0.5
#'
#' # the second pair
#' delta.syn[nbt[1] + 1:(nbt[2]/4)] <- 0
#' delta.syn[nbt[1] + (nbt[2]/4+1):(nbt[2]/2)] <- 0.5
#' delta.syn[nbt[1] + (nbt[2]/2+1):nbt[2]] <- 1
#'
#' # the third pair
#' delta.syn[sum(nbt[1:2]) + 1:(nbt[3]/4)] <- 0
#' delta.syn[sum(nbt[1:2]) + (nbt[3]/4+1):(nbt[3]/2)] <- -0.5
#' delta.syn[sum(nbt[1:2]) + (nbt[3]/2+1):nbt[3]] <- -1
#'
#' # the forth pair
#' delta.syn[sum(nbt[1:3]) + 1:(nbt[4]/2)] <- 0
#' delta.syn[sum(nbt[1:3]) + (nbt[4]/2+1):nbt[4]] <- -1
#'
#' ##############################################################
#' # batch-specific and gene-specific overdispersion parameters #
#' ##############################################################
#' phi.syn<-matrix(10, G, B)#mean 2 var 0.5
#' phi.syn[(G*0.4 + 1):(G * 0.8),] <- 3
#' phi.syn[(G*0.8 + 1):G,] <- 1
#'
#' rdir <- function(n, xi){
#'   .K <- length(xi)
#'   res <- matrix(NA, n, .K)
#'   for(i in 1:n){
#'     gam <- rgamma(.K, xi)
#'     res[i,] <- gam / sum(gam)
#'   }
#'   return(res)
#' }
#'
#' # the cell-type proportions in each batch
#' pi.syn <- matrix(NA, S, K)
#'
#' xi.syn <- matrix(NA, P, K)
#' xi.syn[1,] <- c(0.3,0.5,0.2)
#' xi.syn[2,] <- c(0.3,0.4,0.3)
#' xi.syn[3,] <- c(0.3,0.5,0.2)
#' xi.syn[4,] <- c(0.3,0.4,0.3)
#' xi.syn <- xi.syn * 1000
#'
#' s_ind <- 0
#' for(p in 1:P){
#'   pi.syn[s_ind + 1:m, ] <- rdir(m, xi.syn[p, ])
#'   s_ind <- s_ind + m
#' }
#'
#' ##############################################
#' # Simulate Latent Varibles and Observed data #
#' ##############################################
#' # the cell-type indicators of each cell
#' w <- NULL
#'
#' for(s in 1:S){
#'   w <- c(w, apply(rmultinom(ns[s],1,pi.syn[s,]),2,function(x){which(x==1)}))
#' }
#'
#' # the indicators for dropout events
#' z <- matrix(NA, G, N)
#'
#' # the underlying true expression levels
#' x <- matrix(NA, G, N)
#'
#' # the observed expression levels
#' y <- matrix(NA, G, N)
#'
#' # the logarithm of mean expreesion level of each gene in each cell
#' log.mu <- matrix(NA, G, N)
#'
#' # generate the latent variable and observed data
#' for(i in 1:N){
#'
#'   # obtain the batch index
#'   b <- b_infor[i]
#'   t <- t_infor[i]
#'   p <- p_infor[i]
#'   s <- s_infor[i]
#'
#'   log.mu[,i] <- alpha.syn + beta.syn[,w[i]] + eta.syn[,(t-1) * K + w[i]] +
#'     nu.syn[,b] + delta.syn[i]
#'
#'   x[,i]<-rnbinom(G, size = phi.syn[,b], mu = exp(log.mu[,i]))
#'
#'   prob_drop <- exp(gamma.syn[b,1] + gamma.syn[b,2] *
#'   x[,i])/(1+exp(gamma.syn[b,1] + gamma.syn[b,2] * x[,i]))
#'   z[,i] <- rbinom(G, size = 1, prob = prob_drop)
#'
#'   y[,i] <- (1-z[,i]) * x[,i]
#'
#' }
#'
#' DIFseq_obj <- CreateDIFseqObject(assay = list(counts = y),
#'                                  batch = b_infor, condition = t_infor, sample = s_infor,
#'                                  n.celltype = 3)
#' }
#' }
"DIFseq_obj"

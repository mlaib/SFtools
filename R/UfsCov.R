#' UfsCov algorithm for unsupervised feature selection
#'
#' Applies the UfsCov algorithm based on the space filling concept,
#' by using a sequatial forward search (SFS).
#' @usage UfsCov(data)
#' @param data Data of class: \code{matrix} or \code{data.frame}.
#'
#' @return A list of two elements:
#'  \itemize{
#'   \item \code{CovD} a vector containing the coverage measure of
#'   each step of the SFS.
#'   \item \code{IdR} a vector containing the added variables during
#'   the selection procedure.
#'   }
#' @author Mohamed Laib \email{Mohamed.Laib@@unil.ch}
#'
#' @note The algorithm does not deal with missing values and constant
#' features. Please make sure to remove them.
#' @details Since the algorithm is based on pairwise distances, and
#' according to the computing power of your machine, large number of
#' data points can take much time and needs more memory.
#' See \code{\link{UfsCov_par}} for parellel computing, or
#' \code{\link{UfsCov_ff}} for memory efficient storage of large data
#' on disk and fast access (by using the \code{ff} and the \code{ffbase} packages).
#'
#' @examples
#' infinity<-Infinity(n=1000)
#' Results<- UfsCov(infinity)
#'
#' cou<-colnames(infinity)
#' nom<-cou[Results[[2]]]
#' par(mfrow=c(1,1), mar=c(5,5,2,2))
#' names(Results[[1]])<-cou[Results[[2]]]
#' plot(Results[[1]] ,pch=16,cex=1,col="blue", axes = FALSE,
#' xlab = "Added Features", ylab = "Coverage measure")
#' lines(Results[[1]] ,cex=2,col="blue")
#' grid(lwd=1.5,col="gray" )
#' box()
#' axis(2)
#' axis(1,1:length(nom),nom)
#' which.min(Results[[1]])
#'
#' \dontrun{
#'
#' #### UfsCov on the Butterfly dataset ####
#' require(IDmining)
#'
#' N <- 1000
#' raw_dat <- Butterfly(N)
#' dat<-raw_dat[,-9]
#'
#' Results<- UfsCov(dat)
#' cou<-colnames(dat)
#' nom<-cou[Results[[2]]]
#' par(mfrow=c(1,1), mar=c(5,5,2,2))
#' names(Results[[1]])<-cou[Results[[2]]]
#'
#' plot(Results[[1]] ,pch=16,cex=1,col="blue", axes = FALSE,
#' xlab = "Added Features", ylab = "Coverage measure")
#' lines(Results[[1]] ,cex=2,col="blue")
#' grid(lwd=1.5,col="gray" )
#' box()
#' axis(2)
#' axis(1,1:length(nom),nom)
#' which.min(Results[[1]])
#'
#' }
#' @references
#' M. Laib and M. Kanevski (2017). Unsupervised Feature Selection Based on Space
#' Filling Concept, \href{https://arxiv.org/abs/1706.08894}{arXiv:1706.08894}.
#'
#' @import wordspace
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export

UfsCov<-function(data){
  if (!is.matrix(data) & !is.data.frame(data)) {
    stop('X must be a matrix or a data frame.')
  }
  so<-so1<-as.matrix(1:ncol(data))
  sf<-0
  ind<-vlu<-c()
  dat<- apply(data , MARGIN = 2,FUN = function(X) (X - min(X))/diff(range(X)))
  bprog <- txtProgressBar(min = 0, max = ncol(data))

  for (i in 1:ncol(data)){
    disA<-apply(so, 1,FUN=function(a) (measr(dat[,c(a,sf)])))
    #if (all(is.na(disA))==TRUE) {disA<-rep(0,length(disA))}
    ind<-c(ind,so[which.min(disA)])
    vlu<-c(vlu,min(disA, na.rm = TRUE))
    sf<-ind

    so<-as.matrix(so1[-sf])
    setTxtProgressBar(bprog,i)

    }

  return(list(CovD=vlu, IdR=sf))}


measr<-function (design, verbose=F) {
  X <- as.matrix(design)
  n <- nrow(X)
  Distance <- dist.matrix(X, method = "euclidean")
  diag(Distance) <- 1000
  Dmin <- as.matrix(apply(Distance, 2, min))
  gammabar <- mean(Dmin)
  s <- sum(apply(Dmin, 2, FUN= function(a)((a-gammabar)^2)))
  cov <- (1/gammabar) * ((1/n) * s)^(1/2)
  return(cov)}


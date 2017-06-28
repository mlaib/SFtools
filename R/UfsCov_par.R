#' UfsCov algorithm for unsupervised feature selection
#'
#' Applies the UfsCov algorithm based on the space filling concept,
#' by using a sequatial forward search (SFS).This function
#' offers a parellel computing.
#' @usage UfsCov_par(data, ncores=2)
#' @param data Data of class: \code{matrix} or \code{data.frame}.
#' @param ncores Number of cores to use (by default: \code{ncores=2}).
#' @return A list of two elements:
#'  \itemize{
#'   \item \code{CovD} a vector containing the coverage measure of
#'   each step of the SFS.
#'   \item \code{IdR} a vector containing the added variables during
#'   the selection procedure.
#'   }
#' @author Mohamed Laib \email{Mohamed.Laib@@unil.ch}
#' @note The algorithm does not deal with missing values and constant
#' features. Please make sure to remove them. Note that it is not recommanded to
#' use this function with small data, it takes more time than using the
#' standard \code{\link{UfsCov}} function.
#'
#' @details Since the algorithm is based on pairwise distances, and
#' according to the computing power of your machine, large number of
#' data points needs more memory. See \code{\link{UfsCov_ff}} for memory
#' efficient storage of large data on disk and fast access (by using the
#' \code{ff} and the \code{ffbase} packages).
#'
#'
#' @seealso \code{\link{UfsCov}}, \code{\link{UfsCov_ff}}
#'
#' @examples
#' N <- 1000
#' dat<-Infinity(N)
#' Results<- UfsCov_par(dat,ncores=2)
#'
#' cou<-colnames(dat)
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
#' N<-5000
#' dat<-Infinity(N)
#'
#' ## Little comparison:
#' system.time(Uf<-UfsCov(dat))
#' system.time(Uf.p<-UfsCov_par(dat, ncores = 4))
#'
#' }
#' @references
#' M. Laib and M. Kanevski (2017). Unsupervised Feature Selection Based on Space
#' Filling Concept, \href{https://arxiv.org/abs/1706.08894}{arXiv:1706.08894}.
#'
#' @import wordspace doParallel parallel
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'

UfsCov_par<-function(data, ncores=2){
  so<-so1<-as.matrix(1:ncol(data))
  sf<-0
  ind<-c()
  vlu<-c()
  dat<- apply(data , MARGIN = 2,FUN = function(X) (X - min(X))/diff(range(X)))
  cl <- makeCluster(ncores)
  bprog <- txtProgressBar(min = 0, max = ncol(data))
  clusterExport(cl, varlist = c("measr","dat","sf","so"), envir=environment())
  for (i in 1:ncol(data)){
    disA<-parRapply(cl, so, function(X) (measr(dat[,c(X,sf)])))
    ind<-c(ind,so[which.min(disA)])
    vlu<-c(vlu,min(disA, na.rm = TRUE))
    sf<-ind

    so<-as.matrix(so1[-sf])
    setTxtProgressBar(bprog,i)}
  stopCluster(cl)

  return(list(CovD=vlu, IdR=sf))}


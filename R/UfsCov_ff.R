#' UfsCov for unsupervised features selection
#'
#' Applies the UfsCov algorithm based on the space filling concept,
#' by using a sequatial forward search (for memory efficient storage
#' of large data on disk and fast access).
#' @usage UfsCov_ff(data, blocks=2)
#' @param data Data of class: \code{matrix} or \code{data.frame}.
#' @param blocks Number of splits to facilitate the computation of the
#' distance matrix (by default: blocks=2).
#' @return A list of two elements:
#'  \itemize{
#'   \item \code{CovD} a vector containing the coverage measure of
#'   each step of the SFS.
#'   \item \code{IdR} a vector containing the added variables during
#'   the selection procedure.
#'   }
#' @author Mohamed Laib \email{Mohamed.Laib@@unil.ch}
#' @note This function is still under developement.
#' @examples
#'
#' \dontrun{
#' #### Infinity dataset ####
#' N <- 1000
#' dat<-Infinity(N)
#' Results<- UfsCov_ff(dat)
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
#' #### Butterfly dataset ####
#'
#' require(IDmining)
#' N <- 1000
#' raw_dat <- Butterfly(N)
#' dat<-raw_dat[,-9]
#'
#' Results<- UfsCov_ff(dat)
#'
#' }
#' @references
#' M. Laib and M. Kanevski (2017). Unsupervised Feature Selection Based on Space
#' Filling Concept, \href{https://arxiv.org/abs/1706.08894}{arXiv:1706.08894}.
#'
#' @import wordspace ff
#' @export

UfsCov_ff<-function(data, blocks=2){
  if (!is.matrix(data) & !is.data.frame(data)) {
    stop('X must be a matrix or a data frame.')
  }
  sbl<-blocks
  so<-so1<-as.matrix(1:ncol(data))
  sf<-0
  ind<-vlu<-c()
  dat<- apply(data , MARGIN = 2,FUN = function(X) (X - min(X))/diff(range(X)))

  for (i in 1:ncol(data)){
    disA<-apply(so, 1,FUN=function(a) (measrff(dat[,c(a,sf)],bl=sbl)))
    #if (all(is.na(disA))==TRUE) {disA<-rep(0,length(disA))}
    ind<-c(ind,so[which.min(disA)])
    vlu<-c(vlu,min(disA, na.rm = TRUE))
    sf<-ind
    so<-as.matrix(so1[-sf])}

  return(list(CovD=vlu, IdR=sf))}



measrff<-function (design,  bl) {
  X1 <- as.ffdf(as.data.frame(design))
  n <- nrow(X1)
  dimension <- ncol(design)
  Distance <- distff(X1, nblocks = bl)
  diag(Distance) <- 1000
  Dmin <- as.matrix( ffrowapply(apply(Distance[1:n,], 2, min), X=Distance, RETURN=TRUE, CFUN="pmin"))
  gammabar <- mean(Dmin)
  s <- sum(apply(Dmin, 2, FUN= function(a)((a-gammabar)*(a-gammabar))))
  cov <- (1/gammabar) * ((1/n) * s)^(1/2)
  return(cov)}



distff <- function(tdata, nblocks=2) {
  nro <- nrow(tdata)
  ffmat <- ff(vmode="single", dim=c(nro,nro))
  splen<-rep_len(1:nblocks, nro)
  splt <- split(1:nro, f=splen)
  COMBS <- expand.grid(1:length(splt), 1:length(splt))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i,]
    if (COMB[1] != COMB[2]) {
      g1 <- splt[[COMB[1]]]
      g2 <- splt[[COMB[2]]]
      slj <-dist.matrix(as.matrix(tdata[c(g1,g2),]), method="euclidean")
      ffmat[c(g1,g2), c(g1,g2)] <- slj

    }
  }
  return(ffmat)}




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
#'
#' @examples
#' Sim_Data<-SimData(n=800)
#' Results<- UfsCov(Sim_Data)
#'
#' cou<-colnames(Sim_Data)
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
#' M. Laib, M. Kanevski, 
#' \href{https://www.elen.ucl.ac.be/Proceedings/esann/esannpdf/es2018-57.pdf}{A novel 
#' filter algorithm for unsupervised feature selection based on a space filling measure}. 
#' Proceedings of the 26rd European Symposium on Artificial Neural Networks, Computational 
#' Intelligence and Machine Learning (ESANN), pp. 485-490, Bruges (Belgium), 2018.
#' 
#' M. Laib and M. Kanevski, A new algorithm for redundancy minimisation in 
#' geo-environmental data, 2019.
#' \href{https://www.sciencedirect.com/science/article/pii/S0098300418310975}{Computers & 
#' Geosciences, 133 104328}.
#'
#' @import Biobase
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
    ind<-c(ind,so[which.min(disA)])
    vlu<-c(vlu,min(disA, na.rm = TRUE))
    sf<-ind

    so<-as.matrix(so1[-sf])
    setTxtProgressBar(bprog,i)

  }
  if(!is.null(colnames(dat))){
    sf <- colnames(data)[sf]
  }

  return(list(CoValue=vlu, VarName=sf))}


measr<-function (design, verbose=F) {
  X <- as.matrix(design)
  n <- nrow(X)
  Dmin <- as.matrix(matchpt(X)[,2]) 
  gammabar <- mean(Dmin)
  s <- sum(apply(Dmin, 1, FUN = function(a) ((a - gammabar)^2)))
  cov <- (1/gammabar) * ((1/n) * s)^(1/2)
  return(cov)}


#' Simulated dataset
#'
#' Generates a simulated dataset (the Infinity dataset)
#' @usage Infinity(n=1000)
#' @param n Number of generated data points (by default: \code{n=1000}).
#' @return A \code{data.frame} of simulated dataset, with \eqn{7} features
#' (\eqn{4} of them are redundants)
#' @author Mohamed Laib \email{Mohamed.Laib@@unil.ch}
#'
#' @examples
#' infinity<-Infinity(n=1000)
#' plot(infinity$x1,infinity$x2)
#'
#' \dontrun{
#'
#' #### Visualisation of the infinity dataset (3D) ####
#' require(rgl)
#' require(colorRamps)
#'
#' c <- cut(infinity$z,breaks=100)
#' cols <- matlab.like(100)[as.numeric(c)]
#' plot3d(infinity$x1,infinity$x2,infinity$z,radius=0.01, col=cols,
#' type="s",xlab="x1",ylab="x2",zlab="z",box=F)
#' grid3d(c("x","y","z"),col="black",lwd=1)
#'
#' }
#' @importFrom stats runif
#' @export

Infinity<-function(n=1000){
  z1 <-runif(n)
  p  <-runif(n, min=0, max=2*pi)
  x1 <-16*sin(p)^3
  x1 <-(x1-min(x1))/diff(range(x1))
  x2 <-sin(p)-sin(2*p)-sin(3*p)
  x2 <-(x2-min(x2))/diff(range(x2))
  z  <- z1*2
  #### redundant features ####
  r1 <-sin(x1+x2)
  r2 <-log(x2+x1)*2
  r3 <-(x1*x2)^2
  r4 <-2*x2+4*x1^2
  return(Infinity=data.frame(x1,x2,z,r1,r2,r3,r4))}

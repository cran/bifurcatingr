#' Bifurcating Autoregressive Tree generator
#'
#' This function generate bifurcating autoregressive (BFA) trees of any size
#' based on a BFA model of any order.
#'
#' @param n tree size (integer)
#' @param p an integer determining the order of bifurcating autoregressive model
#' @param s1 standard deviation of the errors distribution
#' @param s2 standard deviation of the second component of the mixture normal
#'   distribution generating contaminated errors. s2 should be greater than s1.
#'   s2 is only effective when g>0.
#' @param r1 correlation between pairs of errors
#' @param r2 is used in combination with \code{r1} to compute the correlation
#'   between pairs of errors in the second component of the mixture normal
#'   distribution generating the contaminated errors. r2 is only effective when
#'   g>0.
#' @param g proportion of contamination. Defaults to zero producing
#'   non-contaminated multivariate normal errors for the tree generation.
#' @param intercept the intercept in the BAR model generating the tree
#' @param ar.coef a vector of length p giving the autoregressive coefficients in
#'   the BAR model generating the tree
#' @return A numeric vector representing a bifurcating autoregressive (BFA) tree
#'   with \code{n} observations.
#' @export
#' @examples
#' #Non-contaminted BAR(1) tree:
#' bfa.tree.gen(127, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' #Non-contaminted BAR(2) tree:
#' bfa.tree.gen(127, 2, 1, 1, 0.5, 0.5, 0, 10, c(0.5, 0.3))
#' #Contaminted BAR(1) tree:
#' bfa.tree.gen(127, 1, 1, 2, 0.5, 0.5, 0.2, 10, c(0.7))
bfa.tree.gen <- function(n, p, s1, s2, r1, r2, g, intercept, ar.coef){
  if (sum(ar.coef)==1) {stop("Error: the sum of the coefficients must not equal 1")}
  eps <- (g*(s2^2))/(((1-g)*(s1^2)) + g*(s2^2))
  if (g != 0) r=(r1/eps)-r2*((1-eps)/eps)
  if (g == 0) r=0
  #Use start-up BAR(1) tree
  y <- rep(NA,n)
  y[1] <- intercept/(1-sum(ar.coef))
  e <- rcontmnorm((n-1)/2, d=2, sigma1=s1^2*matrix(c(1,r1,r1,1),nrow=2), sigma2=s2^2*matrix(c(1,r,r,1),nrow=2), g=g)
  error <- rep(NA,n)
  for(i in 1:((n-1)/2)){
    error[(i*2)] <- e[i,1]
    error[(i*2)+1] <- e[i,2]
  }
  for(j in 2:n){
    y[j] <- intercept+ar.coef[1]*y[floor(j/(2^1))]+error[j]
  }
  #Generate BAR(p) tree using last subtree with (2^p-1) y-values as starting values
  term.matrix <- matrix(rep(NA,p*n),ncol=p)
  tree <- rep(NA,n)
  tree[1:((2^p)-1)] <- y[bfa.subtree(n,p)]
  e <- rcontmnorm((n-1)/2, d=2, sigma1=s1^2*matrix(c(1,r1,r1,1),nrow=2), sigma2=s2^2*matrix(c(1,r,r,1),nrow=2), g=g)
  error <- rep(NA,n)
  for(i in 1:((n-1)/2)){
    error[(i*2)] <- e[i,1]
    error[(i*2)+1] <- e[i,2]
  }
  for(i in (2^p):n){
    for(k in 1:p){
      term.matrix[,k] <- ar.coef[k]*tree[floor(i/(2^k))]
    }
    tree[i] <- intercept+sum(term.matrix[i,])+error[i]
  }
  return(tree)
}



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
#' @param g proportion of contamination when contaminated normal distribution is
#'    selected, or a positive value representing the degrees
#'    of freedom when skew t-student distribution is selected. Defaults to zero
#'    producing non-contaminated multivariate normal errors.
#' @param intercept the intercept in the BAR model generating the tree
#' @param ar.coef a vector of length p giving the autoregressive coefficients in
#'   the BAR model generating the tree
#' @param dist determine the distribution of the error. Three distributions are
#'  available; Contaminated normal distribution "cnorm", Skew normal distribution
#'  "snorm", and Skew t-student distribution "st".
#' @param a an integer which regulates the the slant of the density when skew
#'    normal distribution or skew t-student distribution is selected. Defaults to zero
#'    producing non-skewed multivariate normal errors, and non-skewed multivariate
#'    t-student errors for the tree generation.
#' @return A numeric vector representing a bifurcating autoregressive (BFA) tree
#'   with \code{n} observations.
#' @importFrom fMultivar rmvsnorm rmvst
#' @export
#' @examples
#' # Non-contaminated normal BAR(1) tree:
#' bfa.tree.gen(127, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' # Non-contaminated normal BAR(2) tree:
#' bfa.tree.gen(127, 2, 1, 1, 0.5, 0.5, 0, 10, c(0.5, 0.3))
#' # Contaminated normal BAR(1) tree:
#' bfa.tree.gen(127, 1, 1, 2, 0.5, 0.5, 0.2, 10, c(0.7))
#' # BAR(1) tree with error generated from skewed-t distribution with skewness equals to -3:
#' bfa.tree.gen(127, 1, 1, 2, 0.5, 0.5, 0, 10, c(0.7),dist="snorm",-3)
#' # BAR(2) tree with error generated from skewed-t distribution with skewness equals to 3:
#' bfa.tree.gen(127, 2, 1, 2, 0.5, 0.5, 0, 10, c(0.7,0.5),dist="snorm",3)
#' # BAR(1) tree with error generated from skewed-t distribution with skewness equals
#' # to -3 and df equals to 10:
#' bfa.tree.gen(127, 1, 1, 2, 0.5, 0.5, 10, 10, c(0.7),dist="st",-3)
#' # BAR(2) tree with error generated from skewed-t distribution with skewness equals
#' # to 3 and df equals to 1:
#' bfa.tree.gen(127, 2, 1, 2, 0.5, 0.5, 10, 10, c(0.7,0.5),dist="st",3)
bfa.tree.gen <- function(n, p, s1, s2, r1, r2, g, intercept, ar.coef, dist="cnorm", a){
  if (sum(ar.coef)==1) {stop("Error: the sum of the coefficients must not equal 1")}
  if (dist=="cnorm"){
  eps <- (g*(s2^2))/(((1-g)*(s1^2)) + g*(s2^2))
  if (g != 0) r=(r1/eps)-r2*((1-eps)/eps)
  if (g == 0) r=0
  #Use start-up BAR(1) tree
  y <- rep(NA,127)
  y[1] <- intercept/(1-sum(ar.coef))
  e <- rcontmnorm((127-1)/2, d=2, sigma1=s1^2*matrix(c(1,r1,r1,1),nrow=2), sigma2=s2^2*matrix(c(1,r,r,1),nrow=2), g=g)
  error <- rep(NA,127)
  for(i in 1:((127-1)/2)){
    error[(i*2)] <- e[i,1]
    error[(i*2)+1] <- e[i,2]
  }
  for(j in 2:127){
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
  if (dist=="snorm"){
      #Use start-up BAR(1) tree
      y <- rep(NA,127)
      y[1] <- intercept/(1-sum(ar.coef))
      e<-fMultivar::rmvsnorm((127-1)/2, dim = 2, mu = rep(0, 2), Omega = matrix(c(s1,r1,r2,s2), 2, 2), alpha = rep(a, 2))
      error <- rep(NA,127)
      for(i in 1:((127-1)/2)){
        error[(i*2)] <- e[i,1]
        error[(i*2)+1] <- e[i,2]
      }
      for(j in 2:127){
        y[j] <- intercept+ar.coef[1]*y[floor(j/(2^1))]+error[j]
      }
      #Generate BAR(p) tree using last subtree with (2^p-1) y-values as starting values
      term.matrix <- matrix(rep(NA,p*n),ncol=p)
      tree <- rep(NA,n)
      tree[1:((2^p)-1)] <- y[bfa.subtree(n,p)]
      e<-fMultivar::rmvsnorm((n-1)/2, dim = 2, mu = rep(0, 2), Omega = matrix(c(s1,r1,r2,s2), 2, 2), alpha = rep(a, 2))
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
  if (dist=="st"){
    #Use start-up BAR(1) tree
    y <- rep(NA,127)
    y[1] <- intercept/(1-sum(ar.coef))
    e <- fMultivar::rmvst(127, dim = 2, mu = rep(0, 2), Omega = matrix(c(s1,r1,r2,s2), 2, 2), alpha = rep(a, 2), df = g)
    error <- rep(NA,n)
    for(i in 1:((127-1)/2)){
      error[(i*2)] <- e[i,1]
      error[(i*2)+1] <- e[i,2]
    }
    for(j in 2:127){
      y[j] <- intercept+ar.coef[1]*y[floor(j/(2^1))]+error[j]
    }
    #Generate BAR(p) tree using last subtree with (2^p-1) y-values as starting values
    term.matrix <- matrix(rep(NA,p*n),ncol=p)
    tree <- rep(NA,n)
    tree[1:((2^p)-1)] <- y[bfa.subtree(n,p)]
    e<-mt <- fMultivar::rmvst(n, dim = 2, mu = rep(0, 2), Omega = matrix(c(s1,r1,r2,s2), 2, 2), alpha = rep(a, 2), df = g)
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
}



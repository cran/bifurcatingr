#' Subtree Extractor
#'
#' This function extracts subtree of size \eqn{(2^p-1)} from the end of a given
#' bifurcating autoregressive tree (model) of order p.
#' @param n subtree size (integer)
#' @param p an integer determining the order of bifurcating autoregressive model
#' @return A numeric vector representing a subtree of size \eqn{(2^p-1)} from
#'   the end of a given bifurcating autoregressive tree.
#' @export
#' @examples
#' bfa.subtree(31, 1)
#' bfa.subtree(31, 2)
bfa.subtree <- function(n, p){
  Rs=floor((n)/(2^(p-1)))
  subtree=c()
  for(R in Rs){
    indx=NULL
    for (i in 0:10){
      ind=((2^i)*R):((2^i)*R+(2^i-1))
      indx=c(indx,if(utils::tail(ind, n=1)<=n){ind})
    }
    subtree <- append(subtree, c(indx))
  }
  return(subtree)
}

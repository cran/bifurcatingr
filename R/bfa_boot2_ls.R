#' Double Bootstrap of Least Squares Estimators of BAR(p) Models
#'
#' This function performs double bootstrapping of the least squares estimators
#' of the autoregressive coefficients in a bifurcating autoregressive (BAR)
#' model of any order \code{p} as described in Elbayoumi and Mostafa (2020).
#'
#' @param z a numeric vector containing the tree data
#' @param p an integer determining the order of bifurcating autoregressive model
#'   to be fit to the data
#' @param burn number of tree generations to discard before starting the
#'   bootstrap sample (replicate)
#' @param B1 number of bootstrap samples (replicates) used in first round of
#'   bootstrapping
#' @param B2 number of bootstrap samples (replicates) used in second round of
#'   bootstrapping
#' @return \item{boot_est}{a matrix containing the first-stage bootstrapped
#' least squares estimates of the autoregressive coefficients} \item{boot2}{a
#' matrix containing the second-stage bootstrapped least squares estimates of
#' the autoregressive coefficients}
#' @references Elbayoumi, T. M. & Mostafa, S. A. (2020). On the estimation bias
#'   in bifurcating autoregressive models. \emph{Stat}, 1-16.
#' @export
#' @examples
#' z <- bfa_tree_gen(31, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' bfa_boot2_ls(z, p=1, B1=99, B2=9)
bfa_boot2_ls <- function(z, p, burn = 5, B1, B2){
  #step 1
  boot1 <- bfa_boot1_ls(z, p, burn = burn, B1, boot_est=TRUE, boot_data=TRUE)
  #step 2
  boot2 <- apply(boot1$boot_data, 1, bfa_boot1_ls, p, burn = burn, B2, boot_est=TRUE)
  list(boot1$boot_est, boot2)
}

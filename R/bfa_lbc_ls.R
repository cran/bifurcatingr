#' Linear Function Bias-Corrected Estimators for BAR(p); p=1,2,...,6
#'
#' This function performs bias correction on the least squares estimators of the
#' autoregressive coefficients in a BAR(p) model based on the assumption that
#' the bias of the least squares estimator is approximately linear as a function
#' of the parameter as described in Elbayoumi and Mostafa (2020).
#' @inheritParams bfa_ls
#' @return \item{coef_lbc}{linear-bias-function-based bias-corrected least
#' squares estimates of the autoregressive coefficients}
#' @references Elbayoumi, T. M. & Mostafa, S. A. (2020). On the estimation bias
#'   in bifurcating autoregressive models. \emph{Stat}, 1-16.
#' @export
#' @examples
#' z <- bfa_tree_gen(127, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' bfa_lbc_ls(z, p=1)
#' z <- bfa_tree_gen(127, 2, 1, 1, 0.5, 0.5, 0, 10, c(0.5, 0.3))
#' bfa_lbc_ls(z, p=2)
bfa_lbc_ls <- function(z, p){
  if (p > 6) {stop("Error: p must be less than 7 for method = 'LBC'.")}
  n <- length(z)
  est <- bfa_ls(z,p)
  coef <- est$coef[1,2:(p+1)]
  error_cor <- est$error_cor
  l <- length(coef)
  coef_lbc <- c(rep(NA,length(coef)))
  if (p == 1)  {
  coef_lbc[1] <- (1/(n-3*(1+error_cor)))*(n*coef[1]+(1+error_cor))
  } else if (p == 2){
    coef_lbc[1] <- (1/(n-(1+error_cor)))*(n*coef[1]+(1+error_cor)*(1+coef[2]))
    coef_lbc[2] <- (1/(n-4*(1+error_cor)))*(n*coef[2]+2*(1+error_cor))
  } else if (p == 3){
    coef_lbc[1] <- (1/(n-(1+error_cor)))*(n*coef[1]+(1+error_cor)*(1+2*coef[3]))
    coef_lbc[2] <- (1/(n-4*(1+error_cor)))*(n*coef[2]+(1+error_cor)*(2+coef[1]+coef[3]))
    coef_lbc[3] <- (1/(n-5*(1+error_cor)))*(n*coef[3]+(1+error_cor))
  } else if (p == 4){
    coef_lbc[1] <- (1/(n-(1+error_cor)))*(n*coef[1]+(1+error_cor)*(1+coef[4]))
    coef_lbc[2] <- (1/(n-(1+error_cor)*(2-coef[1]+coef[3]+2*coef[4])))*(n*coef[2]+2*(1+error_cor))
    coef_lbc[3] <- (1/(n-5*(1+error_cor)))*(n*coef[3]+(1+error_cor)*(1-2*coef[1]+coef[4]))
    coef_lbc[4] <- (1/(n-6*(1+error_cor)))*(n*coef[4]+2*(1+error_cor))
  } else if (p == 5){
    coef_lbc[1] <- (1/(n-(1+error_cor)))*(n*coef[1]+(1+error_cor)*(1+2*coef[5]))
    coef_lbc[2] <- (1/(n-2*(1+error_cor)))*(n*coef[2]+(1+error_cor)*(2-coef[1]+2*coef[4]+coef[5]))
    coef_lbc[3] <- (1/(n-5*(1+error_cor)))*(n*coef[3]+(1+error_cor)*(1-2*coef[1]-coef[2]+coef[4]+2*coef[5]))
    coef_lbc[4] <- (1/(n-6*(1+error_cor)))*(n*coef[4]+(1+error_cor)*(2-coef[1]+coef[5]))
    coef_lbc[5] <- (1/(n-7*(1+error_cor)))*(n*coef[5]+(1+error_cor))
  } else if (p == 6){
    coef_lbc[1] <- (1/(n-(1+error_cor)))*(n*coef[1]+(1+error_cor)*(1+coef[6]))
    coef_lbc[2] <- (1/(n-2*(1+error_cor)))*(n*coef[2]+(1+error_cor)*(2-coef[1]+coef[5]+2*coef[6]))
    coef_lbc[3] <- (1/(n-3*(1+error_cor)))*(n*coef[3]+(1+error_cor)*(1-2*coef[1]-coef[2]+coef[4]+2*coef[5]+coef[6]))
    coef_lbc[4] <- (1/(n-6*(1+error_cor)))*(n*coef[4]+(1+error_cor)*(2-coef[1]-2*coef[2]+coef[5]+2*coef[6]))
    coef_lbc[5] <- (1/(n-7*(1+error_cor)))*(n*coef[5]+(1+error_cor)*(1-2*coef[1]+coef[6]))
    coef_lbc[6] <- (1/(n-8*(1+error_cor)))*(n*coef[6]+2*(1+error_cor))
    }
  coef_lbc
}

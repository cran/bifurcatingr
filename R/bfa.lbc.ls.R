#' Linear Function Bias-Corrected Estimators for BAR(p); p=1,2,...,6
#'
#' This function performs bias correction on the least squares estimators of the
#' autoregressive coefficients in a BAR(p) model based on the assumption that
#' the bias of the least squares estimator is approximately linear as a function
#' of the parameter as described in Elbayoumi and Mostafa (2020).
#' @inheritParams bfa.ls
#' @return \item{coef.lbc}{linear-bias-function-based bias-corrected least
#' squares estimates of the autoregressive coefficients}
#' @references Elbayoumi, T. M. & Mostafa, S. A. (2020). On the estimation bias
#'   in bifurcating autoregressive models. \emph{Stat}, 1-16.
#' @export
#' @examples
#' z <- bfa.tree.gen(127, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' bfa.lbc.ls(z, p=1)
#' z <- bfa.tree.gen(127, 2, 1, 1, 0.5, 0.5, 0, 10, c(0.5, 0.3))
#' bfa.lbc.ls(z, p=2)
bfa.lbc.ls <- function(z, p){
  if (p > 6) {stop("Error: p must be less than 7 for method = 'LBC'.")}
  n <- length(z)
  est <- bfa.ls(z,p)
  coef <- est$coef[1,2:(p+1)]
  error.cor <- est$error.cor
  l <- length(coef)
  coef.lbc <- c(rep(NA,length(coef)))
  if (p == 1)  {
  coef.lbc[1] <- (1/(n-3*(1+error.cor)))*(n*coef[1]+(1+error.cor))
  } else if (p == 2){
    coef.lbc[1] <- (1/(n-(1+error.cor)))*(n*coef[1]+(1+error.cor)*(1+coef[2]))
    coef.lbc[2] <- (1/(n-4*(1+error.cor)))*(n*coef[2]+2*(1+error.cor))
  } else if (p == 3){
    coef.lbc[1] <- (1/(n-(1+error.cor)))*(n*coef[1]+(1+error.cor)*(1+2*coef[3]))
    coef.lbc[2] <- (1/(n-4*(1+error.cor)))*(n*coef[2]+(1+error.cor)*(2+coef[1]+coef[3]))
    coef.lbc[3] <- (1/(n-5*(1+error.cor)))*(n*coef[3]+(1+error.cor))
  } else if (p == 4){
    coef.lbc[1] <- (1/(n-(1+error.cor)))*(n*coef[1]+(1+error.cor)*(1+coef[4]))
    coef.lbc[2] <- (1/(n-(1+error.cor)*(2-coef[1]+coef[3]+2*coef[4])))*(n*coef[2]+2*(1+error.cor))
    coef.lbc[3] <- (1/(n-5*(1+error.cor)))*(n*coef[3]+(1+error.cor)*(1-2*coef[1]+coef[4]))
    coef.lbc[4] <- (1/(n-6*(1+error.cor)))*(n*coef[4]+2*(1+error.cor))
  } else if (p == 5){
    coef.lbc[1] <- (1/(n-(1+error.cor)))*(n*coef[1]+(1+error.cor)*(1+2*coef[5]))
    coef.lbc[2] <- (1/(n-2*(1+error.cor)))*(n*coef[2]+(1+error.cor)*(2-coef[1]+2*coef[4]+coef[5]))
    coef.lbc[3] <- (1/(n-5*(1+error.cor)))*(n*coef[3]+(1+error.cor)*(1-2*coef[1]-coef[2]+coef[4]+2*coef[5]))
    coef.lbc[4] <- (1/(n-6*(1+error.cor)))*(n*coef[4]+(1+error.cor)*(2-coef[1]+coef[5]))
    coef.lbc[5] <- (1/(n-7*(1+error.cor)))*(n*coef[5]+(1+error.cor))
  } else if (p == 6){
    coef.lbc[1] <- (1/(n-(1+error.cor)))*(n*coef[1]+(1+error.cor)*(1+coef[6]))
    coef.lbc[2] <- (1/(n-2*(1+error.cor)))*(n*coef[2]+(1+error.cor)*(2-coef[1]+coef[5]+2*coef[6]))
    coef.lbc[3] <- (1/(n-3*(1+error.cor)))*(n*coef[3]+(1+error.cor)*(1-2*coef[1]-coef[2]+coef[4]+2*coef[5]+coef[6]))
    coef.lbc[4] <- (1/(n-6*(1+error.cor)))*(n*coef[4]+(1+error.cor)*(2-coef[1]-2*coef[2]+coef[5]+2*coef[6]))
    coef.lbc[5] <- (1/(n-7*(1+error.cor)))*(n*coef[5]+(1+error.cor)*(1-2*coef[1]+coef[6]))
    coef.lbc[6] <- (1/(n-8*(1+error.cor)))*(n*coef[6]+2*(1+error.cor))
    }
  return(coef.lbc)
}

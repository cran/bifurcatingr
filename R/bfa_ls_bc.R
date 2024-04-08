#' Bias-Corrected Least Squares Estimators for Bifurcating Autoregressive Models
#'
#' This function performs bias correction on the least squares estimators of the
#' autoregressive coefficients in a BAR(p) model using single, double and
#' fast-double bootstrapping, and the linear-bias-function approach as described
#' in Elbayoumi and Mostafa (2021).
#'
#' @param method method of bias correction. Currently, "boot1", "boot2",
#'   "boot2fast" and "LBC" are supported and they implement single bootstrap,
#'   double bootstrap, fast-double bootstrap, and linear-bias-function
#'   bias-correction, respectively.
#' @inheritParams bfa_boot1_ls
#' @inheritParams bfa_boot2_ls
#' @return \item{coef_ls_bc}{bias-corrected least squares estimates of the
#' autoregressive coefficients}
#' @references Elbayoumi, T. M. & Mostafa, S. A. (2021). On the estimation bias
#'   in bifurcating autoregressive models. \emph{Stat}, e342.
#' @export
#' @examples
#' z <- bfa_tree_gen(63, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' bfa_ls_bc(z, p=2, method="boot1")
#' z <- bfa_tree_gen(63, 2, 1, 1, 0.5, 0.5, 0, 10, c(0.5, 0.3))
#' bfa_ls_bc(z, p=2, method="LBC")
bfa_ls_bc <- function (z, p, method = "boot1", burn = 5, B1 = 999, B2 = 499,
                      boot_est = TRUE, boot_data = FALSE)
{
  if (method == "boot1") {
    coef_ls <- bfa_ls(z, p, error_cor = FALSE)
    coef_ls_bc <- 2 * coef_ls$coef[-1] - colMeans(bfa_boot1_ls(z,
                                                               p, burn = burn, B = B1, boot_est = boot_est)$boot_est)[-1]
    return(coef_ls_bc)
  }
  if (method == "boot2") {
    coef_ls <- bfa_ls(z, p, error_cor = FALSE)
    bt2 <- bfa_boot2_ls(z, p, burn = burn, B1 = B1, B2 = B2)
    bias1 <- colMeans(bt2[[1]])[-1] - coef_ls$coef[-1]
    bt2_1 <- matrix(unlist(lapply(lapply(lapply(bt2[[2]],
                                                FUN = unlist), FUN = matrix, nrow = B2, byrow = FALSE),
                                  FUN = colMeans)), nrow = B1, byrow = TRUE)
    bias2 <- bias1 - colMeans(bt2_1)[-1] + colMeans(bt2[[1]])[-1]
    coef_ls_bc <- coef_ls$coef[-1] - bias1 - bias2
    return(coef_ls_bc)
  }
  if (method == "boot2fast") {
    coef_ls <- bfa_ls(z, p, error_cor = FALSE)
    bt2 <- bfa_boot2fast_ls(z, p, burn = burn, B = B1)
    bias1 <- colMeans(bt2[[1]])[-1] - coef_ls$coef[-1]
    bt2_1 <- matrix(unlist(lapply(bt2[[2]], FUN = unlist)),
                    nrow = B1, byrow = TRUE)
    bias2 <- bias1 - colMeans(bt2_1)[-1] + colMeans(bt2[[1]])[-1]
    coef_ls_bc <- coef_ls$coef[-1] - bias1 - bias2
    return(coef_ls_bc)
  }
  if (method == "LBC") {
    coef_ls_bc <- bfa_lbc_ls(z, p)
  coef_ls_bc
  }
}

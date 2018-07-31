#' @title The Negative Log-Likelehood Function in \code{LVGP} Package
#'
#' @description Calculates the negative log-likelihood (excluding all the constant terms) as described in \code{reference 1}.
#'
#' @param hyperparam Hyperparameters of the LVGP model
#' @param p_quant Number of quantative variables
#' @param p_qual Number of qualitative variables
#' @param lvs_qual Levels of each qualitative variable
#' @param n_lvs_qual Number of levels of each qualitative variable
#' @param dim_z Dimensionality of latent variables, usually 1 or 2
#' @param X_quant Input data of quantative variables
#' @param X_qual Input data of qualitative variables
#' @param Y Vector containing the outputs of data points
#' @param min_eig The smallest eigen value that the correlation matrix is allowed to have, which determines the nugget added to the correlation matrix.
#' @param k Number of data points, \code{nrow(X_quant)} or \code{nrow(X_qual)}
#' @param M Vector of ones with length \code{k}
#'
#' @return The negative log-likelihood (excluding all the constant terms) value.
#'
#' @details \code{\link[LVGP]{LVGP_fit}} calls this function as its optimization objective function.
#'
#' @note This function is \strong{NOT} exported once the package is loaded.
#' @export
#'
#' @references
#' \enumerate{
#' \item "A Latent Variable Approach to Gaussian Process Modeling with Qualitative and Quantitative Factors", Yichi Zhang, Siyu Tao, Wei Chen, and Daniel W. Apley (\href{https://arxiv.org/abs/1806.07504}{arXiv})
#' }
#'
#' @seealso
#' \code{\link[LVGP]{LVGP_fit}} to see how a GP model can be fitted to a training dataset.\cr
#' \code{\link[LVGP]{LVGP_predict}} to use the fitted LVGP model for prediction.\cr
#' \code{\link[LVGP]{LVGP_plot}} to plot the features of the fitted model.
#'
#' @examples
#' # see the examples in the documentation of the function LVGP_fit.

neg_log_l <- function(hyperparam, p_quant, p_qual, lvs_qual, n_lvs_qual, dim_z,
                      X_quant, X_qual, Y, min_eig, k, M){
  if (p_qual == 0) { # no qualitative/categorical variables
    R <- corr_mat(X_quant, X_quant, hyperparam)
  } else {
    z_vec <- hyperparam[-c(1:p_quant)]
    X_qual_la <- to_latent(X_qual, lvs_qual, n_lvs_qual, p_qual, z_vec, dim_z, k)
    X_full <- cbind(X_quant, X_qual_la)
    phi_full <- c(hyperparam[1:p_quant], rep(0, p_qual*dim_z))
    R <- corr_mat(X_full, X_full, phi_full)
  }
  R <- (R + t(R))/2

  raw_min_eig = sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values)[1]
  if (raw_min_eig < min_eig){
    R = R + diag(x = 1, k, k)*(min_eig - raw_min_eig)
  }

  LT <- t(chol(R)) # R = t(L)%*%L = LT%*%t(LT)
  MTLinv <- t(solve(LT, M))
  beta_hat <- MTLinv%*%solve(LT, Y)/sum(MTLinv^2)
  temp <- solve(LT, Y - M%*%beta_hat)
  sigma2 <- sum(temp^2)/k

  det_R <- det(R)
  if ( sigma2 <= 1e-300) {
    sigma2 <- 1e-300
  }
  if (det_R <= 1e-300) {
    det_R <- 1e-300
  }
  n_log_l <- k*log(sigma2) + log(abs(det_R))

  return(n_log_l)
}


#' @title The Function for Constructing the Correlation Matrix in \code{LVGP} Package
#'
#' @description Builds the correlation matrix given two datasets, and the type and parameters of the correlation function.
#'
#' @param X1,X2 Matrices containing the data points. The rows and columns of both \code{X1} and \code{X2} denote individual observation settings and dimension, respectively.
#' @param phi_full The vector storing all the scale (aka roughness) parameters of the correlation function. See \code{reference 1}.
#'
#' @return R The Correlation matrix with size \code{nrow(X1)}-by-\code{nrow(X2)}. See \href{https://en.wikipedia.org/wiki/Correlation_matrix}{here}.
#'
#' @note This function is \strong{NOT} exported once the LVGP package is loaded.
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

corr_mat <-  function(X1, X2, phi_full){

  k = nrow(X1)
  p = ncol(X1)
  kk = nrow(X2)
  R = matrix(0, k, kk)
  phi_full = as.vector(phi_full)

  if (p != length(phi_full)){
    stop(paste('    Length of roughness parameters (', toString(length(phi_full)), ') does not match input dimension (',
               toString(p), ')!'))
  }

  if (k >= kk){
    for (i in 1: kk) {
      R[, i] = colSums((t(X1) - X2[i, ])^2*(10^phi_full))
    }
  } else{
    for (i in 1: k) {
      R[i , ] = colSums((t(X2) - X1[i, ])^2*(10^phi_full))
    }
  }

  R = exp(-0.5*R)

  return(R)
}

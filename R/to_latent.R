#' @title The Function for Transforming Qualitative/Categorical Variables into Latent Variables in \code{LVGP} Package
#'
#' @description Transforms qualitative/categorical variables into latent variables.
#'
#' @param X_qual Matrix or data frame containing (only) the qualitative/categorical data.
#' @param lvs_qual List containing levels of each qualitative variable
#' @param n_lvs_qual Number of levels of each qualitative variable
#' @param p_qual Number of qualitative variables
#' @param z_vec Latent variable parameters, i.e., latent variable values for each level of qualitative/categorical variables
#' @param dim_z Dimensionality of latent variables, usually 1 or 2
#' @param k Number of data points, \code{nrow(X_quant)} or \code{nrow(X_qual)}
#'
#' @return Matrix containing transformed data
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

to_latent <-  function(X_qual, lvs_qual, n_lvs_qual, p_qual, z_vec, dim_z, k){

  X_qual_la <- matrix(0, k, p_qual*dim_z)
  # note: the first levels of each variable are zeros in the latent space, no need to touch upon.
  ind_temp <- as.integer(0)
  for (i in 1:p_qual) {
    n_lvs <- n_lvs_qual[i]
    lvs <- lvs_qual[[i]]
    z_i <- z_vec[(dim_z*ind_temp+1) : (dim_z*(ind_temp+n_lvs-1))]
    ind_temp <- ind_temp+n_lvs-1
    for (j in 2:n_lvs) {
      mask <- X_qual[, i] == lvs[j]
      num_row <- sum(mask)
      if (num_row > 0) {
        X_qual_la[mask, ((i-1)*dim_z+1):(i*dim_z)] <-
          matrix(z_i[(dim_z*(j-2)+1):(dim_z*(j-1))], nrow = num_row, ncol = dim_z, byrow = TRUE)
      }

    }
  }
  return(X_qual_la)
}

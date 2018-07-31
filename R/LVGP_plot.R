#' @title The Plotting Function of \code{LVGP} Package
#'
#' @description Plots the qualitative/categorical variable levels in the latent space (only for 1D or 2D cases).
#'     If the qualitative/categorical variables are not specified, all the qualified variables will be plotted.
#'     See \code{Arguments} for more details on the options.
#'
#' @param model The LVGP model fitted by \code{\link[LVGP]{LVGP_fit}}
#' @param ind_qual_plot An array of index (indices) of the qualitative/categorical variable(s) to be plotted.
#'     Default is NULL, in which case all the qualitative/categorical variables will be plotted.
#'
#' @return NULL
#'
#' @import graphics
#'
#' @note This plot function \strong{only} works for 1D or 2D latent spaces.
#' @export
#'
#' @references
#' \enumerate{
#' \item "A Latent Variable Approach to Gaussian Process Modeling with Qualitative and Quantitative Factors", Yichi Zhang, Siyu Tao, Wei Chen, and Daniel W. Apley (\href{https://arxiv.org/abs/1806.07504}{arXiv})
#' }
#'
#' @seealso
#' \code{\link[LVGP]{LVGP_fit}} to fit LVGP model to the datasets.\cr
#' \code{\link[LVGP]{LVGP_predict}} to use the fitted LVGP model for prediction.
#'
#' @examples
#' # see the examples in the documentation of the function LVGP_fit.

LVGP_plot <- function(model, ind_qual_plot = NULL){
  ind_qual_all <- model$data$ind_qual

  ## checking indices of the variables to be plotted
  if (is.null(ind_qual_all)) {
    stop('    No qualitative/categorical variables in the model. No plot.')
  }
  if (is.null(ind_qual_plot)) {
    ind_qual_plot <- ind_qual_all
    n_plot <- length(ind_qual_plot)
    qual_only_ind <- 1:n_plot
  } else {
    if (!is.numeric(ind_qual_plot)) {
      stop('    The specified index (indices) are invalid!')
    }
    ind_qual_plot <- as.vector(as.integer(ind_qual_plot))
    n_init <- length(ind_qual_plot)
    ind_qual_plot <- unique(ind_qual_plot)
    n_uniq <- length(ind_qual_plot)
    if (n_init != n_uniq) {
      warning('    Ignoring duplicate indices specified.')
    }
    checkers <- (ind_qual_plot %in% ind_qual_all)
    if (!all(checkers)) {
      warning('    Ignoring invalid indices.')
    }
    n_plot <- sum(checkers)
    if (n_plot == 0) {
      stop('    All indices input are invalid! No plot.')
    }
    ind_qual_plot <- ind_qual_plot[checkers]
    qual_only_ind <- match(ind_qual_plot, ind_qual_all)
  }

  ## loading more data
  dim_z <- model$qual_param$dim_z
  z <- model$qual_param$z
  lvs_qual <- model$data$lvs_qual
  n_lvs_qual <- model$data$n_lvs_qual

  ## plotting
  for (i in 1:n_plot) {
    qual_only_ind_i <- qual_only_ind[i]
    z_i <- z[[qual_only_ind_i]]
    lvs_qual_i <- signif(lvs_qual[[qual_only_ind_i]], 3)
    if (dim_z == 1) {
      z_i <- as.vector(z_i)
      stripchart(z_i,
                 main = paste('Latent variables for var #', toString(ind_qual_plot[i]), sep = ''),
                 ylab = "", xlab = "Latent variable value")
      text(z_i, 1.1, labels=lvs_qual_i)
    } else if (dim_z == 2) {
      x_range <- max(z_i[,1]) - min(z_i[,1])
      y_range <- max(z_i[,2]) - min(z_i[,2])
      plot(z_i[,1], z_i[,2], type = "p",
           main = paste('Latent variables for var #', toString(ind_qual_plot[i]), sep = ''),
           ylab = "Latent variable 2 value", xlab = "Latent variable 1 value")
      text(z_i[,1]-x_range*0.02, z_i[,2]-y_range*0.02, labels=lvs_qual_i)
    } else {
      stop('    This plot function only works for 1D or 2D latent spaces!')
    }
  }

}


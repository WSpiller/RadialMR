#' plot_radial
#'
#' A function for producing radial IVW and MR-Egger plots either individually or simultaneously. The function allows for a variety of aesthetic and scaling options, utilising the output from the \code{IVW_radial} and \code{egger_radial} functions.
#'
#' @param r_object An object of class \code{"IVW"} or \code{"egger"}. For visualising both estimates simultaneously, both objects should be included as a vector c(\code{A},\code{B}), where \code{A} and \code{B} denote the \code{"IVW"} and \code{"egger"} objects respectively.
#' @param radial_scale Indicates whether to produce a plot including a full radial scale (\code{TRUE}), or a scatterplot showing only the effect estimates (\code{FALSE}).
#' @param show_outliers Indicates whether display only the set of variants identified as outliers (\code{TRUE}) or the complete set of variants (\code{FALSE}). Note that when (\code{show_outliers=TRUE}), non-outlying variants further from the origin than the furthest outlier will cause an error message that one or more points have been omitted. These are non-outlying variants beyond the scale. If no outliers are present, a plot will be produced using the full set of variants, with an accompanying message indicating the absence of outliers.
#' @param scale_match Indicates whether x and y axes should have the same range(\code{TRUE}), or different ranges (\code{FALSE}) This improves the interpretation of the radial scale, and is set to \code{FALSE} when the radial scale is omitted from the plot.
#' @return A ggplot object containing a radial plot of either the IVW, MR-Egger, or both estimates simultaneously.
#' @author Wes Spiller; Jack Bowden; Tom Palmer.
#' @references Bowden, J., et al., Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression. International Journal of Epidemiology, 2018. 47(4): p. 1264-1278.
#' @importFrom stats quantile median
#' @export
#' @examples
#' ldl.dat <- data_radial[data_radial[,10]<5*10^-8,]
#' ldl.fdat <- format_radial(ldl.dat[,6], ldl.dat[,9],
#'                           ldl.dat[,15], ldl.dat[,21], ldl.dat[,1])
#' ivw.object <- ivw_radial(ldl.fdat, 0.05, 1, 0.0001, TRUE)
#' plot_radial(ivw.object)

plot_radial<-function(r_object,radial_scale,show_outliers,scale_match){


  if(missing(radial_scale)) {
    radial_scale <- TRUE
  }

  if(missing(show_outliers)) {
    show_outliers <- FALSE
  }

  if(missing(scale_match)) {
    scale_match <- TRUE
  }


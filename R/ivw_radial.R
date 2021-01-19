#' ivw_radial
#'
#' Fits a radial inverse variance weighted (IVW) model using a range of weighting specifications. Outliers are determined with respect to their contribution to global heterogeneity, quantified by Cochran's Q-statistic, using a significance threshold specified by the user. The `ivw_radial` function returns an object of class \code{"IVW"}, containing effect estimates, total estimated heterogeneity using Cochran's Q-statistic, the individual contribution to overall heterogeneity of each variant, and a data frame for used in the downstream plotting functions \code{plot_radial} and \code{plotly_radial}.
#'
#' @param r_input A formatted data frame using the \code{format_radial} function.
#' @param alpha A value specifying the statistical significance threshold for identifying outliers (\code{0.05} specifies a p-value threshold of 0.05).
#' @param weights A value specifying the inverse variance weights used to calculate IVW estimate and Cochran's Q statistic. By default modified second order weights are used, but one can choose to select first order (\code{1}), second order (\code{2}) or modified second order weights (\code{3}).
#' @param tol A value indicating the tolerance threshold for performing the iterative IVW approach. The value represents the minimum difference between the coefficients of the previous and current iterations required for a further iteration to be performed (default= \code{0.0001}).
#' @param summary A logical argument (\code{TRUE} or \code{FALSE}) indicating whether a summary of results and heterogeneity should be presented (default= \code{TRUE}).
#' @return An object of class \code{"IVW"} containing the following components:\describe{
#' \item{\code{coef}}{The estimated coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
#' \item{\code{qstatistic}}{Cochran's Q statistic for overall heterogeneity.}
#' \item{\code{df}}{Degrees of freedom. This is equal to the number of variants -1 when fitting the radial IVW model.}
#' \item{\code{outliers}}{A data frame containing variants identified as outliers, with respective Q statistics, chi-squared tests and SNP identification.}
#' \item{\code{data}}{A data frame containing containing SNP IDs, inverse variance weights, the product of the inverse variance weight and ratio estimate for each variant, contribution to overall heterogeneity with corresponding p-value, and a factor indicator showing outlier status.}
#' \item{\code{confint}}{A vector giving lower and upper confidence limits for the radial IVW effect estimate.}
#' \item{\code{it.coef}}{The estimated iterative coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
#' \item{\code{it.confint}}{A vector giving lower and upper confidence limits for the iterative radial IVW effect estimate.}
#' \item{\code{fe.coef}}{The estimated fixed effect exact coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
#' \item{\code{fe.confint}}{A vector giving lower and upper confidence limits for the fixed effect exact radial IVW effect estimate.}
#' \item{\code{re.coef}}{The estimated random effect exact coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
#' \item{\code{re.confint}}{A vector giving lower and upper confidence limits for the random effect exact radial IVW effect estimate.}
#' \item{\code{mf}}{The mean F statistic for the set of genetic variants, indicative of instrument strength.}
#'
#'}
#' @author Wes Spiller; Jack Bowden; Tom Palmer.
#' @references Bowden, J., et al., Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression. International Journal of Epidemiology, 2018. 47(4): p. 1264-1278.

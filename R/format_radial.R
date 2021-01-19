#' format_radial
#'
#' A function which restructures summary GWAS data for downstream two-sample Mendelian randomization analyses. Where variant identification numbers are not provided, an index vector is generated corresponding to the ordering of variants provided.
#'
#' @param BXG A numeric vector of beta-coefficient values for genetic associations with the first variable (exposure).
#' @param BYG A numeric vector of beta-coefficient values for genetic associations with the second variable (outcome).
#' @param seBXG The standard errors corresponding to the beta-coefficients \code{BXG}.
#' @param seBYG The standard errors corresponding to the beta-coefficients \code{BYG}.
#' @param RSID A vector of names for genetic variants included in the analysis. If variant IDs are not provided (\code{RSID = "NULL"}), a vector of ID numbers will be generated.
#'
#' @return The function provides a data frame containing the following columns:\describe{
#'  \item{\code{SNP}}{The identification number for each variant}
#'  \item{\code{beta.exposure}}{The association estimate for the genetic variant with respect to the exposure}
#'  \item{\code{beta.outcome}}{The association estimate for the genetic variant with respect to the outcome}
#'  \item{\code{se.exposure}}{The standard error for the variant-exposure association \code{beta.exposure}}
#'  \item{\code{se.outcome}}{The standard error for the variant-outcome association \code{beta.outcome}}
#'}
#'
#' @author Wes Spiller; Jack Bowden; Tom Palmer.
#' @references Bowden, J., et al., Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression. International Journal of Epidemiology, 2018. 47(4): p. 1264-1278.
#' @export

#' Two Sample Summary GWAS data published in Do et al (2014).
#'
#' A dataset containing summary data GWAS data for 185 independent SNPs with respect to lipid fractions and coronary heart disease (CHD), previously published in Do et al (2013).
#' Lipid fractions include low-density lipoprotein (LDL-C),high-density lipoprotein (HDL-C), and triglycerides, with association estimates
#' obtained from the Global Lipids Genetics Consortium (GLGC) (Willer et al, 2013). Summary data for CHD was obtained from the CARDIoGRAM study (Schunkert et al, 2011).
#' Association estimates were obtained by regressing each phenotype upon the genetic variant and additional adjusted covariates, with links to further information
#' related to each study presented below.
#'
#' data_radial
#'
#' @format A data frame with 185 rows and 21 variables. Specifically this includes the following information:
#' \describe{
#'  \item{\code{rsid}}{The identification number for each variant}
#'  \item{\code{a1}}{The reference allele for each variant}
#'  \item{\code{a2}}{The other allele for each variant}
#'  \item{\code{chr}}{The chromosome number in which each variant is located}
#'  \item{\code{pos}}{The genomic position for the genetic variant relative to chromosome number}
#'  \item{\code{ldlcbeta}}{The association estimate for the genetic variant obtained by regressing LDL-C upon the genetic variant}
#'  \item{\code{hdlcbeta}}{The association estimate obtained by regressing HDL-C upon the genetic variant}
#'  \item{\code{tgbeta}}{The association estimate obtained by regressing triglycerides upon the genetic variant}
#'  \item{\code{chdbeta}}{The association estimate for CHD obtained by regressing CHD upon the genetic variant}
#'  \item{\code{ldlcp2}}{The p-value corresponding to association estimate `ldlcbeta`}
#'  \item{\code{hdlcp2}}{The p-value corresponding to association estimate `hdlcbeta`}
#'  \item{\code{tgp2}}{The p-value corresponding to association estimate `tgbeta`}
#'  \item{\code{chdp2}}{The p-value corresponding to association estimate `chdbeta`}
#'  \item{\code{ldlcz}}{The z-score corresponding to association estimate `ldlcbeta`}
#'  \item{\code{ldlcse}}{The standard error corresponding to association estimate `ldlcbeta`}
#'  \item{\code{hdlcz}}{The z-score corresponding to association estimate `hdlcbeta`}
#'  \item{\code{hdlcse}}{The standard error corresponding to association estimate `hdlcbeta`}
#'  \item{\code{tgz}}{The z-score corresponding to association estimate `tgbeta`}
#'  \item{\code{tgse}}{The standard error corresponding to association estimate `tgbeta`}
#'  \item{\code{chdz}}{The z-score corresponding to association estimate `chdbeta`}
#'  \item{\code{chdse}}{The standard error corresponding to association estimate `chdbeta`}
#' }
#'
#' @source
#' \itemize{
#' \item \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3904346/}
#' \item \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3838666/}
#' \item \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3119261/}
#' }
#' @author Wes Spiller; Jack Bowden; Tom Palmer.
#' @examples
#' head(data_radial)
"data_radial"

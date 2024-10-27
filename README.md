# RadialMR

<!-- badges: start -->
[![R-CMD-check](https://github.com/WSpiller/RadialMR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WSpiller/RadialMR/actions/workflows/R-CMD-check.yaml)
[![r-universe](https://mrcieu.r-universe.dev/badges/RadialMR)](https://mrcieu.r-universe.dev/RadialMR)
<!-- badges: end -->

## Installation

RadialMR can be installed from the [MRCIEU R-Universe](https://mrcieu.r-universe.dev/) with

```r
install.packages("RadialMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
```

To install `RadialMR` directly from the GitHub repository, first make sure you have the `remotes` package installed:

```r
install.packages("remotes")
```

Then the `RadialMR` package can be installed using:

```r
remotes::install_github("WSpiller/RadialMR")
```

To update the package just run the command above again.

## Description

We have written the `RadialMR` R package to produce radial plots and to perform radial
regression for inverse variance weighted and MR-Egger regression models. The package contains a total of five functions:

1. The `format_radial()` function is used to convert a data frame containing summary data into a set format for radial analyses.

2. The `ivw_radial()` function fits a radial inverse variance weighted (IVW) model using either first order, second order, or modified
second order weights. It provides an effect estimate and allows for outliers to be identified using Cochran's Q-statistic. This function 
now also includes iterative and exact IVW estimation, as described in: Improving the accuracy of two-sample summary data Mendelian randomization: 
moving beyond the NOME assumption(https://www.biorxiv.org/content/early/2018/07/02/159442).

3. The `egger_radial()` function fits a radial MR-Egger model using either first order, second order, or modified
second order weights. It provides an effect estimate and allows for outliers to be identified using Rucker's Q-statistic.

4. The `plotly_radial()` function produces interactive radial plots corresponding to the output of the `ivw_radial()` and `egger_radial()` functions.

5. The `plot_radial()` function produces a radial plot corresponding to the output of the `ivw_radial()` and `egger_radial()` functions. The
function provides a range of scaling and aesthetic options showing either an IVW estimate, MR-Egger estimate, or both estimates simultaneously.

Radial plots are produced by many existing R packages such as `metafor`, `numOSL`, and `Luminescence`. Care will need to be taken, however, to input data from an
MR-analysis appropriately into these generic platforms. For this reason we will also continue to develop our own `RadialMR` package to produce radial plots and conduct
radial plot regression for the MR-setting.

## Citation

The paper has been published in the International Journal of Epidemiology:

Bowden, J., et al., Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression. International Journal of Epidemiology, 2018. 47(4): p. 1264--1278. <https://doi.org/10.1093/ije/dyy101>

## License

This project is licensed under GNU GPL v3.

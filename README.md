# RadialMR

## Installation

To install `RadialMR` directly from the GitHub repository, first make sure you have the `devtools` package installed:

    install.packages("devtools")

Then the `RadialMR` package can be installed using:

    library(devtools)
    install_github("WSpiller/RadialMR")
    
To update the package just run the `install_github("WSpiller/RadialMR")` command again.

## Description

We have written the `RadialMR` R package to produce radial plots and to perform radial
regression for inverse variance weighted and MR-Egger regression models. The package contains a total of four functions:

1. The 'format_radial' function is used to convert a data frame containing summary data into a set format for radial analysis.

2. The 'ivw_radial' function fits a radial inverse variance weighted (IVW) model using either first order, second order, or modified
second order weights. It provides an effect estimate and allows for outliers to be identified using Cochran's Q-statistic.

3. The 'egger_radial' function fits a radial MR-Egger model using either first order, second order, or modified
second order weights. It provides an effect estimate and allows for outliers to be identified using Rucker's Q-statistic.

4. The 'plot_radial' function produces a radial plot corresponding to the output of the 'ivw_radial' and 'egger_radial' functions. The
function provides a range of scaling and aesthetic options showing either an IVW estimate, MR-Egger estimate, or both estimates simultaneously. 

Radial plots are produced by many existing R packages such as `metafor`, `numOSL`, and `Luminescence`. Care will need to be taken, however, to input data from an
MR-analysis appropriately into these generic platforms. For this reason we will also continue to develop our own `RadialMR` package to produce radial plots and conduct
radial plot regression for the MR-setting.

## Citation

The preprint is on BioRxiv:

[Improving the visualisation, interpretation and analysis of two-sample summary data Mendelian randomization via the radial plot and radial regression. bioRxiv. doi: https://doi.org/10.1101/200378](https://www.biorxiv.org/content/early/2017/10/11/200378)

## License

This project is licensed under GNU GPL v3.





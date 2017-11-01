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
regression for inverse variance weighted and MR-Egger regression models.

Each plot shows either an IVW estimate, MR-Egger estimate, or both estimates simultaneously. When IVW is plotted alone, a radial curve is included to highlight the ratio estimate for each genetic
variant, as well as the overall IVW estimate. Data points with large contributions to
Cochran’s Q statistic are highlighted in the produced plots, and the significance level for identifying these points as outliers can be defined by the user. 

A table of the exact Q and Q0 contributions for each variant is given as an output for the
researcher to conduct a more detailed analysis. Radial plots are produced by many existing R packages such as `metafor`, `numOSL`, and `Luminescence`. Care will need to be taken, however, to input data from an
MR-analysis appropriately into these generic platforms. For this reason we will also
continue to develop our own `RadialMR` package to produce radial plots and conduct
radial plot regression for the MR-setting.

## Citation

The preprint is on BioRxiv:

[Improving the visualisation, interpretation and analysis of two-sample summary data Mendelian randomization via the radial plot and radial regression. bioRxiv. doi: https://doi.org/10.1101/200378](https://www.biorxiv.org/content/early/2017/10/11/200378)

## License

This project is licensed under GNU GPL v3.





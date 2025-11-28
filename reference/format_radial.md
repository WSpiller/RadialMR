# format_radial

A function which restructures summary GWAS data for downstream
two-sample Mendelian randomization analyses. Where variant
identification numbers are not provided, an index vector is generated
corresponding to the ordering of variants provided.

## Usage

``` r
format_radial(BXG, BYG, seBXG, seBYG, RSID)
```

## Arguments

- BXG:

  A numeric vector of beta-coefficient values for genetic associations
  with the first variable (exposure).

- BYG:

  A numeric vector of beta-coefficient values for genetic associations
  with the second variable (outcome).

- seBXG:

  The standard errors corresponding to the beta-coefficients `BXG`.

- seBYG:

  The standard errors corresponding to the beta-coefficients `BYG`.

- RSID:

  A vector of names for genetic variants included in the analysis. If
  variant IDs are not provided (`RSID = "NULL"`), a vector of ID numbers
  will be generated.

## Value

The function provides a data frame containing the following columns:

- `SNP`:

  The identification number for each variant

- `beta.exposure`:

  The association estimate for the genetic variant with respect to the
  exposure

- `beta.outcome`:

  The association estimate for the genetic variant with respect to the
  outcome

- `se.exposure`:

  The standard error for the variant-exposure association
  `beta.exposure`

- `se.outcome`:

  The standard error for the variant-outcome association `beta.outcome`

## References

Bowden, J., et al., Improving the visualization, interpretation and
analysis of two-sample summary data Mendelian randomization via the
Radial plot and Radial regression. International Journal of
Epidemiology, 2018. 47(4): p. 1264-1278.

## Author

Wes Spiller; Jack Bowden; Tom Palmer.

## Examples

``` r
ldl.dat <- data_radial[data_radial[,10]<5e-8,]
ldl.fdat <- format_radial(ldl.dat[,6], ldl.dat[,9],
                          ldl.dat[,15], ldl.dat[,21], ldl.dat[,1])
head(ldl.fdat)
#>          SNP beta.exposure beta.outcome se.exposure se.outcome
#> 1 rs10903129        -0.033       -0.012 0.003692528 0.01366904
#> 2  rs1998013        -0.380       -0.150 0.021953470 0.09647707
#> 3  rs4587594        -0.049        0.017 0.003842235 0.01509245
#> 4  rs6603981         0.034        0.012 0.004444080 0.01698989
#> 5   rs646776         0.160        0.094 0.004375672 0.01724356
#> 6  rs1010167        -0.025       -0.028 0.003969023 0.01897288
class(ldl.fdat)
#> [1] "data.frame" "rmr_format"
```

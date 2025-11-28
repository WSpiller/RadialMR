# Two Sample Summary GWAS data published in Do et al (2014).

A dataset containing summary data GWAS data for 185 independent SNPs
with respect to lipid fractions and coronary heart disease (CHD),
previously published in Do et al (2013). Lipid fractions include
low-density lipoprotein (LDL-C),high-density lipoprotein (HDL-C), and
triglycerides, with association estimates obtained from the Global
Lipids Genetics Consortium (GLGC) (Willer et al, 2013). Summary data for
CHD was obtained from the CARDIoGRAM study (Schunkert et al, 2011).
Association estimates were obtained by regressing each phenotype upon
the genetic variant and additional adjusted covariates, with links to
further information related to each study presented below.

## Usage

``` r
data_radial
```

## Format

A data frame with 185 rows and 21 variables. Specifically this includes
the following information:

- `rsid`:

  The identification number for each variant

- `a1`:

  The reference allele for each variant

- `a2`:

  The other allele for each variant

- `chr`:

  The chromosome number in which each variant is located

- `pos`:

  The genomic position for the genetic variant relative to chromosome
  number

- `ldlcbeta`:

  The association estimate for the genetic variant obtained by
  regressing LDL-C upon the genetic variant

- `hdlcbeta`:

  The association estimate obtained by regressing HDL-C upon the genetic
  variant

- `tgbeta`:

  The association estimate obtained by regressing triglycerides upon the
  genetic variant

- `chdbeta`:

  The association estimate for CHD obtained by regressing CHD upon the
  genetic variant

- `ldlcp2`:

  The p-value corresponding to association estimate `ldlcbeta`

- `hdlcp2`:

  The p-value corresponding to association estimate `hdlcbeta`

- `tgp2`:

  The p-value corresponding to association estimate `tgbeta`

- `chdp2`:

  The p-value corresponding to association estimate `chdbeta`

- `ldlcz`:

  The z-score corresponding to association estimate `ldlcbeta`

- `ldlcse`:

  The standard error corresponding to association estimate `ldlcbeta`

- `hdlcz`:

  The z-score corresponding to association estimate `hdlcbeta`

- `hdlcse`:

  The standard error corresponding to association estimate `hdlcbeta`

- `tgz`:

  The z-score corresponding to association estimate `tgbeta`

- `tgse`:

  The standard error corresponding to association estimate `tgbeta`

- `chdz`:

  The z-score corresponding to association estimate `chdbeta`

- `chdse`:

  The standard error corresponding to association estimate `chdbeta`

## Source

- <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3904346/>

- <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3838666/>

- <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3119261/>

## Details

data_radial

## Author

Wes Spiller; Jack Bowden; Tom Palmer.

## Examples

``` r
head(data_radial)
#>         rsid a1 a2 chr      pos ldlcbeta hdlcbeta  tgbeta chdbeta  ldlcp2
#> 1 rs10903129  A  G   1 25641524   -0.033  -0.0009 -0.0080 -0.0120 4.0e-19
#> 2  rs4660293  A  G   1 39800767   -0.011   0.0350 -0.0200 -0.0110 1.4e-02
#> 3  rs1998013  T  C   1 55730618   -0.380   0.0350  0.0089 -0.1500 4.0e-67
#> 4 rs10493326  A  G   1 62725961    0.021  -0.0013  0.0310 -0.0087 6.0e-07
#> 5  rs4587594  A  G   1 62906518   -0.049  -0.0150 -0.0690  0.0170 3.0e-37
#> 6  rs6603981  T  C   1 92766395    0.034   0.0039  0.0072  0.0120 2.0e-14
#>    hdlcp2    tgp2 chdp2     ldlcz      ldlcse     hdlcz      hdlcse        tgz
#> 1 7.9e-01 1.7e-02  0.38  8.936966 0.003692528 0.2663106 0.003379512  2.3867077
#> 2 7.0e-19 1.0e-07  0.48  2.457263 0.004476525 8.8748941 0.003943709  5.3267239
#> 3 7.7e-02 6.6e-01  0.12 17.309336 0.021953470 1.7683644 0.019792300  0.4399132
#> 4 7.4e-01 1.0e-15  0.65  4.991217 0.004207391 0.3318533 0.003917393  8.0268589
#> 5 5.0e-05 3.0e-87  0.26 12.752994 0.003842235 4.0556270 0.003698565 19.7996112
#> 6 3.5e-01 7.6e-02  0.48  7.650628 0.004444080 0.9345893 0.004172956  1.7743819
#>          tgse      chdz      chdse
#> 1 0.003351898 0.8778963 0.01366904
#> 2 0.003754653 0.7063026 0.01557406
#> 3 0.020231265 1.5547736 0.09647707
#> 4 0.003862034 0.4537622 0.01917304
#> 5 0.003484917 1.1263911 0.01509245
#> 6 0.004057751 0.7063026 0.01698989
```

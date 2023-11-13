
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TRexDAD

A pipeline for the design of mutagenizing oligonucleotides required for
Tile Region Exchange Mutagenesis (T-Rex), which is a crucial step for
deep mutation scanning (DMS), by harnessing data-optimized assembly
design (DAD).

<!-- badges: start -->

![GitHub](https://img.shields.io/github/license/yunyicheng/TRexDAD)
![GitHub language
count](https://img.shields.io/github/languages/count/yunyicheng/TRexDAD)
![GitHub last commit
(branch)](https://img.shields.io/github/last-commit/yunyicheng/TRexDAD/main)
![GitHub commit activity
(branch)](https://img.shields.io/github/commit-activity/w/yunyicheng/TRexDAD/main)
![GitHub
issues](https://img.shields.io/github/issues/yunyicheng/TRexDAD)
![Discord](https://img.shields.io/discord/1172273306162958441)

<!-- badges: end -->

## Description

Designing oligonucleotides is a significant part of Tile Region Exchange
mutagenesis (T-Rex). However, derive the assembly design by hand is
often error-prone and time-consuming since the fidelity of resulting
oligonucleotides does not merely depend on mere rules. Therefore, a
data-driven solution is required for designing oligonucleotides with
high correct assembly rate. Although there are web tools currently
available for this procedure, few of them utilized Data-optimized
Assembly Design.`TRexDAD` is a streamline tool that facilitates
oligonucleotide design for Tile Region Exchange Mutagenesis by
harnessing Data-optimized Assembly Design (DAD) techniques, along with
visualization of the results. The package is targeted for researchers
who perform Tile Region Exchange Mutagenesis for some genes of interest
and need an efficient and reliable streamline to produce mutagenizing
oligonucleotides. The scope of the R package is to improve the work flow
in Tile Region Exchange Mutagenesis, facilitating open-source researches
for interested bioinformatitians and computational biologists. The
`TRexDAD` package was developed using `R version 4.3.1 (2023-06-16)`,
`Platform: aarch64-apple-darwin20 (64-bit)` and
`Running under: macOS Sonoma 14.1.1`.

## Installation

To install the latest version of the package:

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("yunyicheng/TRexDAD", build_vignettes = TRUE)
library("TestingPackage")
```

<!-- To run the Shiny app: -->
<!-- ``` r -->
<!-- runTestingPackage() # not for Assessment 4; only for Assessment 5 -->
<!-- ``` -->

## Overview

Provide the following commands, customized to your R package. Then
provide a list of user accessible functions within the package and a
brief description of each. Include one image illustrating the overview
of the package that shows the inputs and outputs. Ensure the image is
deposited in the correct location, as discussed in class. Point the user
to vignettes for a tutorial of your package. E.g.,

``` r
ls("package:TRexDAD")
data(package = "TRexDAD") 
browseVignettes("TRexDAD")
```

`TRexDAD` contains 3 functions.

1.  ***InfCriteriaCalculation*** for calculating information criteria
    given dataset dimensions, log-likelihood and probability.

2.  ***NormFactors*** for calculating normalization factors via via
    trimmed mean of M-values (TMM).

3.  ***InfCriteriaPlot*** for plotting information criteria values as a
    scatter plot.

The package also contains two RNA sequencing datasets, called GeneCounts
and GeneCounts2. Refer to package vignettes for more details. An
overview of the package is illustrated below.

![](./inst/extdata/SILVA_A_A1.png)

## Contributions

Provide a paragraph clearly indicating the name of the author of the
package and contributions from the author. Outline contributions from
other packages/sources for each function. Outline contributions from
generative AI tool(s) for each function. Include how the tools were used
and how the results from AI tools were incorporated. Remember your
individual contributions to the package are important. E.g., <br> <br>
<br>

The author of the package is Yunyi Cheng. The author wrote the
*InfCriteriaCalculation* function, which calculates the information
criteria values given data specifications. Here, the Bayesian
information criterion (BIC), Akaike information criterion (AIC) and
Integrated Complete Likelihood (ICL) are calculated. The
*InfCriteriaCalculation* function makes use of map function from
`mclust` R package to generate information criteria values. The `stats`
R package is used for generating multinomially distributed random number
vectors. Part of the code for *InfCriteriaCalculation* function has been
taken from `<NamePackage>` R package. (Section of the borrowed code
should be clearly indicated and referenced in the InfCriteriaCalculation
R script). The *InfCriteriaPlot* is written by the author and generates
a plot of information criteria values. The *InfCriteriaPlot* function
makes use of the `graphics` R package.

*NormFactors* is a function that calculates normalization factors via
Trimmed Mean of M-values (TMM). *NormFactors* function uses Trimmed Mean
of M-values (TMM) as implemented in `edgeR` R package. No generative AI
tools were used in the development of this package.

## References

Provide full references for all sources used, including for the packages
and tools mentioned under ‘Contributions’, in one format. E.g., <br>
<br>

- Akaike, H. (1973). Information theory and an extension of the maximum
  likelihood principle. In *Second International Symposium on
  Information Theory*, New York, USA, 267–281. Springer Verlag.
  <https://link.springer.com/chapter/10.1007/978-1-4612-1694-0_15>.

- Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture
  model for clustering with the integrated classification likelihood.
  *IEEE Transactions on Pattern Analysis and Machine Intelligence* 22.
  <https://hal.inria.fr/inria-00073163/document>

- BioRender. (2020). Image created by Silva, A. Retrieved October 30,
  2020, from <https://app.biorender.com/>

- McCarthy, D. J., Chen Y. and Smyth, G. K. (2012). Differential
  expression analysis of multifactor RNA-Seq experiments with respect to
  biological variation. *Nucleic Acids Research* 40. 4288-4297.
  <https://pubmed.ncbi.nlm.nih.gov/22287627/>

- R Core Team (2023). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna, Austria.
  <https://www.R-project.org/>

- Schwarz, G. (1978). Estimating the dimension of a model. *The Annals
  of Statistics* 6, 461–464.
  <https://projecteuclid.org/euclid.aos/1176344136>.

- Scrucca, L., Fop, M., Murphy, T. B. and Raftery, A. E. (2016) mclust
  5: clustering, classification and density estimation using Gaussian
  finite mixture models. *The R Journal* 8(1), 289-317.
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5096736/>

- Wickham, H. and Bryan, J. (2019). *R Packages* (2nd edition). Newton,
  Massachusetts: O’Reilly Media. <https://r-pkgs.org/>

## Acknowledgements

This package was developed as part of an assessment for 2023 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA.

`TRexDAD` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/yunyicheng/TRexDAD/issues). Many thanks to
those who provided feedback to improve this package.

## Package Structure

The package structure is illustrated below:

![](./inst/extdata/SILVA_A_A2.png) <br> <br> The package tree structure
is provided below.

``` r
- TestingPackage
  |- TestingPackage.Rproj
  |- DESCRIPTION
  |- NAMESPACE
  |- LICENSE
  |- README
  |- data
    |- GeneCounts.rda
    |- GeneCounts2.rda
  |- inst
    CITATION
    |- extdata
      |- SILVA_A_A1.png
      |- GeneCountsData2.csv
    |- shiny-scripts 
        |- app.R
  |- man
    |- GeneCounts.Rd
    |- InfCriteriaCalculation.Rd
    |- NormFactors.Rd
    |- InfCriteriaPlot.Rd
  |- R
    |- data.R
    |- InfCriteriaCalculation.R
    |- InfCriteriaPlot.R
    |- NormFactorCalculation.R
  |- vignettes
    |- TestingPackageVignette.Rmd
  |- tests
    |- testthat.R
    |- testthat
      |- test-InfCriteriaCalculation.R
```

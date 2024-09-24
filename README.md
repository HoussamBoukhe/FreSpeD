# FreSpeD

**Frequency specific change point detection (FreSpeD)** is a method proposed by [Schröder & Ombao, 2019](https://doi.org/10.1080/01621459.2018.1476238) for detecting change points in EEG signals. By locally estimating spectral density matrices and applying CUSUM statistics, FreSpeD efficiently identifies significant changes in the spectral profile. The method has demonstrated speed and effectiveness in detecting epileptic seizures.


## Installation

You can install the development version of <tt>FreSpeD</tt> from
[GitHub](https://github.com) with:

``` r
devtools::install_github("HoussamBoukhe/FreSpeD")
```

## Dependencies
R version: 4.4.1

doParallel 1.0.17
foreach 1.5.2
reshape 1.4.4
ggplot2 3.5.1

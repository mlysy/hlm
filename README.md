# hlm: Inference for a Heteroscedastic Linear Model

*Martin Lysy, Tian You, Yifan Wang*

---

### Description

Provides efficient algorithms to compute the MLE of a heteroscedastic regression model, with linear mean and log-linear variance covariates.  An Expectation-Maximization algorithm is provided for right-censored responses, which are expected in the context of survival analysis for the corresponding heteroscedastic accelerated failure-time model of [Wang et al (2019)](https://arxiv.org/abs/1508.05137).

### Installation

Install the R package [**devtools**](https://CRAN.R-project.org/package=devtools) and run
```r
devtools::install_github("mlysy/hlm")
```

### Usage

The main function in the package is `hlm()`, which you would typically call like this:
```r
M <- hlm(y ~ x1 + x2 | x1)
```
This would fit a heteroscedastic model with mean linear in `x1` and `x2`, and variance log-linear in `x1`.  For more examples please see `?hlm`.

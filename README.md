# Approximate Bayesian Inference for Cumulative Probit Regression Models

This repository accompanies the manuscript [Aliverti (2025), *Approximate Bayesian Inference for Cumulative Probit Regression Models*](http://arxiv.org/abs/2511.06967), and implements three scalable variational algorithms for Bayesian inference under the cumulative probit model.

The main package, `viord`, can be installed using:

```
devtools::install_github("emanuelealiverti/EPCP", subdir = "viord")
```

The core function, `viord()`, performs approximate inference for the cumulative probit model using **Expectation Propagation (EP)** or **Variational Bayes** with either **Mean-Field** or **Partially-Factorized Mean-Field** assumptions; all routines are implemented in `c++` and available using such unified interface. The package also performs thresholds optimization via (approximate) maximum-marginal likelihood. 

The package provides a unified interface and standard methods, including: `summary()`, `predict()` (via predictive probabilities) and `simulate()` (to sample from the approximate posterior distribution).

The `tutorials` folder includes illustrative examples that show the use of the paper in practice. In particular, [`tutorials/brazilian.md`](https://github.com/emanuelealiverti/EPCP/blob/main/tutorials/brazilian.md) reproduces the analysis of the **Brazilian Bank dataset** (Section 4.1 of the paper).

# Quantile Regression using *I*-spline and Neural Network

QUINN is a novel Bayesian framework for simultaneous estimation of non-crossing, non-linear quantile curves. It models the conditional CDF non-parametrically using *I*-spline basis expansion, where the spline coefficients depend on the covariates through a neural network. 
We assign prior to weight parameters, so that uncertainty of the model can be analyzed. Simulation study shows that QUINN, compared to existing models, can recover the underlying quantile process more accurately, especially when the sample size is small. Although QUINN is
considered a "black-box" model, some interpretability can be obtained by estimating its quantile marginal effects using model-agnostic tools. We extended the recently proposed accumulative local effects plot (**ALEplot**) [[1]](#1) to quantile regression and estimated main and second-
order interaction effects of QUINN. We applied QUINN to analyze tropical cyclone intensity and U.S. birthweight data.

## Usage 

- R folder

The **quinn.R** contains the main functions to fit a QUINN model. The `quinn_samp` function draws posterior samples of the weight parameters using no-U-turn sampling (NUTS); the specific implementation is largely based on the **adnuts** package [[2]](#2). The `quinn_pred`
function estimates the quantile process using the posterior samples.

The **utils.R** contains helper functions fot implementing the NUTS.

The **ALEPlot_qt.R** can be used to estimate the quantile marginal effects of QUINN

The **example_1d.R** and **example_2d.R** provide two examples, one univariate and one bivariate, of using QUINN. 

- Data folder

Example datasets to which QUINN can be applied. The folder currenly contains **hurricane.csv**, which is the tropical cyclone intensity data analyzed by Elsner et. al [[3]](#3).

## References

<a id="1">[1]</a> 
Apley, D.W. and Zhu, J., 2020. Visualizing the effects of predictor variables in black box supervised learning models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 82(4), pp.1059-1086.

<a id="2">[2]</a> 
Monnahan, C.C. and Kristensen, K., 2018. No-U-turn sampling for fast Bayesian inference in ADMB and TMB: Introducing the adnuts and tmbstan R packages. PloS one, 13(5), p.e0197954.

<a id="3">[3]</a> 
Elsner, J.B., Kossin, J.P. and Jagger, T.H., 2008. The increasing intensity of the strongest tropical cyclones. Nature, 455(7209), pp.92-95.
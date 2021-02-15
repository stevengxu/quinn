# Quantile Regression using *I*-spline and Neural Network (QUINN)

QUINN is a novel Bayesian framework for simultaneous estimation of non-crossing, non-linear quantile curves. It models the conditional CDF non-parametrically using *I*-spline basis expansion, where the spline coefficients depend on the covariates through a neural network. 
We assign prior to weight parameters, so that uncertainty of the model can be analyzed. Simulation study shows that QUINN, compared to existing models, can recover the underlying quantile process more accurately, especially when the sample size is small. Although QUINN is
considered a "black-box" model, some interpretability can be obtained by estimating its quantile marginal effects using model-agnostic tools. We extended the recently proposed accumulative local effects plot (ALEplot) to quantile regression and estimated main and second-
order interaction effects of QUINN. We applied QUINN to analyze tropical cyclone intensity and U.S. birthweight data.

##Usage 

- R folder

The **quinn** file contains the main functions to fit a QUINN model. The `quinn_samp` function draws posterior samples of the weight parameters using no-U-turn sampling (NUTS); the specific implementation was largely based on the **adnuts** package (Monnahan and Kristensen, 2018). The `quinn_pred`
function estimates the quantile process using the posterior samples.

The **utils** file contains helper functions fot implementing the NUTS.

The **example_1d** and **example_2d** files provide two examples, one univariate and one bivariate, of using QUINN. 

- data foler

Example datasets to which QUINN can be applied. The folder currenly contains `hurricane.csv`, which is the tropical cyclone intensity data analyzed by Elsner et. al (2008).

## References
<a id="1">[1]</a> 
Monnahan, C.C. and Kristensen, K., 2018. No-U-turn sampling for fast Bayesian inference in ADMB and TMB: Introducing the adnuts and tmbstan R packages. PloS one, 13(5), p.e0197954.

<a id="2">[2]</a> 
Elsner, J.B., Kossin, J.P. and Jagger, T.H., 2008. The increasing intensity of the strongest tropical cyclones. Nature, 455(7209), pp.92-95.
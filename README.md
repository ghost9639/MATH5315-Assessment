
# Table of Contents

1.  [Project Specifications](#orgd5575ce)
2.  [Pearson Family Fitting (Battery Lifespans)](#orgb7f013d)
3.  [ARIMA and VAR (Stock Close Prices)](#org9bd1778)
4.  [ARIMA (COVID Mortality)](#org70bab1f)
5.  [GARCH Modelling (Stock Returns)](#org7b36a89)



<a id="orgd5575ce"></a>

# Project Specifications

This is a project exploring a wide variety of datasets, from financial to industry production, and evaluating them using standard time series methods. All datasets have been included in data/, and the code is linked and explained in [a markdown file](src/Main.rmd). A full [report](Project.pdf) has been compiled as well.


<a id="orgb7f013d"></a>

# Pearson Family Fitting (Battery Lifespans)

The first section predicts the likelihood of a battery failing within a certain period using the Pearson family of distributions. After examining the Kolmogorov-Smirnov test, AIC, and BIC values, it was concluded that the exponential distribution was the most accurate, but plots are made for both.

<images/babelHardDriveLifeFit.pdf>


<a id="org9bd1778"></a>

# ARIMA and VAR (Stock Close Prices)

ARIMA is applied to stock close price datasets in order to estimate near future stock prices. After concluding that serial correlation was too much of a problem for estimating the individual stocks, the choice was made to fit an Engle-Granger model instead, with greater explanatory power for a set of cointegrated stocks and more datapoints to use to avoid serial correlation. Predictions were then generated and visualised for the portfolio, and no evidence of serial correlation or inappropriate residuals was discovered in post regression testing.

<images/babelEuroVAR.pdf>


<a id="org70bab1f"></a>

# ARIMA (COVID Mortality)

ARIMA and log ARIMA models were fitted and compared with an extreme focus on model fit and testing. Keen attention paid to evidence of continued non-stationarity through ACF plots and augmented Dickey-Fuller testing. Projections were generated and compared for near future mortality due to COVID.

<images/babelDeathPlot.pdf>


<a id="org7b36a89"></a>

# GARCH Modelling (Stock Returns)

Conversion of stock prices to stock returns was followed for ARCH testing. After non-random volatility was detected and an aDF test confirmed non-unit variance, a GARCH model was fitted. Post regression testing examined a suite of tests to ensure parameter selection was appropriate, and model alpha and beta carefully evaluated and explained in simple terms.


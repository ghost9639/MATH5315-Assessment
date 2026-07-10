
# Table of Contents

1.  [Project Specifications](#org439d352)
2.  [Pearson Family Fitting (Battery Lifespans)](#org4414cb1)
3.  [ARIMA and VAR (Stock Close Prices)](#org95a4102)
4.  [ARIMA (COVID Mortality)](#org8f8083d)
5.  [GARCH Modelling (Stock Returns)](#org3a55989)



<a id="org439d352"></a>

# Project Specifications

This is a project exploring a wide variety of datasets, from financial to industry production, and evaluating them using standard time series methods. All datasets have been included in data/, and the code is linked and explained in [a markdown file](src/Main.rmd). A full [report](Project.pdf) has been compiled as well.


<a id="org4414cb1"></a>

# Pearson Family Fitting (Battery Lifespans)

The first section predicts the likelihood of a battery failing within a certain period using the Pearson family of distributions. After examining the Kolmogorov-Smirnov test, AIC, and BIC values, it was concluded that the exponential distribution was the most accurate, but plots are made for both.

![img](images/babelHardDriveLifeFit.png "Hard Drive Lifetimes with Fitted Gamma and Exponential distributions")


<a id="org95a4102"></a>

# ARIMA and VAR (Stock Close Prices)

ARIMA is applied to stock close price datasets in order to estimate near future stock prices. After concluding that serial correlation was too much of a problem for estimating the individual stocks, the choice was made to fit an Engle-Granger model instead, with greater explanatory power for a set of cointegrated stocks and more datapoints to use to avoid serial correlation. Predictions were then generated and visualised for the portfolio, and no evidence of serial correlation or inappropriate residuals was discovered in post regression testing.

![img](images/babelEuroVAR.png "VAR Stock Price Projections")


<a id="org8f8083d"></a>

# ARIMA (COVID Mortality)

ARIMA and log ARIMA models were fitted and compared with an extreme focus on model fit and testing. Keen attention paid to evidence of continued non-stationarity through ACF plots and augmented Dickey-Fuller testing. Projections were generated and compared for near future mortality due to COVID.

![img](images/babelDeathPlot.png "Autocorrelation Tests")


<a id="org3a55989"></a>

# GARCH Modelling (Stock Returns)

Conversion of stock prices to stock returns was followed for ARCH testing. After non-random volatility was detected and an aDF test confirmed non-unit variance, a GARCH model was fitted. Post regression testing examined a suite of tests to ensure parameter selection was appropriate, and model alpha and beta carefully evaluated and explained in simple terms.


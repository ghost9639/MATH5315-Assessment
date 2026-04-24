install.packages("tidyverse") # plots and variable handling
install.packages("here")      # reproducible file searching
install.packages("aTSA")      # stationarity hypothesis testing
install.packages("forecast")  # ARIMA modelling
install.packages("vars")      # VAR model
install.packages("fGarch")    # GARCH model
install.packages("FinTS")     # ArchTest function


here::here() # should echo back your project root

suppressMessages(library(tidyverse))
suppressMessages(library(here))
hard_drives <- read.table(here("data", "hard_drive.txt"))

hard_drives <- tibble(hard_drives[1])

hard_drives$V1 <- as.numeric(hard_drives$V1)

hard_drives
summary(hard_drives)
sd(hard_drives$V1)

suppressMessages(library(patchwork))
plt <- ggplot(data = hard_drives, aes(x=V1))
plt <- plt + geom_histogram(binwidth = 0.5)
plt <- plt + labs(title = "Lifetimes of batteries", x = "Years")

plt2 <- ggplot (data = hard_drives, aes(x=V1))
plt2 <- plt2 + geom_boxplot()
plt2 <- plt2 + labs(title = "Boxplot of battery lives", x = "Years")

(plt / plt2) + plot_layout(heights = c(4, 5))

suppressMessages(require(MASS))

exp_mod <- fitdistr(hard_drives$V1, "exponential")
gam_mod <- fitdistr(hard_drives$V1, "gamma")

exp_rate <- exp_mod$estimate['rate']
gam_shape <- gam_mod$estimate['shape']
gam_rate <- gam_mod$estimate['rate']
exp_lk <- exp_mod$loglik
gam_lk <- gam_mod$loglik

print("MASS stores key model fit information that we can access and read.")
sprintf("The exponential distribution rate is %.3f, and the gamma
distribution follows shape %.3f and rate %.3f", exp_rate,
gam_shape, gam_rate)
sprintf("The log-likelihood for the exponential distribution is %f, and
the log-likelihood for the gamma distribution is %f.",
exp_lk, {gam_mod$loglik})


hist(x=hard_drives$V1, freq = FALSE, breaks = 10,
     main = "Histogram of hard drive lifetimes with fitted distributions",
     xlab = "Hard Drive Lifetimes (years)")

curve(dexp(x, rate = exp_mod$estimate["rate"]),
      add = TRUE, col = "red", lwd = 2)

curve(dgamma(x, rate = gam_mod$estimate["rate"],
             shape = gam_mod$estimate["shape"]),
      add = TRUE, col = "blue", lwd= 2)

ks.test(hard_drives$V1, "pexp", exp_rate)
ks.test(hard_drives$V1, "pgamma", shape = gam_shape, rate = gam_rate)

## calculating AIC BIC values for exp and gam

n <- dim(hard_drives) 
AIC_exp <- 2 * 1 - 2 * exp_mod$loglik
AIC_gam <- 2 * 2 - 2 * gam_mod$loglik

BIC_exp <- 1 * log(n) - 2 * exp_mod$loglik
BIC_gam <- 2 * log(n) - 2 * gam_mod$loglik

fit_stats <- data.frame (
    Criterion = c("AIC", "BIC"),
    `Exponential Model` = c(AIC_exp, BIC_exp[1]),
    `Gamma Model` = c(AIC_gam, BIC_gam[1])
)

fit_stats

set.seed(1234)
exp_sim <- rexp(n = 200, rate = exp_rate)
gam_sim <- rgamma(n = 200, shape = gam_shape, rate = gam_rate)
sims <- data.frame(exp_sim, gam_sim)

## gamma histogram
plt3 <- ggplot(data = sims, aes(x=gam_sim))
plt3 <- plt3 + geom_histogram(binwidth = 0.5)
plt3 <- plt3 + labs(title = "Gamma simulated lifetimes", x = "Years")

## gamma boxplot
plt4 <- ggplot (data = sims, aes(x=gam_sim))
plt4 <- plt4 + geom_boxplot()
plt4 <- plt4 + labs(title = "Gamma simulated summary", x = "Years")

## exponential histogram
plt5 <- ggplot(data = sims, aes(x=exp_sim))
plt5 <- plt5 + geom_histogram(binwidth = 0.5)
plt5 <- plt5 + labs(title = "Exponential simulated lifetimes", x = "Years")

## gamma boxplot
plt6 <- ggplot (data = sims, aes(x=exp_sim))
plt6 <- plt6 + geom_boxplot()
plt6 <- plt6 + labs(title = "Exponential simulated summary", x = "Years")

## plt  plt2
## plt3 plt4
## plt5 plt6

(plt | plt2) / (plt3 | plt4) / (plt5 | plt6)


dclose <- tibble(read.csv(here("data", "EuroIndices_Close.csv")))

## lubridate supplies dmy parsing
dclose$Date <- strptime(dclose$Date, "%b %d, %y")

plt <- ggplot(data = dclose, aes(x = as.POSIXct(Date)))
plt <- plt + geom_line(aes(y = AEX_Close, colour = "AEX")) +
    geom_line(aes(y = BEL20_Close, colour = "BEL20"))+
    geom_line(aes(y = CAC40_Close, colour = "CAC40")) +
    geom_line(aes(y = DAX_Close, colour = "DAX")) +
    labs(title = "European Indices",
         y = "Closing Price (€)",
         x = "Time")
plt


par(mfrow = c(2,1))

acf(dclose$DAX_Close, main = "ACFs for DAX Close")
pacf(dclose$DAX_Close, main = "Partial ACFs for DAX Close")

par(mfrow = c(1,1))

suppressMessages(library(aTSA))
## augmented dickey-fuller can test stationarity in the presence
## of serial-correlation

adf.test(x = dclose$DAX_Close, nlag = 2)

d1dclose <- dclose[-1, ]
d1dclose[,-1] <- lapply(dclose[,-1], diff)
d1dclose

adf.test(x = d1dclose$DAX_Close, nlag = 2)


suppressMessages(library(forecast))

print("First we consider an AR(p) model.")
mod_AR_DAX <- arima(dclose$DAX_Close, c(1,0,0))

print("Second, we consider an MA(q) model.")
mod_MA_DAX <- arima(dclose$DAX_Close, c(0,0,1))

print("Finally, we select an ARIMA model.")
mod_arima_DAX <- auto.arima(dclose$DAX_Close)

AIC_AR <- mod_AR_DAX$aic
AIC_MA <- mod_MA_DAX$aic
AIC_ARIMA <- mod_arima_DAX$aic

sprintf("The AR model has an AIC of %.2f, the MA model has an AIC
of %.2f, and the ARIMA model has an AIC of %.2f", AIC_AR, AIC_MA,
AIC_ARIMA)

proj <- forecast(mod_arima_DAX, h = 10)

autoplot(proj, ylab = "DAX20 Closing Price")

adf_checks <- (lapply(d1dclose[,-1], function(x)
    adf.test(x, nlag = 2)$p.value))

coint.test(d1dclose$DAX_Close, d1dclose$AEX_Close, nlag = 2)
coint.test(d1dclose$DAX_Close, d1dclose$BEL20_Close, nlag = 2)
coint.test(d1dclose$DAX_Close, d1dclose$CAC40_Close, nlag = 2)
coint.test(d1dclose$AEX_Close, d1dclose$BEL20_Close, nlag = 2)
coint.test(d1dclose$AEX_Close, d1dclose$CAC40_Close, nlag = 2)
coint.test(d1dclose$BEL20_Close, d1dclose$CAC40_Close, nlag = 2)

cointegration <- data.frame(Indices = c("AEX", "BEL20", "CAC40", "DAX"),
                        AEX = c(NA, 0.01, 0.01, 0.01),
                        BEL20 = c(0.01, NA, 0.01, 0.01),
                        CAC40 = c(0.01, 0.01, NA, 0.01),
                        DAX = c(0.01, 0.01, 0.01, NA))

cointegration

suppressMessages(library(vars))

mod_var <- VARselect(d1dclose[,-1], lag.max = 5, type = c("const"))
mod_var

temp_mod <- VAR(d1dclose[,-1], p = 2)
var_proj <- predict(temp_mod, h = 10)
par(mar = c(4, 4, 2, 1))
plot(var_proj)
par(mfrow = c(1,1))

Q2_arima_mod <- arima(dclose$DAX_Close, c(0,1,1))
Box.test(Q2_arima_mod$residuals, lag = 2, type = c("Ljung-Box"))
shapiro.test(Q2_arima_mod$residuals)

ts.diag(arima(dclose$DAX_Close, c(0,1,1)))

VAR_temp_mod <- VAR(d1dclose[,-1], p = 2)

## serial correlation
serial.test(VAR_temp_mod)

## normality tests
normality.test(VAR_temp_mod)

deaths <- tibble(read.csv(here("data", "deaths.csv")))

deaths$date <- strptime(deaths$date, "%Y-%m-%d")

par(mfrow = c(3,1))

ts.plot(deaths$wk.deaths, main = "Deaths over time",
        ylab = "Deaths")
acf(deaths$wk.deaths, main = "ACFs for deaths")
pacf(deaths$wk.deaths, main = "PACFs for deaths")

par(mfrow = c(1,1))

adf.test(deaths$wk.deaths, nlag = 3)

acf(deaths$wk.deaths)

d.deaths <- deaths[-1,]
d.deaths$wk.deaths <- diff(deaths$wk.deaths)

adf.test(d.deaths$wk.deaths, nlag = 3)


par(mfrow = c(3,1))

ts.plot(d.deaths$wk.deaths, main = "Change in deaths over time",
        ylab = "Deaths")
acf(d.deaths$wk.deaths, main = " ")
pacf(d.deaths$wk.deaths, main = " ")

par(mfrow = c(1,1))

death_mod <- auto.arima(deaths$wk.deaths, d = 1)
death_mod # fixed d = 1 or default auto.arima picks a worse model

proj_deaths <- forecast::forecast(death_mod, h = 2)
autoplot(proj_deaths, main = "Death projections for next two weeks",
         ylab = "Deaths")

sqdeaths <- deaths
sqdeaths$wk.deaths <- sqrt(deaths$wk.deaths)

mod_sqd <- auto.arima(sqdeaths$wk.deaths, d = 1)

sqd_proj <- forecast::forecast(mod_sqd, h = 2)

sqd_proj$x <- sqd_proj$x^2
sqd_proj$fitted <- sqd_proj$fitted^2
sqd_proj$lower <- sqd_proj$lower^2
sqd_proj$upper <- sqd_proj$upper^2
sqd_proj$mean <- sqd_proj$mean^2
sqd_proj$residuals <- sqd_proj$residuals^2


autoplot(sqd_proj,
         main = "Back-fitted projections from rooted deaths",
         ylab = "Deaths")

temp_death_mod <- arima(deaths$wk.deaths, c(3, 1, 0))

ts.diag(temp_death_mod)

temp_mod_sqd <- arima(sqdeaths$wk.deaths, c(2,1,2))

ts.diag(temp_mod_sqd)

Box.test(death_mod$residuals, lag = 2, type = c("Ljung-Box"))
Box.test(mod_sqd$residuals, lag = 2, type = c("Ljung-Box"))


shapiro.test(death_mod$residuals)
shapiro.test(mod_sqd$residuals)

ex_rat <- tibble(read.csv(here("data", "USD-GBP.csv")))

ex_rat$Date <- strptime(ex_rat$Date, "%d/%m/%Y")

ex_returns <- ex_rat[-1,]
ex_returns$USD.GBP.Close <- diff(log(ex_rat$USD.GBP.Close))

ts.plot(ex_returns$USD.GBP.Close)

adf.test(ex_returns$USD.GBP.Close, nlag = 2)

## We can use dplyr to filter by criteria
## I think tidyverse masks filter but we call it explicitly
pre_crash <- ex_returns %>%
    dplyr::filter(Date < as.POSIXct("2007-08-01"))

post_crash <- ex_returns %>%
    dplyr::filter(Date >= as.POSIXct("2007-08-01"))

pre_dur <- range(pre_crash$Date)
post_dur <- range(post_crash$Date)

sprintf("The pre-crash data ranges from %s to %s, and the
post-crash data ranges from %s to %s.",
as.character(pre_dur[1,]), as.character(pre_dur[2]),
as.character(post_dur[1,]), as.character(post_dur[2]))

auto.arima(pre_crash$USD.GBP.Close)
auto.arima(post_crash$USD.GBP.Close)

arima_mod_pre <- arima(pre_crash$USD.GBP.Close, c(0,0,1))
arima_mod_post <- arima(post_crash$USD.GBP.Close, c(0,0,1)) 

shapiro.test(arima_mod_pre$residuals)
shapiro.test(arima_mod_post$residuals)

Box.test(arima_mod_pre$residuals, lag = 2, type = c("Ljung-Box"))
Box.test(arima_mod_post$residuals, lag = 2, type = c("Ljung-Box"))

ts.diag(arima_mod_pre) # all residuals look fine

ts.diag(arima_mod_post)

suppressMessages(library(FinTS))

ArchTest(pre_crash$USD.GBP.Close)
ArchTest(post_crash$USD.GBP.Close)

suppressMessages(library(fGarch))

## iterating over models for pre-crisis
grm_pre <- garchFit(formula = USD.GBP.Close ~ arma(1,0) + garch(1,1),
                data = pre_crash,
                trace = FALSE) # -6.692718 -6.525458, NaNs produced

grm_pre <- garchFit(formula = USD.GBP.Close ~ arma(0,0) + garch(1,1),
                data = pre_crash,
                trace = FALSE) # -5.956245 -5.822436, NaNs produced

grm_pre <- garchFit(formula = USD.GBP.Close ~  arma(1,2) + garch(1,1),
                data = pre_crash,
                trace = FALSE) # -7.667488 -7.433323


## iterating over models for post-crisis
grm_post <- garchFit(formula = USD.GBP.Close ~  + garch(1,1),
                data = post_crash,
                trace = FALSE) # -1.665750 -1.657883

grm_post <- garchFit(formula = USD.GBP.Close ~ arma(2,1)+ garch(1,1),
                data = post_crash,
                trace = FALSE) # -6.434903 -6.421137


summary(grm_pre)
summary(grm_post)


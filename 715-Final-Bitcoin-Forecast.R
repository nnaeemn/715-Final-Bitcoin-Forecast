
setwd("~/Desktop/bitcoin")
library("tidyverse")
library("flipTime")
library("forecast")
library("tseries")
library("dplyr")
library("ggplot2")
library("scales")
library("lubridate")
library("xts")
library("dplyr")


# Create the new data from the new file.
price0 <- read.csv("~/Desktop/bitcoin/coindesk-bpi-USD.csv")
price0$Date <- AsDate(price0$Date)
price <- price0[-1,]

par(mfrow=c(2,2))
ggplot(data = price, aes(x=Date, y=Close.Price))+ geom_line(color="blue") + 
  ggtitle("Plot 1. Bitcoin Daily Market Price: January 2014 - May 2018") + 
  scale_x_date(labels = date_format("%d-%m-%Y"), breaks = date_breaks("months")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  theme(plot.title=element_text(size=12, vjust=5, hjust = 0.4)) + 
  labs(y="Market Price (USD)") + theme(plot.margin = unit(c(1.5,1.5,0.75,0.75), "cm"))



# Plot the series

ggplot(data = price, aes(x=Date, y=log(Close.Price)))+ geom_line(color="blue") + 
  ggtitle("Plot 2. Log-Transformed Bitcoin Daily Market Price: January 2014 - April 2018") + 
  scale_x_date(labels = date_format("%d-%m-%Y"), breaks = date_breaks("months")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  theme(plot.title=element_text(size=12, vjust=5, hjust = 0.4)) + 
  labs(y="log(Market Price (USD))") + theme(plot.margin = unit(c(1.5,1.5,0.75,0.75), "cm"))

# We log-transform the series to contain the instability due to its volatile and explosive 
# growth behavior, and look at the ACF and PACFs of both the transformed and original series. 

par(mfrow=c(2,2), oma = c(0,0,2,0))
acf(price$Close.Price, lag.max = 60, main = "ACF - Original Series")
pacf(price$Close.Price, lag.max = 60, main = "PACF - Original Series")
acf(log(price$Close.Price), lag.max = 60, main = "ACF - Log-transformed Series")
pacf(log(price$Close.Price), lag.max = 60, main = "PACF - Log-transformed Series")
title("Plot 3", outer = TRUE)


# Both of the ACF plots show extreme positive correlations at many lags, decreasing linearly and very slowly. 
# However, this might be the strong effect at lag 1 that is propagating back to the previous lags. The PACF 
# plot confirms this; it shows a significant spike (equal almost to 1) at lag 1 and then suddenly cuts off. 
# These add more evidence against stationarity. We draw a final conclusion by checking the ADF test for stationarity.     


xts.series <- xts(price$Close.Price, order.by = price$Date)
time.series<-msts(price$Close.Price, seasonal.periods = c(7,30,365.25), start=c(2014))
time.series.log <- log(time.series)

adf.test(time.series)

# The tests null hypothesis of non-stationarity cannot be rejected since the p-value is large assuming a significance level 
# of 0.05, confirming our assumption of non-stationarity.  

# We also look at the decomposition plot of the log-transformed series. The results in the following plot are based on a 
# frequency of 365.25(to count for leap years). We see the general upward trend after the first year and some pattern is 
# discernable in the residuals. 


price.decomp <- stl(time.series.log, s.window = "periodic")
plot(price.decomp, main = "Plot 4. Log-transformed Series Decomposition")
deseasonal_price <- seasadj(price.decomp)
price$deseasoned <- exp(deseasonal_price)
plot(price$deseasoned)
# We proceed with disregarding seasonality for now, as it requires further analysis.

# We now stationarizing the log-transformed series by taking differences. Before that we need to tentatively determine the 
# order of differencing. The ACF plot above showed large correlations at many lags and PACF had a high spike at lag 1. 
# Thus we proceed by differencing once and look at some plots to determine if further differencing is needed.   

fit.arima000 <- Arima(log(price$Close.Price), order = c(0,1,0))
tsdisplay(residuals(fit.arima000), lag.max=60, main='Plot 5. 
Residuals, ACF and PACF of the differenced log-transformed series')

# The plot of the residuals has the appearance of a stationary process, with a mean that appears to be constant around 0. 
# However, the series exhibits significant volatility clustering as larger fluctuations cluster together, particularly in 
# the begining and in the later parts of the series. This indicates that the conditional variance is changing, implying 
# that we may not have independence and identical distribution. We can check this with the ACF and PACF for the absolute 
# and squared differences, since if the differences are iid, then so are the absolute and squared ones.   

par(mfrow=c(2,2), oma=c(0,0,2,0))
acf(abs(diff(log(price$Close.Price))), main = "ACF of the  abs(diff(log(series)))")
pacf(abs(diff(log(price$Close.Price))), main = "PACF of the abs(diff(log(series)))") 
acf(diff(log(price$Close.Price))^2, main = "ACF of the diff(log(series))^2")
pacf(diff(log(price$Close.Price))^2, main = "PACF of the diff(log(series))^2")
title("Plot 6", outer = T)


# We see that there are significant correlations for multiple lags in all four plots, suggesting that the differences are not 
# independent and identically distributed. But we may still have uncorrelatedness, which is necessary for a series to be white
# noise. The ACF and PACF plots of the differenced series do not show anything concerning as there are
# no statistically significant correlations at any lag. Our empirical conclusion based on these plots is that the differenced 
# log-transformed series is stationary. We check the ADF test as well, and it confirms our conclusion. These suggest that a 
# white-noise model works. We may now proceed to finding appropriate models for our series.   

adf.test(diff(time.series))


# We need to specify the orders of AR and MA for our ARIMA model by looking at the ACF and PACF plots. However, neither of them
# really give us a clue. They neither have any significant lags, particularly at lag 1 they are both very close to 0, nor do 
# they show any strong patterns, such as cutting off at a particular lag. Although the very small value at lag 1 is negative 
# in both, suggesting an MA term for our model. Since there is no conclusive judgement, we will try a few combinations manually
# and by using the auto.arima function. 

fit.auto.arima <- auto.arima(price$Close.Price, lambda = 0)
tsdisplay(residuals(fit.auto.arima), lag.max=60, main='Auto Arima(1,2,0) Residuals')
summary(fit.auto.arima)

# The auto.arima function suggests an ARIMA(1,2,0) model as the best model, contrary to both of our specifications of first order
# differencing and addition of an MA term. The AR coefficient given is -0.4791 with a standard error of 0.0220 and a sigma^2 of 
# about 0.0025. The RMSE is quite high at 303.1. We check the residulas, ACF and PAC plots.  


fit.auto.arima <- auto.arima(price$Close.Price, lambda = 0)
tsdisplay(residuals(fit.auto.arima), lag.max=60, main='Plot 7. ARIMA(1,2,0)')


# The residuals plot exhibits some volatility cluster, overal the mean seem to be constant around zero, with possible one or two
# outliers. The ACF plot of the residuals show significant negative autocorrelation at lags 1 and 2, and the PACF plot shows negative
# spikesat the first several lags, cutting off only after about lag 20. These suggests that the model does not fit the data very well
# and we can model the residuals. Below we look at the results and plots for the model with our specifications, i.e., ARIMA(0,1,1).   

fit.arima000 <- Arima(price$deseasoned, lambda=0, order = c(0,1,1))
summary(fit.arima000)
tsdisplay(residuals(fit.arima000), lag.max=60, main='Plot 8. ARIMA(0,1,1)')


# There are no big differences between the two residual plots, maybe one more outlier in the later model. However, the ACF and PACF plots 
# look much better, with almost no significant lags at all. We will compare the forecast results from the two models for a 11-day 
# period(since after trying different periods this one gave the best result).   

library("forecast")
fit09 <- Arima(time.series[1:1589], order=c(1,2,0), seasonal = F, lambda = 0)
tsdisplay(residuals(fit09), lag.max=60, main='Residuals')
fcast0 <- forecast(fit09, h=11)
plot(fcast0, main = "Plot 9. Forecasts From ARIMA(1,2,0)")
de <- price$Close.Price
lines(de)

train.error0 <- sqrt(mean((price[1:1589,]$Close.Price - as.matrix(as.numeric(fcast0$fitted)))^2))


preds <- data.frame(price[1590:1600,]$Close.Price, as.numeric(fcast0$mean), (price[1590:1600,]$Close.Price - as.numeric(fcast0$mean)))

# Get the test error for the ARIMA(1,2,0) model.
test.error0 <- sqrt(mean((price[1590:1600,]$Close.Price - as.matrix(as.numeric(fcast0$mean)))^2))



fit01 <- Arima(time.series[1:1589], order=c(0,1,1), seasonal = F, lambda = 0)
tsdisplay(residuals(fit09), lag.max=60, main='Residuals')
fcast01 <- forecast(fit01, h=11)
plot(fcast01, main="Plot 10. Forecasts From ARIMA(0,1,1)")
de <- price$Close.Price
lines(de)

# Get the test error for the ARIMA(0,1,1) model.
train.error1 <- sqrt(mean((price[1:1589,]$Close.Price - as.matrix(as.numeric(fcast01$fitted)))^2))
test.error1 <- sqrt(mean((price[1590:1600,]$Close.Price - as.matrix(as.numeric(fcast01$mean)))^2))


# The model suggested by auto.arima performs better as the test error is smaller. Although it is 
# still quite high at 588.321, it  predicts the direction and overal slop correctly. However, the confidence interval in
# this case is much larger than the (0,1,1) model. The forecast from the ARIMA(0,1,1) is almost simply a straight line, 
# that is the mean. The training errors are smaller at about 303 for the ARIMA(1,2,0) model and 244 for the ARIMA(0,1,1) 
# model. 

# We lastly plot the fitted series by our best model against the original series, including the forecasted part. 

par(mfrow=c(1,1))
plot(price$Close.Price , col = "blue", type = "l", lwd=1, ylim=c(0,22000), main= "Plot 11. ARIMA(1,2,0) fitted series against original
fitted - orange curve
original - blue curve", ylab="Daily Price")
lines(fitted(fcast0), col = "orange", lwd = 0.6)
lines(fcast0$mean, col = "red", lwd = 1.5)


# We next want to determine if there are ARCH effects in the series. We look at the ACF and PACF plots of the squared and 
# absolute returns (differences) back in plot 6. All of the ACF and PACF plots exhibit significant correlations at multiple
# lags. Thus, the residuals do not seem to be independent and so we can try to fit a GARCH model in order to model the 
# volatility clustering of the series. We confirm this using the McLeod-Li test below. We see that the p-values for all
# lags are significant, suggesting presence of ARCH effects. The qq-norm plot also indicates non-normality.  


library("TSA")
r <- diff(log(price$Close.Price)*100)
par(mfrow=c(1,2))
McLeod.Li.test(y=r, main = "Plot13. Mc. Leod-Li Test")
qqnorm(r, main = "Plot 14. Normal Q-Q Plot") ; qqline(r)


# We would now like to determine the orders of our GARCH model. We check the EACF table for the squared differences and it is
# not very clear, it could be suggesting a (1,2) or a (2,2) GARCH model. The same table for the absolute differences has a more
# clear cut shape and suggests a (1,2) model. We will try both.  

eacf(r^2)
eacf(abs(r))

library(fGarch)
model = garchFit(formula = ~ garch(1, 2), data = diff(time.series.log)[1:1590], cond.dist = "norm", include.mean = TRUE, trace = FALSE)
fcst=predict(model,n.ahead=10)
tsdisplay(residuals(model), lag.max = 60, main = "Plot 15. GARCH(1,2) Residuals, ACF and PACF")
plot(fcst)  
mean.fcst=fcst$meanForecast
sqrt(mean((diff(time.series.log)[1590:1599]-mean.fcst)^2))
qqnorm(residuals(model), main="Plot 16. Normal-QQ Plot GARCH(1,2) Residuals") ; qqline(residuals(model))
summary(model)
gBox(model1, method = "square")

model1 = garchFit(formula = ~ garch(2, 2), data = diff(time.series.log)[1:1589], cond.dist = "norm", include.mean = TRUE)
tsdisplay(residuals(model1), lag.max = 60, main="Plot17. GARCH(2,2) Model")
fcst1=predict(model1,n.ahead=11)

plot(fcst)  
mean.fcst=fcst$meanForecast
qqnorm(residuals(model)) ; qqline(residuals(model))
summary(model)
gBox(model, method = "square")


# We see that the residuals plot of the GARCH(1,2) model does not look very different from our ARIMA models. However, the ACF and
# PACF plots af the residuals are better, without significant autocorrelations at any lag. p-values in the summary are small as 
# well. Issues are encountered when using the GARCH model and in performing the gBox test (perhaps because 
# I used a different package to fit the GARCH model, since the garch function from the Tseries package did not work properly), and
# so I am not including them here. The plots we went over suggested that a white noise model is appropriate for this series, 
# specially the ACF and PACF of the absolute and squared differences. The Mc.Leod-Li test also suggested the presence of an 
# ARCH effect in our time series.


















################
#HongKong analysis
################

#Import data
# read it as a ts object
# analyse stationarity
# Transform the data to reach stationarity
# identify the orders of the SARIMA model
HKtraf=read.csv2("Hongkong-1998-2016.csv")
str(HKtraf)
traf=ts(HKtraf[,2],start=c(1998,1), frequency=12)

plot.ts(traf, xlab = "Dates (monthly)", ylab = "Traffic", main = "Traffic at HongKong Airport")
abline(reg=lm(traf~time(traf)), col = "red")
plot(decompose(traf))

monthplot(traf)

#The series is not stationary: trend, seasonal effects
# Idea: split the sample around 2008, January to avoid the change of trend

#subsample selection
#sub.traf=window(traf,start=c(1998,1), end=c(2008,8))
#plot(decompose(sub.traf))
#View(sub.traf)

#Stationarity for the sample

par(mfrow=c(2,2))
plot(traf, main="")
acf(ts(traf,frequency = 1), main='ACF')
pacf(ts(traf,frequency = 1), main='PACF')

# The series is non stationary: we observe a trend and seasonal effects
# The ACF has stochastic trend setting


# Turning it into log
ltraf= log(traf)
par(mfrow=c(2,2))
plot(ltraf, main="")
acf(ts(ltraf,frequency = 1), main = "ACF")
pacf(ts(ltraf,frequency = 1), main = "PACF")


#1st diff of log data (to remove the trend)
dltraf=diff(ltraf,1)
par(mfrow=c(2,2))
plot(dltraf)
acf(ts(dltraf,frequency = 1), main = "ACF")
pacf(ts(dltraf,frequency = 1), main = "PACF")

# Remove seasonality
par(mfrow=c(2,2))
dltraf12=diff(dltraf,12)
plot(dltraf12)
acf(ts(dltraf12,frequency = 1), main = "ACF")
pacf(ts(dltraf12,frequency = 1), main = "PACF")

#-----------------------
# Estimation of a first model
#----------------------------
#q=1
#p=2
#Q=1
#P=1

mod1 = arima(ltraf, c(2,1,1), seasonal = list(order=c(1,1,1), period=12), method = "ML")
mod1
# AIC=-706.88

library('TSA')
fit1=fitted(mod1)
par(mfrow=c(1,1))
plot(cbind(ltraf,fit1), plot.type = "single", col=c("black","red"))

# Validation
#-----------

#Significance of coefficients
# Residuals analysis: gaussian white noise
# Prediction accuracy


mod1
mod1$coef #value of coeff
mod1$var.coef #variance of coeff
tstat = mod1$coef / sqrt(diag(mod1$var.coef)) # Student test statistic
tstat
pvalue = 2*(1-pnorm(abs(tstat))) #pvalue
pvalue

#         ar1         ar2         ma1        sar1        sma1 
# 0.041439865 0.005267540 0.230601043 0.304787131 0.004721735 


# We verify a 2nd model with p=2, P=0, q=1, Q=1 because the coefficient of P in the previous model
# was not significant, even though the q value is not relevant we will discard 1 by 1

mod2 = arima(ltraf, c(2,1,1), seasonal = list(order=c(0,1,1), period=12), method = "ML")
mod2
tstat = mod2$coef / sqrt(diag(mod2$var.coef)) # Student test statistic
tstat
pvalue = 2*(1-pnorm(abs(tstat))) #pvalue
pvalue
#        ar1         ar2         ma1        sma1 
#0.017920714 0.001735534 0.230381894 0.000000000

# We verify a 3rd model with p=2, P=0, q=0, Q=1 because the coefficient of q in the previous model
# was not significant

mod3 = arima(ltraf, c(2,1,0), seasonal = list(order=c(0,1,1), period=12), method = "ML")
mod3
tstat = mod3$coef / sqrt(diag(mod3$var.coef)) # Student test statistic
tstat
pvalue = 2*(1-pnorm(abs(tstat))) #pvalue
pvalue

#         ar1          ar2         sma1 
#0.000000e+00 2.078052e-09 3.552714e-15  
# pvalues for AR1, AR2 and SMA1 are smaller than 5%

mod3
# AIC -708.32
# model3 is better since AIC is smaller


# Residuals analysis
#-------------------
library("forecast")
checkresiduals(mod3)

# The ACf of the residuals has no significant coefficient

# The residuals behave like a white noise

# The histogram has bell shape close to a normal distribution


# Standardized residuals
res_norm=mod3$residuals/sqrt(mod3$sigma2)

plot(res_norm)
abline(a=-2, b=0, col="red")
abline(a=2, b=0, col="red")
# we have four points outside the interval, so we will use 4 dummy variables to get rid of them







qqnorm(mod3$residuals)
qqline(mod3$residuals)
# The distribution of residuals is very close to a Gaussian distribution
# Prediction accuracy
#--------------------
# Confidence bounds around the fitted value


library("TSA")
fit3 = fitted(mod3)
cb80 = mod3$sigma2^.5*qnorm(0.9)
plot(cbind(ltraf,fit3-cb80,fit1+cb80), 
     plot.type="single", lty=c(1,2,2))

#proportion of points in the confidence interval
indi = (ltraf-(fit3-cb80))>0&(fit3+cb80-ltraf)>0
prop = 100*sum(indi)/length(indi)
prop

#In-sample / out-of sample analysis
#----------------------------------


#Idea: split the sample into 2 subsamples: training set and test set

data.train=window(ltraf,start=c(1998,1), end=c(2014,1))
str(data.train) #193 observations

data.test = window(ltraf, start=c(2014,2), end=c(2016,12))
str(data.test) #35 obs
mod3.train = arima(data.train, c(2,1,0), seasonal=list(order=c(0,1,1), period=12), method="ML")

pred3.test = predict(mod3.train, n.ahead=31)
accuracy(pred3.test$pred, data.test)

# On average, we observe 38% of error with respect to the true value (MAPE result)
# The best model would be the model with the lowest value for all these parameters

plot(traf, xlim=c(2012,2016),ylim=c(300,800))
lines(2.718^(pred3.test$pred), col="red")
lines(2.718^(pred3.test$pred-1.96*pred3.test$se), col=4,lty=2)
lines(2.718^(pred3.test$pred+1.96*pred3.test$se), col=4,lty=2)


#outlier identification
out1= which(mod3$residuals < -0.1) #identification of outlier
out1
#the outlier corresponds to the observation 50 in the dataset

install.packages("zoo")
library("zoo")
index(mod3$residuals)[out1] #date of the outlier
#View(HKtraf)


#Adding a external variable X=fitting a SARIMAX model
HKtraf$dum1=0
HKtraf$dum1[out1]=1

mod4 =arima(ltraf, c(2,1,0), seasonal=list(order=c(0,1,1), period=12), method="ML",
            xreg = HKtraf$dum1)
mod4
#AIC=-786.92
res_norm=mod4$residuals/sqrt(mod4$sigma2)

plot(res_norm)
abline(a=-2, b=0, col="red")
abline(a=2, b=0, col="red")


library("forecast")
checkresiduals(mod4)
mod4.train = arima(
  data.train,
  c(2, 1, 0),
  seasonal = list(order = c(0, 1, 1), period = 12),
  method = "ML"
  ,
  xreg = HKtraf$dum1[1:193]
)

pred4.test = predict(mod4.train, n.ahead=30, newxreg = 0)
accuracy(pred4.test$pred, data.test)
#MAPE=34% instead of 38%, the predictive power of the model has improved!

plot(traf, xlim=c(2012,2018),ylim=c(300,800))
lines(2.718^(pred4.test$pred), col="red")
lines(2.718^(pred4.test$pred-1.96*pred4.test$se), col=4,lty=2)
lines(2.718^(pred4.test$pred+1.96*pred4.test$se), col=4,lty=2)


# Prediction for future dates (one year ahead)
pred = predict(mod3, n.ahead=12)
2.718^pred$pred
lines(2.718^pred$pred, col="blue")

plot(pred$pred)


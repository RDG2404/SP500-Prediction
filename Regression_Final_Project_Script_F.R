library(quantmod)
library(tidyquant)
library(tidyverse)
library(fredr)
library(bizdays)
library(ggplot2)
library(Hmisc)
library(corrplot)
library(dplyr)
library(ggpubr)
library(bestglm)
library(leaps)
library(rio)
library(glmnet)
library(car)
library(MASS)

fredr_set_key("19bf118e8e53e9ffe49461f2b4c9648e")
create.calendar("Brazil/ANBIMA", holidaysANBIMA, weekdays=c("saturday", "sunday"))

start_date = "2015-01-01"
end_date = "2019-12-31"

#SP500 DATA
final_data = getSymbols("^GSPC", from = start_date, to = end_date, warnings = FALSE, auto.assign =FALSE)
final_data = data.frame(final_data)[4]
date = as.Date(row.names(final_data))
final_data = cbind(date, final_data)
row.names(final_data) = 1:nrow(final_data)
colnames(final_data) = c("date","SP500")


#Collecting Data from Fred
fred_symbols = c("DGS10", "DTB6", "DGS30", "BAMLH0A0HYM2", "BAMLC0A3CA", 
                 "T10YIE", "EFFR", "T5YIE")
for(i in 1:length(fred_symbols)){
  x = fredr(fred_symbols[i], observation_start = as.Date(start_date), 
            observation_end = as.Date(end_date))
  x = x[c("date", "value")]
  colnames(x) = c("date", fred_symbols[i])
  final_data = merge(final_data, x, by = "date")
  rm(x)
}

#Collecting Data from Yahoo Finance
yfin_symbols = c("USDEUR=X", "USDGBP=X", "^FVX")
for(i in 1:length(yfin_symbols)){
  data = getSymbols(yfin_symbols[i], from = start_date, to = end_date, warnings = FALSE, auto.assign =FALSE)
  data = data.frame(data)[4]
  Date = as.Date(row.names(data))
  data = cbind(Date, data)
  row.names(data) = 1:nrow(data)
  colnames(data) = c("date",yfin_symbols[i])
  final_data = merge(final_data, data, by = "date")
  rm(data)
}
final_data$date <- as.Date(final_data$date)


for(i in 2:ncol(final_data)){
  final_data[i] = ((final_data[i]/lag(final_data[i]))-1)*100
}
final_data$date <- floor_date(final_data$date, "monthly")
final_data <- aggregate(final_data[, 2:13], by = list(final_data$date), FUN = mean)
names(final_data)[1] <- "date"

#Getting Monthly Data
monthly_symbols = c("INDPRO", "CPIAUCSL", "UNRATE", "FEDFUNDS", "UMCSENT",
                    "DSPIC96", "M1SL", "CURRCIR", "M2SL")
for(i in 1:length(monthly_symbols)){
  x = fredr(monthly_symbols[i], observation_start = as.Date(start_date), 
                       observation_end = as.Date(end_date))
  x = x[c("date", "value")]
  colnames(x) = c("date", monthly_symbols[i])
  final_data = merge(final_data, x, by = "date")
  rm(x)
}

for(i in 14:ncol(final_data)){
  final_data[i] = ((final_data[i]/lag(final_data[i]))-1)*100
}

final_data = na.locf(final_data, na.rm = FALSE)
final_data = na.locf(final_data, na.rm = FALSE, fromLast = TRUE)
final_data = final_data[, 2:22]

initial_model = lm(SP500~., data = final_data)
summary(initial_model)

#Because most independent variables are macro-variables, we wanted to check for multicollinearity
print("We want the VIF values to be below")
max(10, (1/(1-0.7659)))

vif(initial_model)

#We can clearly see that the model shows multicollinearity, thus we different methods of variable selection and PCA Analysis

#PCA Analysis
final_data.pca <- prcomp(final_data, center = TRUE, scale = TRUE)

# screeplot
screeplot(final_data.pca, type = "l", npcs = 10, main = "Screeplot of First 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue=1"), col=c("red"), lty=5, cex=0.6)

# cumulative variance
cumpro <- cumsum(final_data.pca$sdev^2/sum(final_data.pca$sdev^2))
plot(cumpro, 
     xlab = "Principal component (PC)", ylab = "Variance explained", main = "Cumulative Variance Plot",
     yaxp = c(0, 1, 10))
abline(v = 10, col="blue", lty=5)
abline(h = cumpro[10], col="blue", lty=5)
principal_components = data.frame(final_data.pca$x[,1:10])
principal_components["SP500"] = final_data$SP500

pca_model = lm(SP500~., data = principal_components)
summary(pca_model)

##The following variable selections are done on reduced variables after accounting for VIF elimination.
#Forward Selection
attach(final_data)

step(lm(SP500~DTB6, data = final_data), scope=list(lower=SP500 ~ DTB6, upper = 
                                                      SP500 ~ DTB6 + BAMLH0A0HYM2 + 
                                                      BAMLC0A3CA + EFFR + `USDEUR=X` + 
                                                      `USDGBP=X` + INDPRO + CPIAUCSL + 
                                                      UNRATE + FEDFUNDS + UMCSENT + DSPIC96 +
                                                      M1SL + CURRCIR + M2SL), direction = "forward")
forward_model = lm(formula = SP500 ~ DTB6 + BAMLH0A0HYM2 + BAMLC0A3CA + M2SL, 
                   data = final_data)
summary(forward_model)

#Backward Selection
full = lm(SP500 ~ DTB6 + BAMLH0A0HYM2 + BAMLC0A3CA + EFFR + `USDEUR=X` + 
            `USDGBP=X` + INDPRO + CPIAUCSL + UNRATE + FEDFUNDS + UMCSENT + 
            DSPIC96 + M1SL + CURRCIR + M2SL) 
minimum = lm(SP500~DTB6)
step(full, scope=list(lower=minimum, upper=full), direction="backward")

backward_model = lm(formula = SP500 ~ DTB6 + BAMLH0A0HYM2 + EFFR + INDPRO + CURRCIR + 
                      M2SL)
summary(backward_model)


## We compare Ridge, Lasso, Elastic Net on Selected Variables.
#Lasso
Xpred = cbind(DTB6, BAMLH0A0HYM2, BAMLC0A3CA, EFFR, `USDEUR=X`, `USDGBP=X`, INDPRO,
              CPIAUCSL, UNRATE, FEDFUNDS, UMCSENT, DSPIC96, M1SL, CURRCIR, M2SL)
ymodel.cv = cv.glmnet(Xpred, SP500, alpha=1, nfolds=10)
ymodel = glmnet(Xpred, SP500, alpha=1, nlambda=100)

coef(ymodel, s=ymodel.cv$lambda.min)


#Ridge
attach(final_data)
predictors = scale(Xpred)
y.scaled = scale(SP500)

## Apply ridge regression for a range of penalty constants
lambda = seq(0, 10, by=0.25)
out = lm.ridge(y.scaled~predictors, lambda=lambda)
round(out$GCV, 4)
which(out$GCV == min(out$GCV))
round(out$coef[,10], 4)

#Elastic Net
ymodel.elastic = glmnet(Xpred, SP500, alpha=0.5, nlambda = 100)

## Extract coefficients at optimal lambda
coef(ymodel.elastic, s=ymodel.cv$lambda.min)

# Mallow Cp
out=leaps(Xpred, SP500, method="Cp")
cbind(as.matrix(out$which), out$Cp)
best.model=which(out$Cp==min(out$Cp))
cbind(as.matrix(out$which), out$Cp)[best.model,]


###########Predicted vs Actual for all regressions
model.lasso=lm(SP500~DTB6+BAMLH0A0HYM2+BAMLC0A3CA)
model.elasticnet=lm(SP500~DTB6+BAMLC0A3CA+BAMLH0A0HYM2+EFFR+`USDGBP=X`+INDPRO+UNRATE+FEDFUNDS+CURRCIR+M2SL)
model.ridge=lm(SP500~ 
                 DTB6 + BAMLH0A0HYM2 + 
                 BAMLC0A3CA + EFFR + `USDEUR=X` + 
                 `USDGBP=X` + INDPRO + CPIAUCSL + 
                 UNRATE + FEDFUNDS + UMCSENT + DSPIC96 +
                 M1SL + CURRCIR + M2SL)
model.mcp=lm(SP500~DTB6+BAMLH0A0HYM2+M2SL)
summary(model.lasso)
summary(model.elasticnet)
summary(model.ridge)

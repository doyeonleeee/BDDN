library(Rcpp)
library(RcppArmadillo)
library(coda)
library(dplyr)

 
## load data ##
XdatS = read.csv("data/XdatS.csv", stringsAsFactors = FALSE)
Y = as.matrix(read.csv("data/Y.csv", row.names = 1, check.names = FALSE))

## use centered square root of X as a predictor ##
X2=sqrt(XdatS[,-1])
m=mean(X2$hours)
X2$hours=X2$hours-m


model=model.frame(Y~X2[,1]+X2[,2]+X2[,1]:X2[,2]+I(X2[,2]^2)+X2[,1]:I(X2[,2]^2))
Y=model.response(model)
x2=model.matrix(Y~X2[,1]+X2[,2]+X2[,1]:X2[,2]+I(X2[,2]^2)+X2[,1]:I(X2[,2]^2),model)

colnames(x2)[colnames(x2) == "X2[, 1]"] = "trt"
colnames(x2)[colnames(x2) == "X2[, 2]"] = "time"
colnames(x2)[colnames(x2) == "X2[, 1]:X2[, 2]"] = "trt:time"
colnames(x2)[colnames(x2) == "I(X2[, 2]^2)"] = "time2"
colnames(x2)[colnames(x2) == "X2[, 1]:I(X2[,2]^2)"] = "trt:time2"

x2 = x2[order(x2[, "trt"], x2[, "time"]), ]

sourceCpp("code/covreg.cpp") 
fit_mcmc_cpp2 = MCMC(Y, x2, 1, 200000, 10)


B.psamp2=fit_mcmc_cpp2$B.psamp
A.psamp2=fit_mcmc_cpp2$A.psamp
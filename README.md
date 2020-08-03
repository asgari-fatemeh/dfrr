# dfrr [![Build Status](https://travis-ci.com/asgari-fatemeh/dfrr.svg?branch=master)](https://travis-ci.com/asgari-fatemeh/dfrr)
## Dichotomized functional response regression 
Implementing Function-on-Scalar Regression model, in which the response function is dichotomized and observed sparsely.

The output is a `dfrr-object`, which then can be injected into other methods/functions to postprocess the fitted model, including: `coef.dfrr`,`fitted.dfrr`, `residuals.dfrr`, `predict.dfrr`, `fpca.dfrr`, `summary.dfrr`, `model.matrix.dfrr`, `plot.dfrr`, `plot.coef.dfrr`, `plot.fitted.dfrr`, `plot.residuals.dfrr`, `qq.dfrr`, `plot.predict.dfrr`, `plot.fpca.dfrr`

## Installation

After installing `devtools' first, type in the following in R

```r
require(devtools)
install_github("asgari-fatemeh/dfrr")
```

## Usage
```r
dfrr(
  formula,
  yind = NULL,
  data = NULL,
  ydata = NULL,
  method = c("REML", "ML"),
  rangeval = NULL,
  basis = NULL,
  ...
)
```

## Example
```r
set.seed(2000)
N<-50;M<-24
X<-rnorm(N,mean=0)
time<-seq(0,1,length.out=M)
Y<-simulate_simple_dfrr(beta0=function(t){cos(pi*t+pi)},
                        beta1=function(t){2*t},
                        X=X,time=time)
dfrr_fit<-dfrr(Y~X,yind=time)
plot(dfrr_fit)

##### Fitting dfrr model to the Madras Longitudinal Schizophrenia data
data(madras)

ids<-unique(madras$id)
N<-length(ids)

ydata<-data.frame(.obs=madras$id,.index=madras$month,.value=madras$y)
xdata<-data.frame(Age=rep(NA,N),Gender=rep(NA,N))
for(i in 1:N){
  dt<-madras[madras$id==ids[i],]
  xdata[i,]<-c(dt$age[1],dt$gender[1])
}
rownames(xdata)<-ids

madras_dfrr<-dfrr(Y~Age+Gender+Age*Gender, data=xdata, ydata=ydata, J=11)
coefs<-coef(madras_dfrr)
plot(coefs)

fpcs<-fpca(madras_dfrr)
plot(fpcs)
```

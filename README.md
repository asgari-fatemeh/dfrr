# dfrr
## Dichotomized functional response regression 
Implementing Function-on-Scalar Regression model, in which the response function is dichotomized and observed sparsely.

The output is a `dfrr-object`, which then can be injected into other methods/functions to postprocess the fitted model, including: coefs.dfrr,fitted.dfrr, residuals.dfrr, predict.dfrr, eigen.dfrr, summary.dfrr, qq.dfrr, model.matrix.dfrr, plot.coefs.dfrr, plot.fitted.dfrr, plot.residuals.dfrr, plot.predict.dfrr, plot.eigen.dfrr, plot.residuals.dfrr

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
Y<-simulate.simple.dfrr(beta0=function(t){cos(pi*t+pi)},
                        beta1=function(t){2*t},
                        X=X,time=time)
dfrr_fit<-dfrr(Y~X,yind=time)
summary(dfrr_fit)

##### Fitting dfrr model to the Madras Longitudinal Schizophrenia data
data(madras)
ydata<-data.frame(.obs=madras$id,.index=madras$month,.value=madras$y)
ids<-unique(madras$id)
q<-4
N<-length(ids)
xData<-data.frame(Age=rep(NA,N),Gender=rep(NA,N))
for(i in 1:N){
  dt<-madras[madras$id==ids[i],]
  xData[i,]<-c(dt$age[1],dt$gender[1])
}
rownames(xData)<-ids

madras_dfrr<-dfrr(Y~Age+Gender+Age*Gender, data=xData, ydata=ydata, J=11,T_E=5)
coefs<-coef(madras_dfrr)
plot(coefs)

fpcs<-fpca(madras_dfrr)
plot(fpcs)
```

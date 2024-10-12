library(copula)
library(VineCopula)
library(MASS)
source("C:/Users/97623/Desktop/biFSD/GitHub/BI_FSDfunctions.R")
num=c(50,500)
#####
##Example B.1
for(n0 in num){
  result=c()
  for(i in 1:1000){
    mean1=c(0.85,0.85)
    sigma1=cbind(c(0.36,0.2),c(0.2,0.36))
    mean2=c(0.85,0.85)
    sigma2=cbind(c(0.36,0.2),c(0.2,0.36))
    X=mvrnorm(n0,mean1,sigma1) # generate F1
    Y=mvrnorm(n0,mean2,sigma2) # generate F2
    result=rbind(result,test(i,X=X,Y=Y))
  }
  print(colMeans(result))
}

##Example B.2
for(n0 in num){
  result=c()
  for(i in 1:1000){
    mean1=c(0.6,0.6)
    sigma1=cbind(c(0.64,0.2),c(0.2,0.64))
    mean2=c(0.8,0.8)
    sigma2=cbind(c(0.36,0.2),c(0.2,0.36))
    X=mvrnorm(n0,mean1,sigma1) # generate F1
    Y=mvrnorm(n0,mean2,sigma2) # generate F2
    result=rbind(result,test_cn(i,X=X,Y=Y))
  }
  print(colMeans(result))
}

##Example B.3
for(n0 in num){
  result=c()
  for(i in 1:1000){
    mean1=c(0.85,0.85)
    sigma1=cbind(c(0.36,-0.2),c(-0.2,0.36))
    mean2=c(0.85,0.85)
    sigma2=cbind(c(0.36,0.2),c(0.2,0.36))
    X=mvrnorm(n0,mean1,sigma1) # generate F1
    Y=mvrnorm(n0,mean2,sigma2) # generate F2
    result=rbind(result,test_cn(i,X=X,Y=Y))
  }
  print(colMeans(result))
}

##Example B.4
for(n0 in num){
  result=c()
  for(i in 1:1000){
    mean1=c(0.65,2.1)
    sigma1=cbind(c(0.36,0.2),c(0.2,0.36))
    mean2=c(0.85,0.85)
    sigma2=cbind(c(0.36,0.2),c(0.2,0.36))
    X=mvrnorm(n0,mean1,sigma1) # generate F1
    Y=mvrnorm(n0,mean2,sigma2) # generate F2
    result=rbind(result,test_cn(i,X=X,Y=Y))
  }
  print(colMeans(result))
}
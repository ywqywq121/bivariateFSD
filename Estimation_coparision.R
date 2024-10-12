library(VineCopula)  # For BiCopSelect, BiCopCDF, and rCopula
library(MASS)        # For kde2d
library(copula)
library(mvtnorm)
num=c(50,100,200,300,500,1000)
#####
#comparision of Frank, Clayton, Gumbel, Joe copulas
compare=function(i,cop,n){
  bicop_est=function(data,x){
    a1=pobs(data[,1])
    b1=pobs(data[,2])
    selectedCopula1=BiCopSelect(a1,b1,familyset = 1:10)
    Cop1=BiCop(family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
    oneecdf=function(x,X){
      tran=X-x
      sig=as.matrix(tran<=0) ##the number of X<x
      num=sum(sig)
      return(num/length(X))
    }
    u1=oneecdf(x=x[1],X=data[,1])
    u2=oneecdf(x=x[2],X=data[,2])
    result=BiCopCDF(u1=u1,u2=u2,family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
    return(result)
  }
  data=rCopula(n, cop)
  a1=pobs(data[,1])
  b1=pobs(data[,2])
  selectedCopula1=BiCopSelect(a1,b1,familyset = 1:10)
  Cop1=BiCop(family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
  oneecdf=function(x,X){
    tran=X-x
    sig=as.matrix(tran<=0) ##the number of X<x
    num=sum(sig)
    return(num/length(X))
  }
  G=150##the number of grids
  grids=seq(from=0,to=1,length.out=G)
  ##KDE
  h <- c(
    1.06 * sd(data[,1]) * n^(-1/5),
    1.06 * sd(data[,2]) * n^(-1/5)
  )
  kde_result <- kde2d(data[,1], data[,2],h=h, n = G,lims=c(0,1,0,1))
  x_grid <- kde_result$x
  y_grid <- kde_result$y
  z_density <- kde_result$z

  ISE1=0
  ISE2=0
  ISE3=0
  for(i in 1:G){
    for(j in 1:G){
      u1=oneecdf(x=grids[i],X=data[,1])
      u2=oneecdf(x=grids[j],X=data[,2])
      bi_est=BiCopCDF(u1=u1,u2=u2,family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
      ISE1=ISE1+(bi_est-pCopula(u=c(grids[i],grids[j]),Cop))^2
      cdf1 <- sum(z_density[1:i, 1:j]) * (x_grid[2] - x_grid[1]) * (y_grid[2] - y_grid[1])
      ISE2=ISE2+(cdf1-pCopula(u=c(grids[i],grids[j]),Cop))^2
      cdf2=mean(data[,1] <= grids[i] & data[,2] <= grids[j])
      ISE3=ISE3+(cdf2-pCopula(u=c(grids[i],grids[j]),Cop))^2
    }
  }
  return(c(ISE1,ISE2,ISE3))
}

##Gumbel Copula
for(par in c(2,3,4,5)){
  Cop=gumbelCopula(param=par,dim=2)
  for(n0 in num){
    result=c()
    for( i in 1:1000){
      result=rbind(result,compare(i,cop=Cop,n=n0))
      print(colMeans(result))
    }
  }
}
##Joe Copula
for(par in c(2,3,4,5)){
  Cop=joeCopula(param=par,dim=2)
  for(n0 in num){
    result=c()
    for( i in 1:1000){
      result=rbind(result,compare(i,cop=Cop,n=n0))
      print(colMeans(result))
    }
  }
}

##Frank Copula
for(par in c(1,2,3,4)){
  Cop=frankCopula(param=par,dim=2)
  for(n0 in num){
    result=c()
    for( i in 1:1000){
      result=rbind(result,compare(i,cop=Cop,n=n0))
      print(colMeans(result))
    }
  }
}
##clayton Copula
for(par in c(1,2,3,4)){
  Cop=claytonCopula(param=par,dim=2)
  for(n0 in num){
    result=c()
    for( i in 1:1000){
      result=rbind(result,compare(i,cop=Cop,n=n0))
      print(colMeans(result))
    }
  }
}
##### 
##bivariate normal
compare_norm=function(i,mus,sigmas,n){
  bicop_est=function(data,x){
    a1=pobs(data[,1])
    b1=pobs(data[,2])
    selectedCopula1=BiCopSelect(a1,b1,familyset = 1:10)
    Cop1=BiCop(family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
    oneecdf=function(x,X){
      tran=X-x
      sig=as.matrix(tran<=0) ##the number of X<x
      num=sum(sig)
      return(num/length(X))
    }
    u1=oneecdf(x=x[1],X=data[,1])
    u2=oneecdf(x=x[2],X=data[,2])
    result=BiCopCDF(u1=u1,u2=u2,family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
    return(result)
  }
  data=mvrnorm(n,mus,sigmas) # generate F2
  a1=pobs(data[,1])
  b1=pobs(data[,2])
  selectedCopula1=BiCopSelect(a1,b1,familyset = 1:10)
  Cop1=BiCop(family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
  oneecdf=function(x,X){
    tran=X-x
    sig=as.matrix(tran<=0) ##the number of X<x
    num=sum(sig)
    return(num/length(X))
  }
  G=150##the number of grids
  ##upper and lower bound
  min1=mus[1]-8*sqrt(sigmas[1,1])
  max1=mus[1]+8*sqrt(sigmas[1,1])
  min2=mus[2]-8*sqrt(sigmas[2,2])
  max2=mus[2]+8*sqrt(sigmas[2,2])
  grids1=seq(from=min1,to=max1,length.out=G)
  grids2=seq(from=min2,to=max2,length.out=G)
  h <- c(
    1.06 * sd(data[,1]) * n^(-1/5),
    1.06 * sd(data[,2]) * n^(-1/5)
  )
  kde_result <- kde2d(data[,1], data[,2],h=h, n = G,lims=c(min1,max1,min2,max2))
  x_grid <- kde_result$x
  y_grid <- kde_result$y
  z_density <- kde_result$z
  
  ISE1=0
  ISE2=0
  ISE3=0
  for(i in 1:G){
    for(j in 1:G){
      u1=oneecdf(x=grids1[i],X=data[,1])
      u2=oneecdf(x=grids2[j],X=data[,2])
      bi_est=BiCopCDF(u1=u1,u2=u2,family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
      ISE1=ISE1+(bi_est-pmvnorm(upper=c(grids1[i],grids2[j]),mean=mus,sigma=sigmas))^2
      cdf1 <- sum(z_density[1:i, 1:j]) * (x_grid[2] - x_grid[1]) * (y_grid[2] - y_grid[1])
      ISE2=ISE2+(cdf1-pmvnorm(upper=c(grids1[i],grids2[j]),mean=mus,sigma=sigmas))^2
      cdf2=mean(data[,1] <= grids1[i] & data[,2] <= grids2[j])
      ISE3=ISE3+(cdf2-pmvnorm(upper=c(grids1[i],grids2[j]),mean=mus,sigma=sigmas))^2
    }
  }
  return(c(ISE1,ISE2,ISE3))
}
##Norm1
for(n0 in num){
  result=c()
  for(i in 1:1000){
    mean1=c(0.85,0.85)
    sigma1=cbind(c(0.36,0.2),c(0.2,0.36))
    result=rbind(result,compare_norm(i,mus=mean1,sigmas=sigma1,n=n0))
  }
  print(colMeans(result))
}
##Norm2
for(n0 in num){
  result=c()
  for(i in 1:1000){
    mean1=c(0.6,0.6)
    sigma1=cbind(c(0.64,0.2),c(0.2,0.64))
    result=rbind(result,compare_norm(i,mus=mean1,sigmas=sigma1,n=n0))
  }
  print(colMeans(result))
}
##Norm3
for(n0 in num){
  result=c()
  for(i in 1:1000){
    mean1=c(0.85,0.85)
    sigma1=cbind(c(0.36,-0.2),c(-0.2,0.36))
    result=rbind(result,compare_norm(i,mus=mean1,sigmas=sigma1,n=n0))
  }
  print(colMeans(result))
}
##Norm4
for(n0 in num){
  result=c()
  for(i in 1:1000){
    mean1=c(0.65,2.1)
    sigma1=cbind(c(0.36,0.2),c(0.2,0.36))
    result=rbind(result,compare_norm(i,mus=mean1,sigmas=sigma1,n=n0))
  }
  print(colMeans(result))
}

#####
##mixed copula
compare_mix=function(i,cop1,cop2,n){
  bicop_est=function(data,x){
    a1=pobs(data[,1])
    b1=pobs(data[,2])
    selectedCopula1=BiCopSelect(a1,b1,familyset = 1:10)
    Cop1=BiCop(family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
    oneecdf=function(x,X){
      tran=X-x
      sig=as.matrix(tran<=0) ##the number of X<x
      num=sum(sig)
      return(num/length(X))
    }
    u1=oneecdf(x=x[1],X=data[,1])
    u2=oneecdf(x=x[2],X=data[,2])
    result=BiCopCDF(u1=u1,u2=u2,family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
    return(result)
  }
  data1=rCopula(n, cop1)
  data2=rCopula(n, cop2)
  z=sample(c(0,1),n,TRUE)
  data=data1*z+(1-z)*data2
  a1=pobs(data[,1])
  b1=pobs(data[,2])
  selectedCopula1=BiCopSelect(a1,b1,familyset = 1:10)
  Cop1=BiCop(family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
  oneecdf=function(x,X){
    tran=X-x
    sig=as.matrix(tran<=0) ##the number of X<x
    num=sum(sig)
    return(num/length(X))
  }
  G=150##the number of grids
  grids=seq(from=0,to=1,length.out=G)
  h <- c(
    1.06 * sd(data[,1]) * n^(-1/5),
    1.06 * sd(data[,2]) * n^(-1/5)
  )
  kde_result <- kde2d(data[,1], data[,2],h=h, n = G,lims=c(0,1,0,1))
  x_grid <- kde_result$x
  y_grid <- kde_result$y
  z_density <- kde_result$z

  ISE1=0
  ISE2=0
  ISE3=0
  for(i in 1:G){
    for(j in 1:G){
      u1=oneecdf(x=grids[i],X=data[,1])
      u2=oneecdf(x=grids[j],X=data[,2])
      bi_est=BiCopCDF(u1=u1,u2=u2,family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
      ISE1=ISE1+(bi_est-0.5*pCopula(u=c(grids[i],grids[j]),cop1)-0.5*pCopula(u=c(grids[i],grids[j]),cop2))^2
      cdf1 <- sum(z_density[1:i, 1:j]) * (x_grid[2] - x_grid[1]) * (y_grid[2] - y_grid[1])
      ISE2=ISE2+(cdf1-0.5*pCopula(u=c(grids[i],grids[j]),cop1)-0.5*pCopula(u=c(grids[i],grids[j]),cop2))^2
      cdf2=mean(data[,1] <= grids[i] & data[,2] <= grids[j])
      ISE3=ISE3+(cdf2-0.5*pCopula(u=c(grids[i],grids[j]),cop1)-0.5*pCopula(u=c(grids[i],grids[j]),cop2))^2
    }
  }
  return(c(ISE1,ISE2,ISE3))
}
for(n0 in num){
  Cop1=gumbelCopula(param=4,dim=2)
  Cop2=claytonCopula(param=5,dim=2)
  result=c()
  for(i in 1:1000){
    result=rbind(result,compare_mix(i,cop1=Cop1,cop2=Cop2,n=n0))
  }
  print(colMeans(result))
}

for(n0 in num){
  Cop1=gumbelCopula(param=5,dim=2)
  Cop2=frankCopula(param=5,dim=2)
  result=c()
  for(i in 1:1000){
    result=rbind(result,compare_mix(i,cop1=Cop1,cop2=Cop2,n=n0))
  }
  print(colMeans(result))
}

for(n0 in num){
  Cop1=gumbelCopula(param=3,dim=2)
  Cop2=joeCopula(param=4,dim=2)
  result=c()
  for(i in 1:1000){
    result=rbind(result,compare_mix(i,cop1=Cop1,cop2=Cop2,n=n0))
  }
  print(colMeans(result))
}
##grid search
grid_search_2d <- function(f, x_min, x_max, y_min, y_max, n) {
  step_x <- (x_max - x_min) / (n - 1)
  step_y <- (y_max - y_min) / (n - 1)
  x_values <- seq(x_min, x_max, length.out = n)
  y_values <- seq(y_min, y_max, length.out = n)
  
  max_value <- -Inf
  x_max_value <- x_min
  y_max_value <- y_min
  
  for (x in x_values) {
    for (y in y_values) {
      val <- f(x, y)
      if (val > max_value) {
        max_value <- val
        x_max_value <- x
        y_max_value <- y
      }
    }
  }
  
  return(max_value)
}
### grid search on specific set 
modified_gridsearch <- function(x_min, x_max, y_min, y_max, n, 
                                objective, condition_function, cn_vector) {
  step_x <- (x_max - x_min) / (n - 1)
  step_y <- (y_max - y_min) / (n - 1)
  x_values <- seq(x_min, x_max, length.out = n)
  y_values <- seq(y_min, y_max, length.out = n)
  
  max_vals <- rep(-Inf, length(cn_vector))
  
  for (x in x_values) {
    for (y in y_values) {
      for (k in seq_along(cn_vector)) {
        cn <- cn_vector[k]
        cond_val <- condition_function(x, y)
        if (abs(cond_val) < cn) {
          obj_val <- objective(x, y)
          if (obj_val > max_vals[k]) {
            max_vals[k] <- obj_val
          }
        }
      }
    }
  }
  
  return(max_vals)
}
###set different r_mn
cns=function(m,n){
  if(m==n){
    taun=n
  }else{
    taun=m*n/(m+n)
  }
  cn=c(taun^(1/3),taun^(1/5),taun^(1/7),taun^(1/9),log(taun),log(log(taun)))
  return(cn)
}
###comparison of performance of  different r_mn
test_cn=function(i,X,Y){
  F1=X # generate F1
  F2=Y # generate F2
  m=length(F1[,1])
  n=length(F2[,1])
  rhon=cns(m,n)
  B=500
  alpha1=0.1
  alpha2=0.05
  alpha3=0.01
  data=rbind(F1,F2)
  a1=pobs(F1[,1])
  b1=pobs(F1[,2])
  selectedCopula1=BiCopSelect(a1,b1,familyset = 1:10)
  Cop1=BiCop(family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
  a2=pobs(F2[,1])
  b2=pobs(F2[,2])
  selectedCopula2=BiCopSelect(a2,b2,familyset = 1:10)
  Cop2=BiCop(family=selectedCopula2$family,par=selectedCopula2$par,par2=selectedCopula2$par2)
  twocopula=function(x,y,X,family,par,par2){
    
    u1=oneecdf(x=x,X=X[,1])
    u2=oneecdf(x=y,X=X[,2])
    return(BiCopCDF(u1=u1,u2=u2,family=family,par=par,par2=par2))
  }
  oneecdf=function(x,X){
    tran=X-x
    sig=as.matrix(tran<=0) ##the number of X<x
    num=sum(sig)
    return(num/length(X))
  }
  ##Grid Algorithm
  min1=-0.5
  min2=-0.5
  max1=1.5
  max2=1.5
  steps=20
  copula_statistic=function(x,y){
    F1_hat=twocopula(x=x,y=y,X=F1,family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
    F2_hat=twocopula(x=x,y=y,X=F2,family=selectedCopula2$family,par=selectedCopula2$par,par2=selectedCopula2$par2)
    result=sqrt(n*m/(m+n))*(F1_hat-F2_hat)
    return(result)
  }
  copulasta=grid_search_2d(copula_statistic,min1,max1,min2,max2,steps)
  bootcopula=c()
  for(i in 1:B){
    index1=sample(x=1:m,size=m,replace = TRUE)
    index2=sample(x=1:n,size=n,replace = TRUE)
    X3=F1[index1,]
    X4=F2[index2,]
    a3=pobs(X3[,1])
    b3=pobs(X3[,2])
    para1 =BiCopEst(u1=a3,u2=b3,family=selectedCopula1$family)
    a4=pobs(X4[,1])
    b4=pobs(X4[,2])
    para2=BiCopEst(u1=a4,u2=b4,family=selectedCopula2$family)
    Cop3=BiCop(family=para1$family,par=para1$par,par2=para1$par2)
    Cop4=BiCop(family=para2$family,par=para2$par,par2=para2$par2)  
    statistic2=function(x,y){
      boot1=twocopula(x=x,y=y,X=X3,family=para1$family,par=para1$par,par2=para1$par2)
      boot2=twocopula(x=x,y=y,X=X4,family=para2$family,par=para2$par,par2=para2$par2)
      return(sqrt(n*m/(m+n))*(boot1-boot2)-copula_statistic(x,y))
    }
    bootstar=modified_gridsearch(min1,max1,min2,max2,steps,statistic2,copula_statistic,rhon)
    bootcopula=rbind(bootcopula,bootstar)
  }
  p=numeric(length(rhon))
  rej1=numeric(length(rhon))
  rej2=numeric(length(rhon))
  rej3=numeric(length(rhon))
  for(i in 1:length(rhon)){
    p[i]=mean(bootcopula[,i]>=copulasta)
    rej1[i]=as.numeric(p[i]<alpha1)
    rej2[i]=as.numeric(p[i]<alpha2)
    rej3[i]=as.numeric(p[i]<alpha3)
  }
  return(c(p,rej1,rej2,rej3))
}
###choose specific r_mn to test two populations
test=function(i,X,Y){
  F1=X # generate F1
  F2=Y # generate F2
  m=length(F1[,1])
  n=length(F2[,1])
  if(m==n){
    taun=n
  }else{
    taun=m*n/(m+n)
  }
  rhon=taun^(1/3)
  B=500
  alpha1=0.1
  alpha2=0.05
  alpha3=0.01
  data=rbind(F1,F2)
  a1=pobs(F1[,1])
  b1=pobs(F1[,2])
  selectedCopula1=BiCopSelect(a1,b1,familyset = 1:10)
  Cop1=BiCop(family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
  a2=pobs(F2[,1])
  b2=pobs(F2[,2])
  selectedCopula2=BiCopSelect(a2,b2,familyset = 1:10)
  Cop2=BiCop(family=selectedCopula2$family,par=selectedCopula2$par,par2=selectedCopula2$par2)
  twocopula=function(x,y,X,family,par,par2){
    
    u1=oneecdf(x=x,X=X[,1])
    u2=oneecdf(x=y,X=X[,2])
    return(BiCopCDF(u1=u1,u2=u2,family=family,par=par,par2=par2))
  }
  oneecdf=function(x,X){
    tran=X-x
    sig=as.matrix(tran<=0) ##the number of X<x
    num=sum(sig)
    return(num/length(X))
  }
  ##Grid Algorithm
  min1=-0.5
  min2=-0.5
  max1=1.5
  max2=1.5
  steps=20
  copula_statistic=function(x,y){
    F1_hat=twocopula(x=x,y=y,X=F1,family=selectedCopula1$family,par=selectedCopula1$par,par2=selectedCopula1$par2)
    F2_hat=twocopula(x=x,y=y,X=F2,family=selectedCopula2$family,par=selectedCopula2$par,par2=selectedCopula2$par2)
    result=sqrt(n*m/(m+n))*(F1_hat-F2_hat)
    return(result)
  }
  copulasta=grid_search_2d(copula_statistic,min1,max1,min2,max2,steps)
  bootcopula=c()
  for(i in 1:B){
    index1=sample(x=1:m,size=m,replace = TRUE)
    index2=sample(x=1:n,size=n,replace = TRUE)
    X3=F1[index1,]
    X4=F2[index2,]
    a3=pobs(X3[,1])
    b3=pobs(X3[,2])
    para1 =BiCopEst(u1=a3,u2=b3,family=selectedCopula1$family)
    a4=pobs(X4[,1])
    b4=pobs(X4[,2])
    para2=BiCopEst(u1=a4,u2=b4,family=selectedCopula2$family)
    Cop3=BiCop(family=para1$family,par=para1$par,par2=para1$par2)
    Cop4=BiCop(family=para2$family,par=para2$par,par2=para2$par2)  
    statistic2=function(x,y){
      boot1=twocopula(x=x,y=y,X=X3,family=para1$family,par=para1$par,par2=para1$par2)
      boot2=twocopula(x=x,y=y,X=X4,family=para2$family,par=para2$par,par2=para2$par2)
      return(sqrt(n*m/(m+n))*(boot1-boot2)-copula_statistic(x,y))
    }
    bootstar=modified_gridsearch(min1,max1,min2,max2,steps,statistic2,copula_statistic,rhon)
    bootcopula=rbind(bootcopula,bootstar)
  }
  p=numeric(length(rhon))
  rej1=numeric(length(rhon))
  rej2=numeric(length(rhon))
  rej3=numeric(length(rhon))
  for(i in 1:length(rhon)){
    p[i]=mean(bootcopula[,i]>=copulasta)
    rej1[i]=as.numeric(p[i]<alpha1)
    rej2[i]=as.numeric(p[i]<alpha2)
    rej3[i]=as.numeric(p[i]<alpha3)
  }
  return(c(p,rej1,rej2,rej3))
}

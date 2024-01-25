#' @importFrom stats as.formula binomial predict sd
NULL
#' Fit the interval-censored AFT model with quantile linear model
#' 
#' Fit the interval-censored AFT model with quantile linear model with general interval-censored data
#'
#' @param L left-censoring time, having 0 if left-censored.
#' @param R right-censoring time, having \code{Inf} if right-censored.
#' @param x X matrix of baseline covariates.
#' @param beta0 true parameter values, including intercept, with the default being 1 for all parameters.
#' @param type penalized estimating method, default is "scad", or "lasso", "alasso" penalties are applicable.
#' @param selection method to use variable selection method, default is "GCV", or "BIC", "AIC" methods are applicable.
#' @param lamb.len the range of the lambda grid, default is 50.
#' @param lambmin the minimum value of lambda, default is 2^(-10).
#' @param r the value provided in the SCAD penalty. The default value is 3.7.
#' @param outlier a logical value indicating whether to replace outlier data with some of the normal data. If set to "TRUE," outliers will be replaced.
#' @param outlier.rate the outlier rate to set for the normal data, default is 0.2.
#' @param outlier.value The value to add to the normal data in order to generate outlier values, default is 5.
#' @param weight weights of covariate x and response variable y. The default is 1.
#' @param tol tolerance of estimated as 0 and calculated. The default value is 1e-10.
#'
#' @return \code{rankICpen} returns a data frame containing at least the following components:
#' \itemize{
#'   \item \code{est}: regression estimator.
#'   \item \code{cor}: the number of correctly estimated non-zero values for \code{est}.
#'   \item \code{incor}: the number of correctly estimated zero values for \code{est}.
#'   \item \code{amad}: calculted AMAD values related abstract bias value for \code{est}.
#'   \item \code{mrme}: calculted MRME values related squared bias value for \code{est}.
#'   \item \code{optmeasure}: optimalized measured values for \code{est}.
#'   \item \code{optlamb}: optimalized lambda values for \code{est}.
#'   \item \code{beta}: optimalized beta value.
#' }
#'
#' @details
#' see Kim et al., (2024+) for detailed method explanation.
#'
#' @references
#' 
#' Kim, Y., Park, S., Choi, S. (2024+). Rank-based variable selection with interval-censored data.
#' 
#'
#' @examples
#' \dontrun{
#' # Simulations
#' set.seed(111)
#' n = 200
#' x1 = runif(n,-1,1)
#' x2 = rbinom(n,1,0.43)
#' x = cbind(x1,x2)
#' T = 4 + x1 + x2 + rnorm(n)
#' U = 4 + (1 - 0.25*x1)*runif(n, -6, 5)
#' V = U + (1 - 0.1*x2)*runif(n, 6, 20)
#' U = exp(dplyr::case_when(TRUE ~ T, T>V ~ V, T<U ~ -Inf))
#' V = exp(dplyr::case_when(TRUE ~ T, T>V ~ Inf, T<U ~ U))
#' type="scad"; selection="BIC"; outlier=FALSE
#' rankICpen(L=U,R=V,x=x,type=type,selection=selection,outlier=outlier)
#' 
#' 
#' # Data example
#' library(PICBayes)
#' library(tidyverse)
#' data("mCRC")
#' d = with(data.frame(mCRC), data.frame(U = ifelse(y==0,R,L),
#'                                       V = ifelse(y==2,L,R),
#'                                       # Cluster weighted data
#'                                       id=(rep(c(table(SITE)),c(table(SITE)))),
#'                                       # Treatment arm: 0 = FOLFIRI alone, 1 = Panitumumab + FOLFIRI.
#'                                       x1= case_when(TRT_C == 0 ~ 0, #Pan et al data
#'                                                     TRT_C == 1 ~ 1),
#'                                       # Tumor KRAS mutation status: 0 = wild-type, 1 = mutant.
#'                                       x2= case_when(KRAS_C == 0 ~ 1,
#'                                                     KRAS_C == 1 ~ 0),
#'                                       delta = case_when(IC == 0 ~ 1,
#'                                                         IC == 1 ~ 4)
#'));
#' U=(log(d$U));V=log(d$V); type="scad"; selection="BIC"; outlier=FALSE
#' x = cbind(d$x1,d$x2); id=d$id;  tau=0.1;
#' rankICpen(L=U,R=V,x=x,type=type,selection=selection,outlier=outlier)
#' }
#' @export
#'
#'
#'

rankICpen=function(L,R,x,beta0=rep(1,ncol(x)),type=NULL,selection=NULL,lamb.len=50,lambmin=NULL,r=3.7,outlier="without",outlier.rate=0.2, outlier.value=5,weight = rep(1,length(L)),tol=1e-10){
  
  library(glmnet)
  library(pracma)
  library(dplyr)
  library(tidyverse)
  require(quantreg)
  
  aftrq=function(L,R, x, weight = rep(1,length(L))){
    # options(warn=-1)
    n = nrow(x); p=ncol(x); 
    id_i = rep(1:n, each=n)
    id_j = rep(1:n, times=n)
    wi = weight[id_i]
    Li = L[id_i]; Lj = L[id_j]
    Ri = R[id_i]; Rj = R[id_j]
    xi = x[id_i,]; xj = x[id_j,]
    idd = which(Ri != Inf & Lj != -Inf)
    yy = ((Ri - Lj)*wi)[idd]
    xx = ((xi - xj)*wi)[idd,]
    yy.old = c(yy, 1e10); xx.old=rbind(as.matrix(xx), -colSums(xx))
    rq.fit(x=xx.old,y=yy.old)$coef
  }
  
  aftrq_pen=function(L,R,x, beta, type=NULL, lambda = NULL, lambda2 = NULL, weight = rep(1,length(L))){
    # options(warn=-1)
    n = length(L); p=ncol(x); 
    id_i = rep(1:n, each=n); id_j = rep(1:n, times=n); wi = weight[id_i]
    Li = L[id_i]; Lj = L[id_j]
    Ri = R[id_i]; Rj = R[id_j]
    xi = x[id_i,]; xj = x[id_j,]
    idd = which(Ri != Inf & Lj != -Inf)
    yy = ((Ri - Lj)*wi)[idd]
    xx = ((xi - xj)*wi)[idd,]
    if(type=="oracle"){
      pen = penft(type=type,p=p,lambda=lambda)
    }else{
      pen = penft(type=type,p=p,lambda=lambda,beta=(beta))
    }
    yy.old = c(yy, 1e10); xx.old=rbind(as.matrix(xx), -colSums(xx))
    yy.new = c(yy.old/(n^2), rep(0, p)); xx.new = rbind(xx.old/(n^2), (n)*pen)
    rq.fit(x=xx.new,y=yy.new)$coef
  }
  
  penft = function(type,p,lambda=NULL,beta=NULL) {
    if (type == "lasso") res = diag(case_when(beta > 0 ~ lambda,
                                              beta < 0 ~  -lambda,
                                              TRUE ~ 0))
    if (type == "alasso") res = diag(case_when(beta > 0 ~ (lambda/(abs(beta)+0.05)),
                                               beta < 0 ~  (-lambda/(abs(beta)+0.05)),
                                               TRUE ~ 0))
    if (type == "scad") {
      res = diag(case_when(abs(beta) <= lambda ~ (lambda),
                           abs(beta) >= (r*lambda) ~ 0,
                           TRUE ~  ((((r*lambda) - abs(beta)))/((r-1))) ))
    }
    if (type == "oracle") res = diag(p)*0
    res
  }
  
  # Umatrix
  Efunc=function(L,R,beta,x){
    # options(warn=-1)
    n = length(L); p = ncol(x); n=length(L)
    id_i = rep(1:n, each=n)
    id_j = rep(1:n, times=n)
    wi = weight[id_i]
    Li = L[id_i]; Lj = L[id_j]
    Ri = R[id_i]; Rj = R[id_j]
    xi = x[id_i,]; xj = x[id_j,]
    delta1i = ifelse(Ri != Inf,1,0); delta2j = ifelse(Lj != -Inf,1,0);
    idd = which(Ri != Inf & Lj != -Inf)
    yy = ((Ri - Lj)*wi)[idd]
    xx = ((xi - xj)*wi)[idd,]
    ind = ifelse( ((Ri - Lj)*wi) - ((xi - xj)*wi) %*% beta<=0, 1, 0); length(ind)
    (delta.ind=delta1i * delta2j * ind); (delta2.ind=delta2j * ind)
    delta.indmat = matrix(rep(delta.ind,p),ncol = p); delta2.indmat = matrix(rep(delta2.ind,p),ncol = p)
    (num=delta.indmat * xj) #num
    (denum=delta2.indmat) #denum
    
    sum.new.xvec=NULL; 
    for (i in 0:n) {
      if((n+i*n)==n*(n+1)) break
      ximat = (xi*delta1i)[(1+i*n):(n+i*n)]
      denum_sum = c(sum(denum[(1+i*n):(n+i*n)]))+0.05 # n length
      num_eachsum = c(colSums(num[(1+i*n):(n+i*n),])) # nrow=n => 1 X p
      num_sum = matrix(rep(num_eachsum,each=n),ncol=p) # n X p
      newxmat = (ximat - (num_sum/denum_sum)) # n X p
      new.xvec = colMeans(newxmat)
      sum.new.xvec = rbind(sum.new.xvec, new.xvec)
    }
    Emat=colSums(sum.new.xvec)/n
    Emat
  }
  
  Gfunc=function(L,R,beta,x){
    options(warn=-1)
    L=L; R=R
    n = length(L); p = ncol(x); n=length(L)
    id_i = rep(1:n, each=n)
    id_j = rep(1:n, times=n)
    wi = weight[id_i]
    Li = L[id_i]; Lj = L[id_j]
    Ri = R[id_i]; Rj = R[id_j]
    xi = x[id_i,]; xj = x[id_j,]
    delta1i = ifelse(Ri != Inf,1,0); delta2j = ifelse(Lj != -Inf,1,0);
    idd = which(Ri != Inf & Lj != -Inf)
    yy = ((Ri - Lj)*wi)[idd]
    xx = ((xi - xj)*wi)[idd,]
    
    ind = ifelse( ((Ri - Lj)*wi) - ((xi - xj)*wi) %*% beta<=0, 1, 0); length(ind)
    (delta.ind=delta1i * delta2j * ind); (delta2.ind=delta2j * ind)
    delta.indmat = matrix(rep(delta.ind,p),ncol = p); delta2.indmat = matrix(rep(delta2.ind,p),ncol = p)
    (num=delta.indmat * xj) #num
    (num2=delta2.indmat * xj) #num2
    (denum=delta2.indmat) #denum
    
    newxmat4=matrix(0,p,p)
    for (i in 0:n) {
      if((n+i*n)==n*(n+1)) break
      ximat = (xi*delta1i)[(1+i*n):(n+i*n)]
      ximat2 = (xi)[(1+i*n):(n+i*n)]
      denum_sum = c(sum(denum[(1+i*n):(n+i*n)]))+0.05 # n length
      num_eachsum = c(colSums(num[(1+i*n):(n+i*n),])) # nrow=n => 1 X p
      num_eachsum2 = c(colSums(num2[(1+i*n):(n+i*n),])) # nrow=n => 1 X p
      num_sum = matrix(rep(num_eachsum,each=n),ncol=p) # n X p
      num_sum2 = matrix(rep(num_eachsum2,each=n),ncol=p) # n X p
      newxmat = (ximat - (num_sum/denum_sum)) # n X p
      newxmat2 = (ximat2 - (num_sum2/denum_sum)) # n X p
      newxvec = colMeans(newxmat)
      newxvec2 = colMeans(newxmat2)
      newxmat3 = (newxvec)%*%t(newxvec2)
      newxmat4 = newxmat4 + (newxmat3) + diag(p)*0.05
    }
    Gmat=(newxmat4/n)
    Gmat
  }
  
  Lossfunc=function(L,R,beta,x){
    n = length(L); p = ncol(x); 
    id_i = rep(1:n, each=n)
    id_j = rep(1:n, times=n)
    wi = weight[id_i]
    Li = L[id_i]; Lj = L[id_j]
    Ri = R[id_i]; Rj = R[id_j]
    xi = x[id_i,]; xj = x[id_j,]
    delta1i = ifelse(Ri != Inf,1,0); delta2j = ifelse(Lj != -Inf,1,0);
    idd = which(Ri != Inf & Lj != -Inf)
    yy = ((Ri - Lj)*wi)[idd]
    xx = ((xi - xj)*wi)[idd,]
    res = yy- xx %*% beta
    ind = ifelse(res<=0,1,0)
    (sum(abs(res*ind))/(n))
  }
  
  Tstatft=function(L,R,beta,x){
    n=length(L)
    Eft = Efunc(L=L,R=R,beta=beta,x=x)
    Gft = Gfunc(L=L,R=R,beta=beta,x=x)
    (t(Eft)%*%solve(Gft)%*%(Eft))*n
  }
  
  BICft=function(L,R,beta,x){
    p=ncol(x); df=length(beta[(beta)<tol])
    Tstat = Tstatft(L=L,R=R,beta=beta,x=x); n=length(L)
    Tstat + log(n) * df
  }
  
  AICft=function(L,R,beta,x){
    p=ncol(x); df=length(beta[(beta)<tol])
    Tstat = Tstatft(L=L,R=R,beta=beta,x=x); n=length(L)
    Tstat + 2 * df
  }
  
  GCVft=function(L,R,beta,x){
    p=ncol(x); df=sum(beta[(beta)<tol]); n=length(L)
    Lfunc = Lossfunc(L=L,R=R,beta=(beta),x=x)
    ((Lfunc)/((1-(df/n))^2))
  }
  
  lamb_max=function(L,R,x,weight=rep(1,length(L))){
    # options(warn=-1); 
    n = length(L); p = ncol(x) 
    id_i = rep(1:n, each = n); id_j = rep(1:n, times = n)
    wi = weight[id_i]
    Li = L[id_i]; Lj = L[id_j]
    Ri = R[id_i]; Rj = R[id_j]
    xi = x[id_i,]; xj = x[id_j,]
    idd = which(Ri != Inf & Lj != -Inf)
    yy = ((Ri - Lj)*wi)[idd]
    xx = ((xi - xj)*wi)[idd,]
    
    yy.old = c(yy, 1e10); xx.old=rbind(as.matrix(xx), -colSums(xx))
    max(abs(t(xx.old /(n*n)) %*% sign(yy.old)))
  }
  
  aftrq_penlamb=function(L,R,x){
    
    if(is.null(lambmin)!=TRUE){lambmin=lambmin
    }else if(type=="lasso" ){lambmin=2^(-12)
    }else if(type=="alasso"){lambmin=2^(-12)
    }else if(type=="scad"){lambmin=2^(-2)}
    lambmax=lamb_max(L=L,R=R,x=x); 
    lamb=seq(lambmin,lambmax,length=lamb.len); 
    beta=optbic=optaic=optgcv=NULL;
    old_beta = aftrq(L=L,R=R,x=x)
    
    for (i in 1:length(lamb)) {
      new_beta=aftrq_pen(L=L,R=R,x=x,type=type,beta=old_beta,lambda=lamb[i])
      test.beta=new_beta
      # if(type=="scad")test.beta=abs(new_beta) else test.beta=new_beta
      if(selection=="BIC") optbic=rbind(optbic,BICft(L=L,R=R,x=x,beta=test.beta))
      if(selection=="AIC") optaic=rbind(optaic,AICft(L=L,R=R,x=x,beta=test.beta))
      if(selection=="GCV") optgcv=rbind(optgcv,GCVft(L=L,R=R,x=x,beta=test.beta))
      old_beta=new_beta
      beta=rbind(beta,new_beta)
    }
    
    if(selection=="BIC") {rankpen_beta=beta[which.min(optbic),]; measure=as.vector(optbic); optmeasure=as.vector(which.min(optbic))
    }else if(selection=="AIC") {rankpen_beta=beta[which.min(optaic),]; measure=as.vector(optaic); optmeasure=as.vector(which.min(optaic))
    }else if(selection=="GCV") {rankpen_beta=beta[which.min(optgcv),]; measure=as.vector(optgcv); optmeasure=as.vector(which.min(optgcv))
    }
    
    list(rankpen_beta=rankpen_beta, optmeasure=optmeasure, optlamb=cbind(lambda=lamb,measure=measure), beta=beta)
  }
  
  n=length(L); p=ncol(x);
  if(outlier=="with"){
    split1<- sample(c(rep(0, (1-outlier.rate) * length(L)), rep(1, outlier.rate * length(L))))
    L[split1 == 1] <- L[split1 == 1] + outlier.value
    R[split1 == 1] <- R[split1 == 1] + outlier.value
  }
  
  options(warn = -1)
  aftrq_result=aftrq_penlamb(L=L,R=R,x=x)
  new_beta = aftrq_result$rankpen_beta
  if(type=="oracle"){ new_beta[-c(beta0!=0)]=0}
  biasabs=abs(new_beta-beta0); amad_beta = mean((biasabs)); 
  mrme_beta <- (t((x)%*%biasabs)%*%(scale(x)%*%biasabs))/(nrow(x)); #mrme_beta <- (sum((biasabs)^2))/(nrow(x));
  new_beta_selected_vars <- abs(new_beta) > tol # new_beta_selected_vars <- new_beta != 0
  cor <- sum(new_beta_selected_vars[beta0!=0]); incor <- sum(new_beta_selected_vars) - cor
  list(est=new_beta, cor=cor, incor=incor, amad=amad_beta, mrme=mrme_beta,
       optmeasure= aftrq_result$optmeasure, optlamb = aftrq_result$optlamb, beta = aftrq_result$beta)
}

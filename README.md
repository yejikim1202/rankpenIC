# Introduction to `rankpenIC` R package




## Introduction
`rankpenIC` is the R package to introduce the $L_1$ regularized Rank estimator tailored for interval-censored data, aiming for simultaneous estimation and variable selection when the data are partially interval-censored that include doubly-censored (DC) data and partly interval-censored (PIC).
Let $T$ and $X$ be the event time of interest and its related $p$-vector of covariates, respectively.
Our main objective is to estimate 
the $p$-dimensional linear coefficient vector ${\boldsymbol{\beta}}_0$
in the following linear linear regression model:

$$T_i = {\bf x}_i^T {\boldsymbol{\beta}}_0 + \epsilon_i,\quad i=1, \ldots ,n, $$

where $\epsilon_i$ is the random error.
When the data are subject to partially interval-censoring, 
left and right endpoints of the censoring time, $L$ and $R$,
are observed instead of $T$ such that $T\in(L,R)$.
Note that double-censoring can also be viewed as 
a special case of partly interval-censoring, 
i.e., $T$ is left-censored if $L=0$ and right-censored if $R=\infty$. 



## Description
This R package `rankpenIC` introduces the $L_1$ regularized Rank estimator tailored for interval-censored data, aiming for simultaneous estimation and variable selection for (cluster-correlated) partially interval-censored data, which includes both double-censoring and partially interval-censoring.

Vignettes is available in [here](http://htmlpreview.github.io/?https://github.com/YejiStat/rankpenIC/blob/main/vignettes/rankpenIC.html).


## Usages 
```{r}
library(PICBayes)
data("mCRC")
d = with(data.frame(mCRC), data.frame(U = ifelse(y==0,R,L),
                                      V = ifelse(y==2,L,R),
                                      # Cluster weighted data
                                      id=(rep(c(table(SITE)),c(table(SITE)))),
                                      # Treatment arm: 0 = FOLFIRI alone, 1 = Panitumumab + FOLFIRI.
                                      x1= case_when(TRT_C == 0 ~ 0, #Pan et al data
                                                    TRT_C == 1 ~ 1),
                                      # Tumor KRAS mutation status: 0 = wild-type, 1 = mutant.
                                      x2= case_when(KRAS_C == 0 ~ 1,
                                                    KRAS_C == 1 ~ 0),
                                      delta = case_when(IC == 0 ~ 1,
                                                        IC == 1 ~ 4)
));
L=(log(d$U));R=log(d$V); delta=d$delta; x = cbind(d$x1,d$x2); id=d$id;
type=c("lasso","alasso","scad","oracle")[3]; selection=c("BIC","AIC","GCV")[3]; outlier=c("without","with")[2]; beta0=c(1,1,0,0,1,0,0,0);
rankICpen(L=U,R=V,x=x,type=type,selection=selection,outlier=outlier,beta0 = beta0)$est
```


## References

* Pan, C. (2021). 
PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. 
https://CRAN.R-project.org/package=PICBayes.

* Pak, D., Langohr, K., Ning, J., Cort ́es Mart ́ınez, J., G ́omez Melis, G., and Shen, Y. (2020). Modeling the coronavirus disease 2019 incubation period: impact on quarantine policy. Mathematics, 8(9):1631.

* Kim et al., (2024+). Rank-based variable selection with interval-censored data.

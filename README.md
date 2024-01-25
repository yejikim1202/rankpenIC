# Introduction to `rankpenIC` R package




## Introduction
`rankpenIC` is the R package to introduce the $L_1$ regularized Rank estimator tailored for interval-censored data, aiming for simultaneous estimation and variable selection when the data are partially interval-censored that include doubly-censored (DC) data and partly interval-censored (PIC) and possibly correlated within the same cluster.
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

```{r}
x=1
y=2
x+y
```

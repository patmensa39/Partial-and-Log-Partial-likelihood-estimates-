# Partial-and-Log-Partial-likelihood-estimates-
Fitting a Weibull distribution, Partial and Log Partial  likelihood estimates 



# Simple way to fit Weibull to gastricXelox data

library(asaur)

library(survival)

survfit.gastric.Xelox <- survfit(Surv(timeWeeks, delta)~1, data = gastricXelox) # Kaplan-Meier Estimates

summary(survfit.gastric.Xelox)





ss <- survfit.gastric.Xelox$surv # Extracts the survival probabilities 

tt <- survfit.gastric.Xelox$time # Extracts the distinct ordered failure times 



loglogss <- log(-log(ss))

logtt <- log(tt)



plot(loglogss ~logtt, main = "Plot of Loglogss against Logtt", pch = 20)



result.lm <- lm(loglogss~logtt) # Fitting a linear model 

abline(result.lm, col = "red") # superimpose the fitted line on the plot 



result.lm





# finding the maximum partial likelihood estimate of ˇ. In R, we may do this by first defining the function 

# l(beta)

plsimple <- function(beta) {

psi <- exp(beta)

result <- log(psi) - log(3*psi + 3) -

log(3*psi + 1) - log(2*psi + 1)

result }

# Finding the Maximum Partial Likelihood estimate

result <- optim(par=0, fn = plsimple, method = "L-BFGS-B",control=list(fnscale = -1),lower = -3, upper = 1)

result$par



# Thus, the m.p.l.e. is beta hat =  -1:326129.





## Homework 5 

# QUESTION 2

library(survival)

setwd("~/Desktop/SPRING 2022/SURVIVAL ANALYSIS/DATA")

geneConfounder <- read.csv("geneConfounder.csv")

result.adj <- coxph(Surv(tt, status)~ trt, data = geneConfounder)

summary(result.adj)





# With treatment as the predictor and tratifyying on genotype

strata_genotype <- factor(geneConfounder$genotype)

result.strata <- coxph(Surv(tt, status) ~ trt + strata(strata_genotype), data = geneConfounder)

summary(result.strata)



# QUESTION 3 



library(survival)

setwd("~/Desktop/SPRING 2022/SURVIVAL ANALYSIS/DATA")

geneConfounder <- read.csv("geneConfounder.csv")

result.trt.gen <- coxph(Surv(tt, status) ~ trt + genotype, data = geneConfounder)

summary(result.trt.gen)





# Simulation of Stratified Survival Data



lambda.mutant.0 <- 0.03

lambda.mutant.1 <- 0.03*0.55

lambda.wt.0 <- 0.03*0.2

lambda.wt.1 <- 0.03*0.2*0.55



set.seed(4321)

tt.control.mutant <- rexp(25, rate=lambda.mutant.0)

tt.treat.mutant <- rexp(125, rate=lambda.mutant.1)

tt.control.wt <- rexp(125, rate=lambda.wt.0)

tt.treat.wt <- rexp(25, rate=lambda.wt.1)

ttAll <- c(tt.control.mutant, tt.treat.mutant, tt.control.wt,

tt.treat.wt)

status <- rep(1, length(ttAll))

genotype <- c(rep("mutant", 150), rep("wt", 150))

trt <- c(rep(0, 25), rep(1, 125), rep(0, 125), rep(1, 25))







# Log Partial Likelihood 



logPartialLikelihood <- function(beta) {

result <- 2*beta - log(2 + 2*exp(beta)) - log(1 + exp(beta))

result

}

beta.vec <- (1:200)/25

logpl <- logPartialLikelihood(beta.vec)

plot(logpl ~ beta.vec, type="l", main = "Partial log Likelihood", col = "red")





# QUESTION 1 Log Partial Likelihood 





logPartialLikelihood <- function(beta) {

result <- 2*beta - log(2 + 2*exp(beta)) - log(1 + exp(beta)) - log(1 + 2*exp(beta)) - log(2 + 3*exp(beta))

result

}

beta.vec <- seq(-1,2, by =0.1)

logpl <- logPartialLikelihood(beta.vec)

plot(logpl ~ beta.vec, type="l", main = "Maximum Likelihood Estimate of Beta", col = "red")

mle <- beta.vec[which((logPartialLikelihood(beta.vec))== max(logPartialLikelihood(beta.vec)))]

abline(v=mle, col = "green")

mle



result.optim.continious <- optim(par = 1.4, fn = logPartialLikelihood, method = "BFGS", 

control = list(fnscale = -1))



result.optim.continious







# QUESTION 1 C 

library(survival)

setwd("~/Desktop/SPRING 2022/SURVIVAL ANALYSIS/DATA")

Clinical_Trial <- read.csv("Clinical_Trial.csv")

attach(Clinical_Trial)

result.ph <- coxph(Surv(time, status) ~ group)

summary(result.ph)



library(survival)

#Loading required package: splines

library(splines)

tt<- c(5,8,7,2,9)

death <- c(1,1,1,1,0)

zz<- c(0,0,1,1,1)

result.ph <- coxph(Surv(tt, death) ~ zz)

summary(result.ph)










library(survival)

#Loading required package: splines

library(splines)

tt <- c(3,4,1,2)

death <- c(1,1,1,0)

zz <- c(1,0,1,0)

result.pl <- coxph(Surv(tt, death) ~ zz)

summary(result.pl)










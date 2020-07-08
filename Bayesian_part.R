
rm(list=ls())
dev.off()

setwd("C:/Users/Emelie/Documents/BayStat/Project")

# project using the olympic dataset. The response label is TotalMedals. Our respons only takes integers and goes from 0 to 104 => poisson distribution    

# full data set
olympic <- read.csv("olympics.csv")

# picking out the important part of the data 
olympic <- olympic[ , c(2 ,3:18)] # ignore the country variable and the last one, it doesn't give us anything
y <- olympic$TotalMedals #our respons, that's given in the project instructions


# Analyzing the data ---------------
#Analyzing the variables and scaling

olympic<- olympics[,-c(1,19)]
y <- olympic$TotalMedals
X <- as.matrix(olympic[,-c(4)])
summary(olympic)
boxplot(X,
        main = "Olympics data set",
        col = "orange",
        border = "brown",
        notch = TRUE
)
y_standar <- (y-mean(y))/sd(y)
X_scale <- scale(X) # put the mean to zero
boxplot(X_scale,
        las=3,
        main = "Olympics data set - scaled",
        col = "orange",
        border = "brown",
        cex.axis=0.75,
        notch = FALSE
)
# correlation matrix 
library(corrplot)
M <- cor(olympic)
corrplot(M, method="shade")

par(mfrow=c(1,2))
boxplot(X_scale,
        las=3,
        main = "Olympics data set - scaled",
        col = "orange",
        border = "brown",
        cex.axis=0.75,
        notch = FALSE
)
corrplot(M, method="shade")

# ------------------------------------------------
# ------------------------------------------------
# Poisson-Gamma model ----------------------------
dev.off()

# Mean and variance about the respons, calc with r 
meanrespons <- mean(y)
varrespons <- var(y)

# gamma parameters calculated from the mean and variance of the responds (TotalMedals)
n <- length(y)
alpha_gamma <- 0.124009
beta_gamma <- 0.0261681
total.samples <-sum(y) 

# Plot the histogram over Total medals
hist(y,prob=T,main="TotalMedals",breaks=50,ylim=c(0,0.37))
lines( dgamma( x=0:100, shape=alpha_gamma, scale=38.2144), col="blue") 

# prior
curve(dgamma(x,shape=alpha_gamma,rate=beta_gamma),0,60,lty=1,ylab=expression(paste(pi,"(",theta,")")),
      xlab=expression(theta),lwd=2,ylim=c(0,2), col=c("blue"))

# Likelihood, m is our lambda, X is the number of medals in n competitions
plot(seq(0,max(y)),dpois(seq(0,max(y)), mean(y)), col=c("red"))

# Posterior
n = nrow(olympic)
alpha_post=alpha_gamma+sum(y)
beta_post=beta_gamma+n
curve(dgamma(x, shape=alpha_post, rate=beta_post),0,100, col=c("green"))

# Plot of them all together
hist(y,prob=T,main="Gamma - Poisson",breaks=50,ylim=c(0,0.6))
curve(dgamma(x,shape=alpha_gamma,rate=beta_gamma),0,100,lty=1,ylab=expression(paste(pi,"(",theta,")")),
      xlab=expression(theta),lwd=2,ylim=c(0,2), col=c("blue"), add=T)
curve(dpois(x, mean(y)),0,100, col=c("red"), add=T)
curve(dgamma(x, shape=alpha_post, rate=beta_post),0,100, add=T, col=c("green"))
legend(80,0.5,c("TotalMedals","Prior","Mle", "Posterior"),cex=0.8,lwd=2,bty="n", col=c("black","blue","red","green"))

#Analyzing the distribution

# Moments of Gamma:
mean.gamma <- function(alpha,beta) alpha/beta
var.gamma <- function(alpha,beta)alpha/beta^2
mode.gamma <- function(alpha, beta) (alpha-1)/beta

mean.prior <- mean.gamma(alpha_gamma,beta_gamma)
var.prior <- var.gamma(alpha_gamma,beta_gamma)
mode.prior <- mode.gamma(alpha_gamma,beta_gamma)

mean.posterior <- mean.gamma(alpha_gamma+total.samples,beta_gamma+n)
var.posterior <- var.gamma(alpha_gamma+total.samples,beta_gamma+n)
mode.posterior <- mode.gamma(alpha_gamma+total.samples, beta_gamma+n)

m.samp <- total.samples/n
mle<-m.samp

resume1<-data.frame(law=c("prior","posterior"), mode=c(mode.prior, mode.posterior),
                    mean=c(mean.prior,mean.posterior),
                    variance=c(var.prior,var.posterior),smv=c(mle,mle)) #check the mode of the prior is a problem

# ------------------------------------------------
# ------------------------------------------------
# Bayesian glm and prediction ----------------

library(dplyr)
library(qpcR)
library(MLmetrics)

# random generates n number between 1 and 203 
set.seed(1234)
n <- 40
x <- sample(1:203, n, replace=T)

# testing on all covariates except the medals, ln and bordapoints
TotalMedals <- olympic$TotalMedals
olympic_2 <- cbind(TotalMedals, olympic[c(6:14)])

# training data
olympic.train <- olympic_2[-c(x), ] 

# pick out test data 
olympic_2 <- olympic[ , c(0 ,6:14)]
olympic.test <- olympic_2 %>% slice(x)

# test responds
y.test <- y[x]

# training responds
y.train <- y[-c(x)] 

library(arm)

baymod = bayesglm(TotalMedals ~., data = olympic.train, family=poisson)
#pred.bay.glm = predict(baymod, newdata = olympic.test, interval = "prediction", level = 0.95)
pred.bay.glm = predict.glm(baymod, newdata = olympic.test, type = "link")
out.bay.glm = cbind(olympic.test$TotalMedals, pred.bay.glm)

# make a list of the predicted values and the real ones 
out.bay.glm
summary(baymod)


par(mfrow=c(1,2))
plot(y.test)
lines(pred.bay.glm, col="blue")
title("Predictions with Bayesglm")
legend(30,50,c("Bayesglm"),cex=0.8,lwd=2,bty="n", col=c("green"))

# plot over the prediction
plot(y.test)
lines(pred.uni$Ypred[1,], col="blue")
lines(pred$Ypred[1,], col="red")
lines(pred.bay.glm, col="green")
title("Predictions with BAS and Bayesglm")
legend(25,50,c("glm predicition","lm prediction","Bayesglm"),cex=0.8,lwd=2,bty="n", col=c("blue","red","green"))



RSS.a <- RSS(baymod)
RSS.a
MSE(pred.bay.glm, olympic.test$TotalMedals)

# Prediction with Zellner's prior----------------------------

library(dplyr)
library(qpcR)
library(MLmetrics)

X <- olympic[ , c(0 ,6:14)]
X <- scale(X)
X <- cbind(X,c(rep(1,203))) #add 1 to get intercept, due to 1*beta_0+x1*beta_1....
y <- olympic$TotalMedals 

# random generates n number between 1 and 203 
set.seed(1234)
n <- 40
x <- sample(1:203, n, replace=T)
x

# pick out test data 
olympic.test <- X[x,] 

# training data
olympic.train <- X[-c(x), ] 

# test responds
y.test <- y[x]

# training responds
y.train <- y[-c(x)] 


betahat=solve(t(olympic.train)%*%(olympic.train),t(olympic.train)%*%y.train)
betahat

c <- 1
betatilde = 0

newbetahat <- 1/(c+1)*(betatilde+c*betahat) #expected value of beta 

newalpha <- newbetahat[length(newbetahat)]

ypred = newalpha + newbetahat[c(1:length(newbetahat)-1)] * olympic.test[ ,c(1:length(newbetahat)-1)]

library(MLmetrics)
MSE(ypred, y.test)

# 2 och 10
# ------------------------------------------------
# ------------------------------------------------
# Model selection Model selection for all except medals, bordapoints and ln in glm: -----------------
# Import library
set.seed(1234)
library(BAS)
library(dplyr)
library(qpcR)
olympic_2 <- olympic[ , c(0 ,6:14)] # take away the log ones, we can use either them or the ordinary ones
olympic_2 <- as.data.frame(scale(olympic_2))

n <- 40 # taking away 20% for testing and using 80% for training
x <- sample(1:203, n, replace=T)
x

# pick out test data 
olympic_test <- olympic_2 %>% slice(x)

# training data
olympic_2 <- olympic_2[-c(x), ] 

# test responds
y.test <- y[x]

# training responds
y.train <- y[-c(x)] 

# Use `bas.lm` to run regression model
cog.BIC.uni = bas.glm(y.train ~ ., data = olympic_2,
                  family=poisson(), betaprior = bic.prior(), modelprior = uniform())
summary(cog.BIC.uni)

cog.coef.uni = coef(cog.BIC.uni)

# We can visualize the coefficients beta1, beta2, bera3, beta4 using the plot function. 
# We use the subset argument to plot only the coefficients of the predictors.

#par(mfrow = c(2, 2), col.lab = "darkgrey", col.axis = "darkgrey", col = "darkgrey")
#plot(cog.coef.uni)

# Find the index of the model with the largest logmarg
best.uni = which.max(cog.BIC.uni$logmarg)

# Retreat the index of variables in the best model, with 0 as the index of the intercept
bestmodel.uni = cog.BIC.uni$which[[best.uni]]

# Create an indicator vector indicating which variables are used in the best model
bestgamma.uni = rep(0, cog.BIC.uni$n.vars) 
bestgamma.uni[bestmodel.uni + 1] = 1 
bestgamma.uni

# Coefficient Estimates Under Reference Prior for Best BIC Model 

# Fit the best BIC model by imposing which variables to be used using the indicators
cog.bestBIC.uni = bas.glm(y.train ~ ., data = olympic_2,
                      family=poisson(), n.models = 1, betaprior = bic.prior(), bestmodel = bestgamma.uni, modelprior = uniform())


# Retreat coefficients information
cog.coef.uni.up = coef(cog.bestBIC.uni)

# Retreat bounds of credible intervals
out.uni = confint(cog.coef.uni.up)[, 1:2]

# Combine results and construct summary table
coef.BIC.uni = cbind(cog.coef.uni.up$postmean, cog.coef.uni.up$postsd, out.uni)
names = c("post mean", "post sd", colnames(out.uni))
colnames(coef.BIC.uni) = names
coef.BIC.uni

pred.uni = predict(cog.bestBIC.uni, newdata = olympic_test, top=1) #calc based on  top 1 model, based on posterior prob
cv.summary.bas(pred.uni$fit, y.test, score="squared-error")

# plot over the prediction
plot(y.test)
lines(pred.uni$Ypred[1,])

MSE(pred.uni$Ypred, y.test)

# Model selection for all except medals, bordapoints and ln in lm: -----------------
dev.off()
# Import library
set.seed(1234)
library(BAS)
library(dplyr)
library(qpcR)

olympic_2 <- olympic[ , c(0 ,6:14)] # take away the log ones, we can use either them or the ordinary ones
olympic_2 <- as.data.frame(scale(olympic_2))

n <- 40 # taking away 20% for testing and using 80% for training
x <- sample(1:203, n, replace=T)
x

# pick out test data 
olympic_test <- olympic_2 %>% slice(x)

# training data
olympic_2 <- olympic_2[-c(x), ] 

# test responds
y.test <- y[x]

# training responds
y.train <- y[-c(x)] 

# Use `bas.lm` to run regression model
cog.BIC = bas.lm(y.train ~ ., data = olympic_2, modelprior = uniform())

summary(cog.BIC.uni)

cog.coef = coef(cog.BIC)

# We can visualize the coefficients beta1, beta2, bera3, beta4 using the plot function. 
# We use the subset argument to plot only the coefficients of the predictors.

#par(mfrow = c(2, 2), col.lab = "darkgrey", col.axis = "darkgrey", col = "darkgrey")
#plot(cog.coef.uni)

# Find the index of the model with the largest logmarg
best = which.max(cog.BIC$logmarg)

# Retreat the index of variables in the best model, with 0 as the index of the intercept
bestmodel = cog.BIC$which[[best]]
bestmodel

# Create an indicator vector indicating which variables are used in the best model
bestgamma = rep(0, cog.BIC$n.vars) 
bestgamma[bestmodel + 1] = 1  
bestgamma


# Fit the best BIC model by imposing which variables to be used using the indicators
cog.bestBIC = bas.lm(y.train ~ ., data = olympic_2, n.models = 1, bestmodel = bestgamma, modelprior = uniform())


cog.coef = coef(cog.bestBIC)
out = confint(cog.coef)[, 1:2]

# Combine results and construct summary table
coef.BIC = cbind(cog.coef$postmean, cog.coef$postsd, out)
names = c("post mean", "post sd", colnames(out))
colnames(coef.BIC) = names
coef.BIC

pred = predict(cog.bestBIC, newdata = olympic_test, top=1)
cv.summary.bas(pred$Ypred, y.test, score="squared-error")

# plot over the prediction
plot(y.test)
lines(pred.uni$Ypred[1,], col="blue")
lines(pred$Ypred[1,], col="red")
title("Predictions with BAS")
legend(30,50,c("glm predicition","lm prediction"),cex=0.8,lwd=2,bty="n", col=c("blue","red"))


# error
MSE <- MSE(pred$Ypred, y.test)
MSE
# ------------------------------------------------
# ------------------------------------------------
# JAGS ZIPois ----------------------------
library(rjags)
library(coda)
n = 203 # number of observations

ZIPois_string <- "model{
  # likelihood
    for (i in 1:n) {
    y[i]~dpois(mu[i])
    mu[i] <- lambda*z[i] + 0.00001  ## 0 is not admitted for poisson
    z[i]~dbern(1-q)
    }
    #predictive 
    yp~dpois(mup)
    mup<-lambda*zp+ 0.00001 
    zp~dbern(1-q)
    # prior 
    lambda ~ dgamma(200,200)
    q ~ dbeta(50,1)
}"

jagsZPois<- jags.model(textConnection(ZIPois_string),
                       data = list("y" = y,"n" = n),
                       n.chains=2,
                       n.adapt=300)

update(jagsZPois, 5000)

outputmcmcZPois=coda.samples(jagsZPois,
                             c('q',"lambda","yp"),
                             n.iter=5000, progress.bar="none")

ypouts=coda.samples(jagsZPois,
                             c("yp"),
                             n.iter=203, progress.bar="none")

plot(outputmcmcZPois)
su=summary(outputmcmcZPois)

ypredictive=as.matrix(outputmcmcZPois[,c("yp")])
lmest=su$statistics[,1][1]
qest=su$statistics[,1][2]
xx=seq(0,104,1)
ddd=dpois(xx,lmest)*(1-qest)+qest*(xx==0)
hist(y,seq(-0.5,max(y)+.5,1),probability = TRUE,main='ZI-Poisson model')
lines(xx,ddd,col="red")

lines(xx,dpois(xx,mean(y)),col="blue")
legend("topright",legend=c("ZI Poisson", "plug-in"),col=c("red", "blue"), lty=1:1, cex=0.8)


# JAGS for beta -------------------

library(rjags)
library(coda)
olympic1<- olympic[,-c(1,2,3,5,15,16,17)]

ZIPois_string <- "model{
  # likelihood
    for (i in 1:n) {
    y[i]~dpois(mu[i])
    mu[i] <- lambda[i]*z[i] + 0.00001  ## 0 is not admitted for poisson
    z[i]~dbern(1-q)
    lambda[i]=inprod(beta[],X[i,]) ##For every covariant
    
    }
    #predictive 
    yp~dpois(mup)
    mup<-lambdap*zp+ 0.00001 
    zp~dbern(1-q)
    
    # prior 
    lambdap ~ dgamma(200,200)
    q ~ dbeta(50,1)
    
    # prior distributions
    beta[1:P] ~ dmnorm( mu.beta[], prior.T[,] )
    tau    ~ dgamma(50,1 )
    
    # prior
    c2 <- n
    
    # prior means
    for (j in 1:P){ mu.beta[j] <- 0.0 }
    
    # calculation of xtx
    for (i in 1:P){ for (j in 1:P){
       inverse.V[i,j] <- inprod( X[,i] , X[,j] )}}
    
    for(i in 1:P){ for (j in 1:P){
      prior.T[i,j] <- inverse.V[i,j] * tau /c2  }}
}"

n <- nrow(olympic1)
total.samples <-y
response=total.samples
X=model.matrix(~., data=olympic1)
jagsZPois<- jagsUI::jags(textConnection(ZIPois_string),
                         data = list("y" =response,"n" =length(response) , "X"=X, "P"=ncol(X)),
                         n.chains=3,
                         n.adapt=4000,
                         inits=NULL,
                         n.burnin = 40000,
                         n.iter=100000,
                         n.thin = 200,
                         parameters.to.save = c('q','yp'))

col_names <- c("intercept",names(olympic1))
#plot(jagsZPois)
summary(jagsZPois)


variable_selection <- data.frame(names = col_names,
                                 mean = jagsZPois$mean$beta,
                                 overlap0 = jagsZPois$overlap0$beta)

variable_selection <- data.frame(names = col_names,
                                 mean = jagsZPois$mean$q,
                                 overlap0 = jagsZPois$overlap0$q)

jagsZPois



# ------------------------------------------------
# ------------------------------------------------
## GST Rheotaxis Experiment 2: Barbel Trials

# Author: Myfanwy Johnston
# Date: April 5, 2016
# Note: all models were run with multiple chains the first time around

# Lines 16 - 37: Data loading and preparation
# Lines 40 - 147: Experimental Models for Barbel Trials
# Lines 148 - 170: Model comparison for Barbel Trials
# Lines 172 - 403: Data and experimental models for before/after barbelectomy comparison
# Lines 403 - end: Model comparison, exploration, prediction simulations, and visualizations for before/after barbelectomy comparison

# Required packages:
library(rethinking) # version 1.55
library(rstan) # version 2.9.0-3
library(beepr)
## Load and prep data ##
#------------------------------------------------------------------------
d <- read.csv("../data_tidy/gst_rheo_all.csv", stringsAsFactors = FALSE, header = TRUE)
## Filter data to include trials without barbels
# d3 <- dplyr::filter(d, Type == "without") #uncomment for dplyr method
d3 <- d[d$Type == "without", ]

## Check for zeros
y <- d3$PropPos
range(y)
sum( dbeta2( y , 0.5 , 10 , log=TRUE ) )
dens(d3$PropPos) # See shape of raw data; highly skewed, not as bi-modal as Exp 1; much more positive rheotaxis displayed in this group

# Create dummy variables for treatments Light and Dark
d3$tL <- ifelse(d3$Treatment == "L", 1, 0)
d3$tD <- ifelse(d3$Treatment == "D", 1, 0)


# Make data list
dlist <- list(y = y, fish_id = coerce_index(d3$FishID), treatment = coerce_index(d3$Treatment),
              light = d3$tL, dark = d3$tD)

## Begin Experiment 2 Models: Barbelectomy Trials Only ##
#------------------------------------------------------------------------
m0b <- map2stan(
  alist(
    y ~ dbeta2(p,theta),
    logit(p) <- a,
    a ~ dnorm(0, 1), 
    theta ~ dcauchy(0, 1)
  ),
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta = 1),
  sample=TRUE , warmup=1000 , iter=1e4 , 
  cores=2 , chains=1 )

# Fixed effects model: fish
#------------------------------------------------------------------------
m1b <- map2stan(
  alist(
    y ~ dbeta2(p,theta),
    logit(p) <- a_fish[fish_id],
    a_fish[fish_id] ~ dnorm(0,1),
    theta ~ dcauchy(0,1)
  ),
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2 )

# Varying intercepts: fish
#------------------------------------------------------------------------
m2b <- map2stan(
  alist(
    y ~ dbeta2(p,theta),
    logit(p) <- a_fish[fish_id],
    a_fish[fish_id] ~ dnorm(a, sigma_fish),
    a ~ dnorm(0, 1),
    sigma_fish ~ dcauchy(0,1),
    theta ~ dcauchy(0,1)
  ),
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2 )

# varying intercepts for fish by (fixed) treatment effect
#------------------------------------------------------------------------
m3b <- map2stan(
  alist(
    y ~ dbeta2( p, theta ),
    logit(p) <- a + a_fish[fish_id] + b_treatment*treatment ,
    a_fish[fish_id] ~ dnorm(0, sigma_fish),
    a ~ dnorm(0, 1),
    b_treatment ~ dnorm(0, 1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dcauchy(0, 1)
  ),
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000 , iter=1e4, cores=2)

# varying intercepts on fish and treatment
#------------------------------------------------------------------------------------
m4b <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    # linear model
    logit(p) <- a + a_fish[fish_id] + b_treat[treatment],
    # adaptive priors
    a_fish[fish_id] ~ dnorm(0,sigma_fish),
    b_treat[treatment] ~ dnorm(0, sigma_treat),
    # fixed priors
    a ~ dnorm(0, 1),
    theta ~ dexp(1),
    sigma_fish ~ dcauchy(0,1),
    sigma_treat ~ dcauchy(0,1)
  ),
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000 , iter=1e4 , cores=2 )

#------------------------------------------------------------------------
# Varying intercepts and slopes by fish and treatment
#-----------------------------
m1bNC <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    
    # linear model
    logit(p) <- (b_light + bl_fish[fish_id])*light + (b_dark + bd_fish[fish_id])*dark,
    
    # adaptive NON-CENTERED priors 
    c(bl_fish, bd_fish)[fish_id] ~ dmvnormNC(sigma_fish, Rho_fish),
    
    # fixed priors
    c(b_light, b_dark) ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dexp(1),
    Rho_fish ~ dlkjcorr(2)
  ),
  # data
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2)


# ------------------------------------------
## Begin Model comparison for Barbel Trials
# -----------------------------------------
gst_models <- compare(m0b, m1b, m2b, m3b, m4b, m1bNC)
gst_models # smaller WAIC indicates better out-of-sample estimated deviance. pWAIC = estimated effective number of parameters. m3b has the majority of the weight, but barely - basically, all the models that include treatment as a fixed effect are in the top 3.

plot(gst_models) 


# compute predicted proportions for each case in data
pred3b <- sim(m3b)

# to summarize
p3b <- apply( pred8b , 2, mean )


# Plot raw data against simulated predicted points from one of the top models
library(ggplot2)
par(mfrow = c(1,1))

m3bcompare <- data.frame(observed = d3$PropPos, predicted = p3b)
ggplot(m3bcompare, aes(x = observed, y = predicted)) + geom_point() + geom_smooth() 
# model does not do as good a job predicting.


#-------------------------------------------------------------------------------------------------
## Begin Comparison of Barbels/Non Barbel Velocity Trials
#------------------------------------------------------------------------

## Load Data ##
## Filter data to include only velocity trials (Treatments L and D), with or without barbels
# d2b <- dplyr::filter(d, Treatment == "L" | Treatment == "D")
library(dplyr)
d2b <- d[d$Treatment == "L" | d$Treatment == "D", ]
length(unique(d2b$FishID))

with <- filter(d, Type == "with")
without <- filter(d, Type == "without")

## Determine which fish were in both experiments 1 & 2
Fishwith <- as.data.frame(sort(unique(with$FishID))) 
Fishwithout <- as.data.frame(sort(unique(without$FishID)))
names(Fishwith) <- "FishID"
names(Fishwithout) <- "FishID"

reps <- semi_join(Fishwithout, Fishwith) #join for all values of x that are in y (Fishwithout that are in Fishwith)
reps
d2b <- filter(d2b, FishID %in% reps$FishID)

## Create 0/1 indicator for barbels present (1) or absent (0)
d2b$Barbels <- ifelse(d2b$Type == "with", 1, 0)
# Create dummy variables for treatment
d2b$tL <- ifelse(d2b$Treatment == "L", 1, 0)
d2b$tD <- ifelse(d2b$Treatment == "D", 1, 0)

## fix the one zero element
y <- d2b$PropPos
range(y)
y[y==0] <- 0.0001
sum( dbeta2( y , 0.5 , 10 , log=TRUE ) )
head(d2b)

dlistB <- list(y = y, fish_id = coerce_index(d2b$FishID), treatment = coerce_index(d2b$Treatment),
               light = d2b$tL, dark = d2b$tD, barbels = d2b$Barbels)

## Barbel Models (w/wo comparison) ##
#------------------------------------------------------------------------

# fishID only
mB0 <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    
    # linear model
    logit(p) <- a_fish[fish_id],
    
    a_fish[fish_id] ~ dnorm(a, sigma_fish),
    a ~ dnorm(0, 1),
    sigma_fish ~ dcauchy(0,1),
    theta ~ dcauchy(0,1)
  ),
  # data
  data = dlistB,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2)

#------------------







# varying intercepts with fish, fixed treatment effect
mB1 <- map2stan(
    alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    
    # linear model
    logit(p) <- a_fish[fish_id] + b_treatment*treatment ,
    
    # adaptive NON-CENTERED priors 
    a_fish[fish_id] ~ dmvnormNC(sigma_fish, Rho_fish),
    
    # fixed priors
    b_treatment ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dexp(1),
    Rho_fish ~ dlkjcorr(2)
  ),
  # data
  data = dlistB,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2)

#------------------------------------------------------------------------

#  allowing treatment slope to vary by fish
mB2 <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    
    # linear model
    logit(p) <- a_fish[fish_id] + (b_treatment + b_fish[fish_id])*treatment ,
    
    # adaptive NON-CENTERED priors 
    c(a_fish, b_fish)[fish_id] ~ dmvnormNC(sigma_fish, Rho_fish),
    
    # fixed priors
    b_treatment ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dexp(1),
    Rho_fish ~ dlkjcorr(2)
  ),
  # data
  data = dlistB,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2)
#------------------------------------------------------------------------

#  varying intercepts on fish, individual treatment fixed effects
mB3 <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    
    # linear model
    logit(p) <- a_fish[fish_id] + b_light*light + b_dark*dark ,
    
    # adaptive NON-CENTERED priors 
    a_fish[fish_id] ~ dnorm(0, sigma_fish),
    
    # fixed priors
    c(b_light, b_dark) ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dexp(1)
  ),
  # data
  data = dlistB,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2)

#------------------------------------------------------------------------
# Varying intercepts and slopes by fish and treatment
#-----------------------------
mB4 <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    
    # linear model
    logit(p) <- a_fish[fish_id] + (b_light+ bl_fish[fish_id])*light + (b_dark + bd_fish[fish_id])*dark ,
    # adaptive NON-CENTERED priors 
    c(a_fish, bl_fish, bd_fish)[fish_id] ~ dmvnormNC(sigma_fish, Rho_fish),
    
    # fixed priors
    c(b_light, b_dark) ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dexp(1),
    Rho_fish ~ dlkjcorr(2)
  ),
  # data
  data = dlistB,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2)

## Adding barbels as a fixed effect
mB5 <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    
    # linear model
    logit(p) <- a_fish[fish_id] + b_barbels*barbels ,
    
    # adaptive NON-CENTERED priors 
    a_fish[fish_id] ~ dmvnormNC(sigma_fish, Rho_fish),
    
    # fixed priors
    b_barbels ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dexp(1),
    Rho_fish ~ dlkjcorr(2)
  ),
  # data
  data = dlistB,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2)

#  Allowing barbels and treatment as fixed effects

mB6 <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    
    # linear model
    logit(p) <- a_fish[fish_id] + b_barbels*barbels + b_treatment*treatment,
    
    # adaptive NON-CENTERED priors 
    a_fish[fish_id] ~ dmvnormNC(sigma_fish, Rho_fish),
    
    # fixed priors
    c(b_barbels, b_treatment) ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dexp(1),
    Rho_fish ~ dlkjcorr(2)
  ),
  # data
  data = dlistB,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2)

#  Allowing barbel effect slope to vary with treatment

mB7 <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    
    # linear model
    logit(p) <- a_fish[fish_id] + (b_barbels + b_treatment[treatment])*barbels ,
    
    # adaptive NON-CENTERED priors 
    a_fish[fish_id] ~ dmvnormNC(sigma_fish, Rho_fish),
    b_treatment[treatment] ~ dmvnormNC(sigma_treatment, Rho_treatment),
    # fixed priors
    b_barbels ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    c(sigma_fish, sigma_treatment) ~ dexp(1),
    c(Rho_fish, Rho_treatment) ~ dlkjcorr(2)
  ),
  # data
  data = dlistB,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2)

#  Allowing barbel effect slope to vary with fish
mB8 <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    
    # linear model
    logit(p) <- a_fish[fish_id] + (b_barbels + b_fish[fish_id])*barbels ,
    
    # adaptive NON-CENTERED priors 
    c(a_fish, b_fish)[fish_id] ~ dmvnormNC(sigma_fish, Rho_fish),
    
    # fixed priors
    b_barbels ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dexp(1),
    Rho_fish ~ dlkjcorr(2)
  ),
  # data
  data = dlistB,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000, iter=1e4, cores=2)
beep(5)
############################################################
## Begin Barbel Trial/non-Barbel Trial Model Comparison
############################################################
barbelmodels <- compare(mB0, mB1, mB2, mB3, mB4, mB5, mB6, mB7, mB8)
barbelmodels


p <- precis(mB7, prob = 0.95)
p <- data.frame(p@output)
p2 <- sapply(p[ , 1:4], rethinking::logistic)
precis(mB7, prob = 0.95)
p2
dbfish <- filter(d, FishID %in% reps$FishID) # get fish in both experiments again

# calculate mean & sd of response variable for Exp 1 for those fish
dbfish %>% 
  filter(Treatment == "D" | Treatment == "L", Type == "with") %>% 
  summarise(mean = mean(PropPos), sd = sd(PropPos))

# calculate mean & sd of response variable for Exp 2 for those fish
dbfish %>% 
  filter(Treatment == "D" | Treatment == "L", Type == "without") %>% 
  summarise(mean = mean(PropPos), sd = sd(PropPos))

# Plot Predictions Across Treatments
gst_ensemble <- ensemble(mB5) 
# dummy data for predictions across treatments
d.pred <- data.frame( Treatment <- c("L", "D"))

# build prediction ensemble of top 3 models
gst_pensemble <- ensemble(mB5, mB6, mB8)

# summarize
pred.ps <- apply(gst_pensemble$sim, 2, mean)
pred.ps.PI <- apply(gst_pensemble$sim, 2, PI)

plot(0, 0, type = "n", xlab = "Treatment",
     ylab = "Proportion Positively Oriented" , ylim = c(0,1), xaxt = "n" ,
     xlim = c(1, 2))
axis(1, at = 1:2, labels = c("L", "D"))

p <- by(d2b$PropPos, 
        list(d2b$Treatment, d2b$FishID), mean ) # calculate the mean of PropPos for each combination of Treatment and FishID

for( FishID in 1:11)
  lines(1:2, as.vector(p[,FishID]), col = rangi2)

# Add prediction mean lines
lines(1:44, pred.ps) # shows the average predicted proportion of positive rheotaxis, across treatments (sim)
shade(pred.ps.PI, 1:44) #95% percentile interval of the mean prediction
title(main = "Simulated Predictions by Treatment, Averaging over the Posterior Distribution Top Model")

## GST Rheotaxis Experiment 1: Velocity vs. Optimotor Belt Trials

# Author: Myfanwy Johnston
# Date: September 30, 2016
# Note: all models were run with multiple chains the first time around

# Required packages:
library(rethinking) # version 1.55
library(rstan) # version 2.9.0-3
library(beepr) # to alert you when models are finished running
# library(dplyr) # Alternatives to all but one dplyr method used are provided


# Explanation of Experiment 1 ---------------------------------------------

# These data come from a total of 96 trials of juvenile green sturgeon rheotaxis (directionality in the face of perceived current). 23 fish underwent four different treatments under two different conditions. These were:
#   
#   Condition 1: velocity present, optimotor belt absent (48 of 96 trials)
#   Condition 2: optimotor belt present, velocity not present (48 of 96 trials)
# 
# Treatment A: optimotor belt placed above experimental tank, no velocity (24 trials)
# Treatment B: optimotor belt placed below experimental tank, no velocity (24 trials)
# Treatment L: velocity during the day (light) (24 trials)
# Treatment D: velocity during the night (dark) (24 trials)
# 
# The goal was to determine which treatment or condition had the greatest influence on proportion of time spent positively oriented (orienting nose-first into the current or visual "current").


# Load Data ---------------------------------------------------------------
d <- read.csv("../data_tidy/gst_rheo_all.csv", stringsAsFactors = FALSE, header = TRUE)
head(d)

# Prep Data for Exp 1 Models
d2 <- d[d$Type == "with", ] # filter trials from Experiment 1, where barbels were intact. 96 rows.
d2$Velocity <- ifelse(d2$Treatment == "L" | d2$Treatment == "D", 1, 0) # Create dummy variable for velocity

## There is one zero in the data; after running all models without this trial, we adjusted its value with the following code and ran models again; we found no difference in outcome.

y <- d2$PropPos
y[y==0] <- 0.0001 
sum( dbeta2( y , 0.5 , 10 , log=TRUE ) ) #double-check that there are no zeros in the data

## Assign dummy variables to columns in the dataframe
d2$tA <- ifelse(d2$Treatment == "A", 1, 0)
d2$tB <- ifelse(d2$Treatment == "B", 1, 0)
d2$tL <- ifelse(d2$Treatment == "L", 1, 0)
d2$tD <- ifelse(d2$Treatment == "D", 1, 0)

dens(d2$PropPos) # See shape of raw data - highly bimodal

# compile relevant predictors/variables into list for running models
dlist <- list(y = y, fish_id = coerce_index(d2$FishID), vel = d2$Velocity, treatment = coerce_index(d2$Treatment),
              above = d2$tA, below = d2$tB, light = d2$tL, dark = d2$tD)

# Note: Models below are structured with single chains.  To check stationarity of HMC trace plots, run models with multiple chains (chains = n) in the constraints list first.


# Begin Experiment 1 Models -----------------------------------------------

## Intercept Model ##
#------------------------------------------------------------------------
m0 <- map2stan(
  alist(
    y ~ dbeta2(p,theta),
    logit(p) <- a,
    a ~ dnorm(0, 1),
    theta ~ dcauchy(0,1)
  ),
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta = 1),
  sample=TRUE , warmup=1000 , iter=1e4 , 
  cores=2, chains = 1 )


# Fixed effects model: fish
#------------------------------------------------------------------------
m1 <- map2stan(
  alist(
    y ~ dbeta2(p, theta),
    logit(p) <- a_fish[fish_id],
    a_fish[fish_id] ~ dnorm(0,1),
    theta ~ dcauchy(0,1)
  ),
  data=dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000 , iter=1e4, cores=2)

# varying intercepts: fish
#------------------------------------------------------------------------
m2 <- map2stan(
  alist(
    y ~ dbeta2(p,theta),
    logit(p) <- a_fish[fish_id],
    a_fish[fish_id] ~ dnorm(a, sigma_fish),
    a ~ dnorm(0, 1),
    sigma_fish ~ dcauchy(0,1),
    theta ~ dcauchy(0,1)
  ),
  data=dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000 , iter=1e4 , cores=2 )

# varying intercepts for fish by (fixed) treatment effect (single estimate for treatment)
#------------------------------------------------------------------------
m3 <- map2stan(
  alist(
    y ~ dbeta2( p, theta ),
    logit(p) <- a + a_fish[fish_id] + b_treatment*treatment ,
    a_fish[fish_id] ~ dnorm(0, sigma_fish),
    a ~ dnorm(0, 1),
    b_treatment ~ dnorm(0, 1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dcauchy(0,1)
  ),
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000 , iter=1e4 , cores=2 )

# Varying intercepts on fish, fixed (individual) treatment effects:
#------------------------------------------------------------------------
m4 <- map2stan(
  alist(
    y ~ dbeta2(p,theta),
    logit(p) <- a_fish[fish_id] + b_above*above + b_below*below + b_light*light + b_dark*dark, # can break them out separately like this so long as you've ommitted the global mean; if not, you have to have k-1 dummy variables.
    a_fish[fish_id] ~ dnorm(0, sigma_fish),
    c(b_above, b_below, b_light, b_dark) ~ dnorm(0, 1),
    sigma_fish ~ dcauchy(0,1),
    theta ~ dcauchy(0,1)
  ),
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta = 1),
  sample=TRUE , warmup=1000 , iter=1e4 , 
  cores=2, chains = 1 )

# varying intercepts on fish and treatment (treatment as a mean, not individual)
#------------------------------------------------------------------------------------
m4.5 <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    # linear model
    logit(p) <- a + a_fish[fish_id] + a_treat[treatment],
    # adaptive priors
    a_fish[fish_id] ~ dnorm(0,sigma_fish),
    a_treat[treatment] ~ dnorm(0, sigma_treat),
    # fixed priors
    a ~ dnorm(0, 1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dcauchy(0,1),
    sigma_treat ~ dcauchy(0,1)
  ),
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000 , iter=1e4 , cores=2 )


#--------------------------------------------------
# fish as varying intercepts, velocity as fixed effect
m5 <- map2stan(
  alist(
    y ~ dbeta2( p, theta ),
    logit(p) <- a + a_fish[fish_id] + b_velocity*vel,
    a_fish[fish_id] ~ dnorm(0, sigma_fish),
    a ~ dnorm(0, 1),
    b_velocity ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dcauchy(0,1)
  ),
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000 , iter=1e4 , cores=2 )
beep(0)

## interestingly, if you just compare(m0, m1, m2, m3, m4, m4.5) up to this point, m4 is top-weighted. So even when the average effect of treatment (in general) is allowed to vary, individual fixed effects on treatment is better.  Suggests real individual treatment effect.  If you add m5, it splits the weight between the models that include treatment and/or velocity.  Between m4 and m5, no clear winner.  Which means that velocity might be as important as individual treatment.

#------------------------------------------------------------------------------------
# varying intercepts: fish and velocity
m6 <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    # linear model
    logit(p) <- a + a_fish[fish_id] + a_velocity[vel],
    # adaptive priors
    a_fish[fish_id] ~ dnorm(0,sigma_fish),
    a_velocity[vel] ~ dnorm(0, sigma_vel),
    # fixed priors
    a ~ dnorm(0, 1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dcauchy(0,1),
    sigma_vel ~ dcauchy(0,1)
  ),
  data = dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000 , iter=1e4 , cores=2 )
beep(2)
#--------------------------------------------------
# varying intercepts (fish) + fixed effect interaction model (fish X velocity)
m7 <- map2stan(
  alist(
    y ~ dbeta2( p, theta ),
    logit(p) <- a + a_fish[fish_id] + b_velocity*vel + bfXv*fish_id*vel,
    a_fish[fish_id] ~ dnorm(0, sigma_fish),
    a ~ dnorm(0, 1),
    b_velocity ~ dnorm(0, 1),
    bfXv ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dcauchy(0,1)
  ),
  data=dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1) , warmup=1000 , iter=1e4 , cores=2)

#---------
# varying intercepts (fish) + fixed treatment correlation model (says that fish vary overall in positive orientation), but
# everything else is constant across fish
#----------------------------------------------------------------------------
m8 <- map2stan(
  alist(
    #liklihood
    y ~ dbeta2( p, theta ),
    #linear model
    logit(p) <- a_fish[fish_id] + (b_above*above + b_below*below)*a_velocity*vel + 
      (b_light*light + b_dark*dark)*b_velocity*vel,
    #adaptive priors
    a_fish[fish_id] ~ dnorm(0, sigma_fish),
    
    #fixed priors
    c(b_above, b_below, a_velocity, b_light, b_dark, b_velocity) ~ dnorm(0, 1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dcauchy(0,1)
  ),
  data=dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1) , warmup=1000 , iter=1e4 , cores=2, chains = 3)

# Mixed model with varying intercepts on fish, varying slopes on fish/velocity
#------------------------------------------------------------------------
m1NC <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    # linear model
    logit(p) <- a + a_fish[fish_id] + (b_velocity + b_fish)*vel,
    # adaptive NON-CENTERED priors 
    c(a_fish, b_fish)[fish_id] ~ dmvnormNC(sigma_fish, Rho_fish),
    # fixed priors
    c(a, b_velocity) ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dexp(1),
    Rho_fish ~ dlkjcorr(2)
  ),
  # data
  data=dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1) , warmup=1000 , iter=1e4 ,cores=2 )

#------------------------------------------------------------------------
# Varying intercepts and slopes by fish and treatment
#-----------------------------
m2NC <- map2stan(
  alist(
    #likelihood
    y ~ dbeta2( p, theta ),
    
    # linear model
    logit(p) <- (b_above + ba_fish[fish_id])*above + (b_below + bb_fish[fish_id])*below + 
      (b_light + bl_fish[fish_id])*light + (b_dark + bd_fish[fish_id])*dark,
    
    # adaptive NON-CENTERED priors 
    c(ba_fish, bb_fish, bl_fish, bd_fish)[fish_id] ~ dmvnormNC(sigma_fish, Rho_fish),
    # fixed priors
    c(b_above, b_below, b_light, b_dark) ~ dnorm(0,1),
    theta ~ dcauchy(0,1),
    sigma_fish ~ dexp(1),
    Rho_fish ~ dlkjcorr(2)
  ),
  # data
  data=dlist,
  constraints=list(theta="lower=0"),
  start=list(theta=1), warmup=1000 , iter=1e4 , cores=2 )
beep(0)
#---------------------------------
# Begin Model Comparison --------------------------------------------------
# -------------------------
gst_models <- compare(m0, m1, m2, m3, m4, m4.5, m5, m6, m7, m8, m1NC, m2NC)
gst_models # smaller WAIC indicates better out-of-sample estimated deviance. pWAIC = estimated effective number of parameters. 
plot(gst_models) # m2NC has by far the smallest out of sample estimated deviance.

precis(m2NC, prob = 0.95)
p <- precis(m2NC, prob = 0.95) # to convert this from log-odds to proportion (as in the manuscript text), take the logistic of coefficient estimates

plot(precis(m2NC, pars = c("b_above", "b_below", "b_light", "b_dark"), depth = 2)) #intercepts have much wider margins for treatments light & dark

post <- extract.samples(m2NC)
dens(post$sigma_fish, xlab = "sigma", xlim = c(0, 4)) # we can think of the sigma parameter for fish as a crude measure of the cluster's relevance for explaining the variation in the outcome - we observe bimodality in its density just as we observe bimodality in the raw data.


#---------------------------------
# Begin Model Exploration
# -------------------------


# Generate Predictions by Indiviudal Fish from Top Model
#-----------------------------------------------------------
# create dummy data for predictions across treatments
d.pred <- data.frame(
  Treatment = c("A", "B", "L", "D"))

# build prediction ensemble of top model
gst_pensemble <- ensemble(m2NC)

# summarize
pred.ps <- apply(gst_pensemble$sim, 2, mean) # items 48-72 = Treatment "L"
pred.ps.PI <- apply(gst_pensemble$sim, 2, PI)

# create empty plot
plot(0, 0, type = "n", xlab = "Treatment",
     ylab = "Proportion Positively Oriented" , ylim = c(0,1), xaxt = "n" ,
     xlim = c(1, 4))
axis(1, at = 1:4, labels = c("A", "B", "L", "D"))

# calculate the mean of PropPos for each combination of Treatment and FishID
p <- by(d2$PropPos, 
        list(d2$Treatment, d2$FishID), mean ) 

for( FishID in 1:23)
  lines(1:4, as.vector(p[,FishID]), col = rangi2) #plot raw data

# Add prediction mean lines
lines(1:96, pred.ps) # shows the average predicted proportion of positive rheotaxis, across treatments (sim)
shade(pred.ps.PI, 1:96) #95% percentile interval of the mean prediction
title(main = "Observed Data and Predicted Means for Positive Rheotaxis Across Treatments")

# there is still a lot of variation among individuals, although the model is supposed to account for this with the varying intercepts, but averages seem particularly bad at the "L" treatment, as most lines tick upwards for the "L,"; this is likely because no fish were near the actual mean, in the raw observations.  There is more uncertainty surrounding treatment A than B, as well.


str(post) # samples for the parameters are held in matrices. Treatments are in 1-dimensional matrices with 9000 rows.  For the fishxtreatment parameters, they're stored in matrices with 9000 rows each and 23 columns each.  Each column is a parameter, and each row is a sample from the posterior distribution.  So to plot the density for fish 1, we ask for all the rows from column 1:

dens(post$ba_fish[, 1])

# Plot individual fish

fish <- 3

d.pred <- list(
  PropPos = rep(0, 4),
  treatment = c("above", "below", "light", "dark"),
  fish_id = rep(fish, 4)
)

link.m2NC <- link(m2NC)



# compute predicted proportions for each case in data
#------------------------------------------------------
pred <- sim(m2NC)
str(pred)

# to summarize
pmv <- apply( pred , 2 , mean )
# create empty plot
plot(0, 0, type = "n", xlab = "Treatment",
     ylab = "Proportion Positively Oriented" , ylim = c(0,1), xaxt = "n" ,
     xlim = c(1, 4))
axis(1, at = 1:4, labels = c("A", "B", "L", "D"))


d.pred <- data.frame("A" = pmv[1:24], "B" = pmv[25:48], "L" = pmv[49:72], "D" = pmv[73:96])
str(d.pred)

points(d.pred$A, col = "red")


# Plot raw data against simulated predicted points
library(ggplot2)
par(mfrow = c(1,1))

mcompare <- data.frame(observed = d2$PropPos, predicted = pmv)
ggplot(mcompare, aes(x = observed, y = predicted)) + geom_point() + geom_smooth(method = "lm") + ggtitle("m2NC: y ~ [fish_id] + b1[fish_id]*Above, etc")

# Draw Posterior density of m2NC
dens(pmv, lty = 2, col = "brown", lwd = 1.5) # posterior predictions
dens(d2$PropPos, add = TRUE, lwd = 1.5) # raw data
title("Density of Observed Data vs. Posterior Predictions of Top Model")

# Plot Parameter Shapes
#----------------------------------------------------------------
library(dplyr)
# Extract Samples from Posterior
post <- extract.samples(m2NC)
str(post)



# Plot separately
dens(treatl, lty = 2, xlim = c(0, 1), ylim = c(0, 3.5), lwd = 2, col = "dodgerblue", adj = 0.5)
dens(treata, col = "red", ylim = c(0, 4), xlim = c(0, 1))
dens(treatb, col = "brown", add = TRUE)
dens(treatl, col = "dodgerblue", add = TRUE)
dens(treatd, col = "blue", add = TRUE)
# we can see the distinct distributions for the two conditions (velocity/nonvelocity).  Very distinct bi-modality; the model is capturing this aspect of the data very well.


# extract simulated samples from posterior
sim_A <- rnorm(1000, mean = post$b_above, sd = post$sigma_fish)
sim_B <- rnorm(1000, post$b_below, post$sigma_fish)
sim_L <- rnorm(1000, post$b_light, post$sigma_fish)
sim_D <- rnorm(1000, post$b_dark, post$sigma_fish)

# Calculate mean/sd/range
mean(c(logistic(sim_A), logistic(sim_B)))
sd(c(logistic(sim_A), logistic(sim_B)))
mean(c(logistic(sim_L), logistic(sim_D)))
sd(c(logistic(sim_L), logistic(sim_D)))
range(logistic(sim_B))
range(logistic(sim_D))
range(logistic(sim_L))

# transform to a proportion and visualize
dens(logistic(sim_D), col = "Red", xlab = "Proportion Positive Rheotaxis", ylim = c(0, 12), xlim = c(0, 1))
dens(logistic(sim_L), add = TRUE, col = "blue")
dens(logistic(sim_A), add = TRUE)
dens(logistic(sim_B), add = TRUE, col = "brown")
title('Simulated Proportion Time Positively Oriented by Treatment')
title(sub = 'black = Above brown = Below red = Dark blue = Light')

sim_A <- rnorm(24, mean = post$b_above, sd = post$sigma_fish)
dens(logistic(sim_A))
library(dplyr)
tA <- d2 %>% 
  filter(Treatment == "A") %>% 
  select(PropPos)
dens(tA, lty = 2, add = TRUE)



# Calculate Contrast Values -----------------------------------------------
diffLD <- logistic(post$b_dark) - logistic(post$b_light)
dens(diffLD)
mean(diffLD)

diffAB <- logistic(post$b_below) - logistic(post$b_above)
mean(diffAB)

dens(diffLD, col = "blue")
dens(diffAB, add = TRUE, col = "brown")

plot(diffLD, col = col.alpha("blue"))
points(diffAB, col = col.alpha("brown"), add = TRUE)

diffAL <- logistic(post$b_above) - logistic(post$b_light)
mean(diffAL)
diffBL <- logistic(post$b_below) - logistic(post$b_light)
mean(diffBL)
plot(diffAL, col = col.alpha("yellow"))
points(diffBL, col = col.alpha("green"))

diffAD <- logistic(post$b_above) - logistic(post$b_dark)
mean(diffAD)
diffBD <- logistic(post$b_below) - logistic(post$b_dark)
mean(diffBD)
points(diffAD, col = col.alpha("dodgerblue"))
points(diffBD, col = col.alpha("gray"))

## Summary of fixed effects

mean(logistic(post$b_dark))


p <- precis(m2NC)
p2 <- data.frame(p@output)
p3 <- sapply(p2[,1:4], rethinking::logistic)
p3

# Building figure 4 - contrasting raw data with repeated samples from the posterior.
library(rethinking)
library(rstan)
## Load Data ##
#------------------------------------------------------------------------
# d <- readxl::read_excel("../data_untidy/gst_rheo_all.xlsx", na = "NA" )
d <- read.csv("../data_tidy/gst_rheo_all.csv", stringsAsFactors = FALSE, header = TRUE)
head(d)

# Prep Exp 1 Data for visualizing
d2 <- d[d$Type == "with", ] # filter trials from Experiment 1, where barbels were intact. 96 rows.
d2[d2$PropPos == 0,] <- 0.0001
range(d2$PropPos)

# load model data
load("R/m2NC.Rdata")

## Panel Plot with Raw Densities as Fig 4.a, Samples from the Posterior as Fig 4.b

extractvalues <- function(data, treatment) {
  v <- dplyr::filter(data, Treatment == treatment)
  v <- v$PropPos
}

treata <- extractvalues(d2, "A")
treatb <- extractvalues(d2, "B")
treatl <- extractvalues(d2, "L")
treatd <- extractvalues(d2, "D")

# Extract Samples from Posterior
post <- extract.samples(m2NC)

par(mfrow = c(2, 1), mar = c(4, 4, 1.5, 1.5))

dens(treata, xlim = c(0, 1), ylim = c(0, 5), col = "black", adj = 0.7, 
     xlab = "Proportion of Time Spent Positively Oriented", ylab = "Density", lwd = 3)
dens(treatb, col = "purple", lwd = 3, add = TRUE, adj = 0.7)
dens(treatl, lty = 1, lwd = 3, col = "dodgerblue", add = TRUE, adj = 0.7)
dens(treatd, lty = 1, lwd = 3, col = "orange", add = TRUE, adj = 0.7)
# title("Density Curves of Treatment from Empirical Observations")
legend(0, 4.75, cex = 1, pt.cex = 1, c("Above", "Below", "Light", "Dark"), lty = c(1, 1,1,1), lwd = c(3,3,3,3), col = c("black", "purple", "dodgerblue", "orange"), bty = "n")
text(1, 4.6, "A", cex = 2)




## Fig 4.b: sample 300 trials in the posterior for each treatment, fixed effects only
plot(NULL, xlim = c(0, 1) , ylim = c(0, 5) ,
     xlab = "Proportion of Time Spent Positively Oriented", ylab = "Density")
for (i in 1:300)
  curve( dbeta2(x, prob = logistic(post$b_above[i]), theta = post$theta[i]), add = TRUE,
         col = col.alpha("black", 0.15), lwd = 2)
for (i in 1:300)
  curve( dbeta2(x, prob = logistic(post$b_below[i]), theta = post$theta[i]), add = TRUE,
         col = col.alpha("purple", 0.1), lwd = 2)
for (i in 1:300)
  curve( dbeta2(x, prob = logistic(post$b_light[i]), theta = post$theta[i]), add = TRUE,
         col = col.alpha("dodgerblue", 0.1), lwd = 2)
for (i in 1:300)
  curve( dbeta2(x, prob = logistic(post$b_dark[i]), theta = post$theta[i]), add = TRUE,
         col = col.alpha("orange", 0.1), lwd = 2)
text(1, 4.6, "B", cex = 2)

dev.print(tiff, "~/Dropbox/GitHubRepos/gst_rheo/manuscript/figures/Fig4.tiff", res=1000, height=5800, width=7480, units="px")
dev.off()

legend(0.0, 4.75, cex = 1, pt.cex = 1, c("Above", "Below", "Light", "Dark"), lwd = c(1.5, 1.5, 1.5, 1.5), lty = c(1, 1,1,1), col = c("black", "purple", "dodgerblue", "orange"), bty = "n")




## Not Included in Manuscript: Simulating a random fish
plot(NULL, xlim = c(0, 1) , ylim = c(0, 6) ,
     xlab = "Proportion of Time Spent Positively Oriented", ylab = "Density")
for (i in 1:300)
curve( dbeta2(x, prob = logistic(post$b_above[i] + rnorm(1,0,post$sigma_fish[i,1])), theta = post$theta[i]), add = TRUE,  col = col.alpha("black", 0.2))
for (i in 1:300)
  curve( dbeta2(x, prob = logistic(post$b_below[i] + rnorm(1,0,post$sigma_fish[i,1])), theta = post$theta[i]), add = TRUE,  col = col.alpha("purple", 0.1))
for (i in 1:300)
  curve( dbeta2(x, prob = logistic(post$b_light[i] + rnorm(1,0,post$sigma_fish[i,1])), theta = post$theta[i]), add = TRUE,  col = col.alpha("dodgerblue", 0.18))
for (i in 1:300)
  curve( dbeta2(x, prob = logistic(post$b_dark[i] + rnorm(1,0,post$sigma_fish[i,1])), theta = post$theta[i]), add = TRUE,  col = col.alpha("orange", 0.15))
title('Simulating a Random Fish: Repeat Posterior Sampling of m2NC by Treatment', add = TRUE)
legend(0.4, 5.75, cex = 1, pt.cex = 1, c("Above", "Below", "Light", "Dark"), lty = c(1, 1,1,1), col = c("black", "purple", "dodgerblue", "orange"), bty = "n")

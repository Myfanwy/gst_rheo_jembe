## This script parses the barbel and non-barbel trials for the fish that were used in both Experiment 1 and 2, and plots the results.

 library(ggplot2) 
 library(dplyr) # Alternatives to all but one dplyr method used are provided
library(grid)

## Load Data ##
#------------------------------------------------------------------------
# d <- readxl::read_excel("../data_untidy/gst_rheo_all.xlsx", na = "NA" )
d <- read.csv("../data_tidy/gst_rheo_all.csv", stringsAsFactors = FALSE, header = TRUE)
# d2 <- dplyr::filter(d, Type == "with")
head(d)

with <- filter(d, Type == "with")
without <- filter(d, Type == "without")

Fishwith <- as.data.frame(sort(unique(with$FishID))) 
Fishwithout <- as.data.frame(sort(unique(without$FishID)))
names(Fishwith) <- "FishID"
names(Fishwithout) <- "FishID"

reps <- semi_join(Fishwithout, Fishwith) #join for all values of x that are in y (Fishwithout that are in Fishwith)
reps

# visually compare
fw <- sort(unique(with$FishID)) 
fwo <- (sort(unique(without$FishID))) 
fw
fwo
reps$FishID

d <- filter(d, FishID %in% reps$FishID)
means <- d %>% 
  group_by(Type, FishID, Treatment) %>% 
  summarise(means = mean(PropPos))

means <- means %>% 
  group_by(FishID, Type) %>% 
  summarise(meanssum = mean(means))
means
# Plot

p <- ggplot(means) + geom_line(aes(x = Type, y = meanssum, group = FishID, color = as.factor(FishID)), size = 1.0) + 
  theme(legend.position = "none") + ggtitle("Mean proportion positive rheotaxis for velocity treatments \nwith and without barbels") + ylim(c(0, 1)) + ylab("Mean Proportion Positive Rheotaxis \nDuring Velocity Trials (Light and Dark)") + xlab("Barbel Status") 
p
+ 
p + theme(axis.title = element_text(size = 12), axis.text = element_text(size = 14), 
plot.margin = unit(c(1,2,1,1), "cm"), title = element_text(size = 12))
means
library(GGally)
?ggparcoord
diamonds.samp <- diamonds[sample(1:dim(diamonds)[1],100),]
head(diamonds.samp)
# basic parallel coordinate plot, using default settings
 ggparcoord(data = diamonds.samp,columns = c(1,5:10))
ggparcoord(data = means, columns = c(2))
means <- as.data.frame(means)
str(diamonds.samp)


## Graphing the difference between with/without
library(tidyr)
library(ggthemes)
means2 <- means %>% spread(Type, meanssum)
means2$diff <- means2$without - means2$with
head(means2)
str(means2)

p <- ggplot(means2, aes(x = factor(FishID), y = diff)) + geom_bar(aes(alpha = 0.75), stat = "identity")  
p <- p + scale_y_continuous(breaks = c(-0.1, 0, 0.1, 0.2,0.3, 0.4, 0.5), limits=c( -0.1, 0.5))
p                            
p + theme(text = element_text(family = "Times", size = 14), 
          #axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          axis.text = element_text(color = "black", size = 14), 
          axis.line = element_line(colour = "gray"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.y = element_text(vjust = .5), 
          axis.title.x = element_text(vjust = -2),
          plot.margin = (unit(c(1, 1, 1, 1), "cm")),
          panel.background = element_blank()) + 
          ggtitle("Difference in Positive Rheotaxis by Subject After Barbelectomy") + 
          xlab("Subject ID") + ylab("Difference in Proportion Positive Rheotaxis \nduring velocity trials")




## Base R graphics
means2$FishID <- as.factor(means2$FishID)
plot(NULL, xlim = c(0, 11), ylim = c(0, 1))
points(means2$without)
points(means2$with, pch = 16)
means2

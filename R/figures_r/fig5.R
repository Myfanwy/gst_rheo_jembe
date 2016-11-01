library(dplyr)
d <- read.csv("../data_tidy/gst_rheo_all.csv", stringsAsFactors = FALSE, header = TRUE)
head(d)
d2 <- d[d$Type == "with", ] # filter trials from Experiment 1, where barbels were intact. 96 rows.
d2$Velocity <- ifelse(d2$Treatment == "L" | d2$Treatment == "D", 1, 0) # Create dummy variable for velocity
d3 <- d[d$Type == "without", ]

barbs <- d3  %>% 
  group_by(FishID, Treatment)  %>% 
  summarise(mean = mean(PropPos), exp = 2)
barbs
head(d2)

exp1 <- filter(d2, Treatment == "L" | Treatment == "D")

exp1 <- exp1 %>% 
  group_by(FishID, Treatment) %>% 
  summarise(mean = mean(PropPos), exp = 1)
exp1

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

comp <- rbind(barbs, exp1)
comp <- arrange(comp, FishID)
comp <- filter(comp, FishID %in% reps$FishID)

compsum <- comp %>% 
  group_by(FishID, exp) %>% 
  summarise(meanr = mean(mean))

ggplot(compsum, aes(x = as.factor(exp), y = meanr, group = FishID)) + geom_line(aes(color = as.factor(FishID)))

range(d2$FishID)

p <- compsum %>% 
  group_by(FishID) %>% 
  summarise(change = meanr[2] - meanr[1]) %>% 
  ggplot(aes(x = as.factor(FishID), y = change)) + geom_bar(stat = "identity")
p

p <- p + labs(x = "Subject ID", y = "Difference in Proportion Positive Rheotaxis \nBetween Experiment 1 and Experiment 2", 
         title = "Difference in Positive Rheotaxis by Subject After Barbelectomy") +
    theme_bw() +
  ylim(c(-0.25, 0.75)) + 
  theme(text = element_text(family = "Arial", size = 82),
              axis.line = element_line(colour = "black", size = 3),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.title = element_text(vjust= 0.5, margin = margin(15, 15, 15, 15)), 
              plot.margin = (unit(c(0.5, 0.5, 0.5, 0.5), "cm")) ,
              panel.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
        text = element_text(family = "Arial"))
 p 
 
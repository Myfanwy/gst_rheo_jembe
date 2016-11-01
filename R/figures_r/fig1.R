# Fig 1

library(devtools)
library(ggplot2)
library(grid)

# d <- readxl::read_excel("../data_untidy/gst_rheo_all.xlsx", na = "NA" )
d <- read.csv("../data_tidy/gst_rheo_all.csv", stringsAsFactors = FALSE, header = TRUE)
d <- d[d$Type == "with", ]

## After Convo with Pete ##


p <-  ggplot(d, aes(x = factor(Treatment), y = PropPos)) + geom_boxplot(aes(group = Treatment, size = 12), outlier.colour = "transparent", position = position_dodge(width = 0), width = 0.45, alpha = 0.75)

p <- p + geom_point(aes(alpha = 0.8), size = 15, position = position_jitter(width = 0.10))  + scale_x_discrete(labels = c("Above", "Below", "Dark", "Light")) + 
  labs(y = "Proportion Positive Rheotaxis During Trial\n") + xlab("\nTreatment") +  ggtitle("Proportion of Time Spent Positively Oriented by Treatment: \nExperiment 1") +
  theme(text = element_text(family = "Arial", size = 82),
    axis.line = element_line(colour = "black", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(vjust= 0.5, margin = margin(15, 15, 15, 15)), 
        plot.margin = (unit(c(0.5, 0.5, 0.5, 0.5), "cm")) ,
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5))  + 
  theme(legend.position = "none") + geom_hline(yintercept = 0.33, lty = 3, lwd = 3)
p

ggsave(p, filename = "~/Dropbox/GitHubRepos/gst_rheo/manuscript/figures/Fig3.tiff", width = 5512,
       height = 3997, dpi = 72, limitsize = FALSE)

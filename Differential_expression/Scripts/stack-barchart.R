
library(ggplot2)
library(reshape)
library(patchwork)

data <- read.csv("~/Postdoc/TCE-paper/NLR-number-allcultivars-alltimepoints.csv", header=TRUE)

mdata<-melt(data)
write.csv(mdata, "nlr-degs-melted.csv")
data <- read.csv("nlr-degs-melted.csv", header=TRUE)

data<-data[data$dpi == "1dpi",]
data$variable <- factor(data$variable ,levels = c("OA", "SO", "SA"))
dpi1 <- ggplot(data, aes(x = variable, y =value, fill=isolate)) + geom_bar(stat = "identity",  width=.8) +
  xlab("") + ylab("Number of DEGs") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle("1 dpi")

dpi1<-dpi1+ scale_y_continuous(limits = c(0, 3000)) +theme(legend.position = "none", plot.title = element_text(size = 14), axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
dpi1



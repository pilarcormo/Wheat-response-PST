library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)

samples <- read.csv("~/Desktop/samples-2013-2015.CSV")

samples <- read.csv("~/Desktop/RNA-seq_samples_1314.csv")

summary(samples)

my.col <- colorRampPalette(brewer.pal(11, "RdBu"))

mypalette<-brewer.pal(10,"Spectral")

########COUNTRIES

df <- data.frame(table(samples$Location))
newdata <- df[order(df$Freq),] 
df2 <- data.frame(newdata$Var1, newdata$Freq)

x <- df2$newdata.Var1
y <- df2$newdata.Freq*100/sum(df2$newdata.Freq)
df3 <- data.frame(x,y)

p_percent <- ggplot(df3, aes(x = reorder(x, y), y = y)) + 
  xlab("") + ylab("%") +theme_bw()  + geom_bar(stat = "identity",  width=.8, fill="darkgreen") + ggtitle(" ") +
  coord_flip()  + theme(plot.title = element_text(size = 20), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"))

p_percent


p_absolute <- ggplot(df2, aes(x = reorder(newdata.Var1, newdata.Freq), y = newdata.Freq)) + 
  xlab("") + ylab("Counts") +theme_bw() + geom_bar(stat = "identity",  width=.8, fill="darkgreen") + ggtitle(" ") + 
  coord_flip()  + theme(axis.ticks = element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"), axis.title=element_text(size=16,face="bold")) 

p_absolute


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 159), guide = "colourbar") 

###colour bar
p_absolute <- ggplot(df2, aes(x = reorder(newdata.Var1, newdata.Freq), y = newdata.Freq, fill = as.integer(newdata.Freq))) + 
  xlab("") + ylab("Counts") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle(" ") + 
  coord_flip()  + theme(axis.ticks = element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"), axis.title=element_text(size=16,face="bold")) + 
  sc + theme(legend.text = element_text(size = 18)) + guides(colour = guide_colorbar(barwidth = 6, barheight = 20))
######


########VARIETIES

#Choose varieties isolated more than once (1) or twice (2)
df_variety <- data.frame(table(samples$Variety))
order <- df_variety[order(df_variety$Freq),] 
x <- order[order$Freq > 2,]$Var1
y <- order[order$Freq > 2,]$Freq

y_3 <- y*100/sum(y)

df2_variety <- data.frame(x, y)

df3_variety <- data.frame(x, y_3)

v_absolute <- ggplot(df2_variety, aes(x = reorder(x, y), y = y)) + 
  xlab("") + ylab("Counts") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle("Varieties 2013-2015") + 
  coord_flip() + theme(plot.title = element_text(size = 20), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold")) 
 

v_percent <- ggplot(df3_variety, aes(x = reorder(x, y_3), y = y_3)) + 
  xlab("") + ylab("Frequency") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle("Varieties 2013-2015") + 
  coord_flip() + theme(plot.title = element_text(size = 20), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold")) 


########VARIETIES for wheat
wheat <- data.frame(table(samples[samples$Species== 'Wheat',]$Variety))
triticale <- data.frame(table(samples[samples$Species== 'Triticale',]$Variety))
#Choose varieties isolated more than once (1) or twice (2)
order <- triticale[order(triticale$Freq),] ###change for wheat, tritricale or sample
x <- order[order$Freq > 0,]$Var1
y <- order[order$Freq > 0,]$Freq*100/sum(order$Freq)
y
y_3 <- order[order$Freq > 2,]$Freq
df2_variety <- data.frame(x, y)
df3_variety <- data.frame(x, y_3)

v_percent <- ggplot(df2_variety, aes(x = reorder(x, y), y = y)) + 
  xlab("") + ylab("%") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle("Varieties 2013-2015") + 
  coord_flip() + theme(plot.title = element_text(size = 20), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"))

v_absolute <- ggplot(df3_variety, aes(x = reorder(x, y), y = y_3)) + 
  xlab("") + ylab("Frequency") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle("Varieties 2013-2015") + 
  coord_flip() + theme(axis.text=element_text(size=12), axis.title=element_text(size=16,face="bold"))






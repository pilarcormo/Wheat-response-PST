library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)

samples <- read.csv("~/Desktop/samples-2013-2015.CSV")

summary(samples)



samples <- read.csv("~/last_samples20131415.csv")

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
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 50), guide = "colourbar") 

###colour bar
p_absolute <- ggplot(df2, aes(x = reorder(newdata.Var1, newdata.Freq), y = newdata.Freq, fill = as.integer(newdata.Freq))) + 
  xlab("") + ylab("Counts") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle(" ") + 
  coord_flip()  + theme(axis.ticks = element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"), axis.title=element_text(size=16,face="bold")) + 
   sc + theme(legend.text = element_text(size = 18)) + guides(colour = guide_colorbar(barwidth = 6, barheight = 20))
######

  coord_flip()  + theme(plot.title = element_text(size = 20), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold")) + 
  geom_text(df2, aes(x = reorder(newdata.Var1, newdata.Freq), y = newdata.Freq), size = 3) +  scale_colour_gradient()

colors = brewer.pal(9, "BrBG")
pal <- colorRampPalette(colors)


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 50), guide = "colourbar") 
sc

p_absolute <- ggplot(df2, aes(x = reorder(newdata.Var1, newdata.Freq), y = newdata.Freq, fill = as.integer(newdata.Freq))) + 
  xlab("") + ylab("Counts") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle(" ") + 
  coord_flip()  + theme(axis.ticks = element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"), axis.title=element_text(size=16,face="bold")) + 
  sc + theme(legend.text = element_text(size = 18)) + guides(colour = guide_colorbar(barwidth = 6, barheight = 20))


p_absolute <- ggplot(df2, aes(x = reorder(newdata.Var1, newdata.Freq), y = newdata.Freq, fill = as.integer(newdata.Freq))) + 
  xlab("") + ylab("Counts") +theme_bw() + geom_tile() + ggtitle(" ") + 
  coord_flip()  + theme(axis.ticks = element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"), axis.title=element_text(size=16,face="bold")) + 
  scale_x_discrete(breaks=NULL) 

p_absolute

grid <- grid.arrange(p_percent , p_absolute, ncol=2, heights=c(1, 10), widths =c(2,1), as.table =TRUE)




########VARIETIES

#Choose varieties isolated more than once (1) or twice (2)
df_variety <- data.frame(table(samples$Variety))
order <- df_variety[order(df_variety$Freq),] 
x <- order[order$Freq > 2,]$Var1
y <- order[order$Freq > 2,]$Freq

y_3 <- y*100/sum(y)

df2_variety <- data.frame(x, y)

df3_variety <- data.frame(x, y_3)

df2_variety

v_absolute <- ggplot(df2_variety, aes(x = reorder(x, y), y = y)) + 
  xlab("") + ylab("Counts") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle("Varieties 2013-2015") + 
  coord_flip() + theme(plot.title = element_text(size = 20), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold")) 
v_absolute 

v_percent <- ggplot(df3_variety, aes(x = reorder(x, y_3), y = y_3)) + 
  xlab("") + ylab("Frequency") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle("Varieties 2013-2015") + 
  coord_flip() + theme(plot.title = element_text(size = 20), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold")) 


df2_variety <- data.frame(x, y)


v_percent <- ggplot(df2_variety, aes(x = reorder(x, y), y = y)) + 
  xlab("") + ylab("Counts") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle("Wheat varieties 2013-2015") + 
  coord_flip() + theme(plot.title = element_text(size = 20), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold")) 
 

v_absolute <- ggplot(df3_variety, aes(x = reorder(x, y), y = y_3)) + 
  xlab("") + ylab("Frequency") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle("Varieties 2013-2015") + 
  coord_flip() + theme(axis.text=element_text(size=12), axis.title=element_text(size=16,face="bold"))


########VARIETIES for wheat

wheat <- data.frame(table(samples[samples$Species== 'Wheat',]$Variety))
triticale <- data.frame(table(samples[samples$Species== 'Triticale',]$Variety))
#Choose varieties isolated more than once (1) or twice (2)
order <- wheat[order(wheat$Freq),] ###change for wheat, tritricale or sample
x <- order[order$Freq > 2,]$Var1
y <- order[order$Freq > 0,]$Freq*100/sum(order$Freq)

y_3 <- order[order$Freq > 2,]$Freq
df2_variety <- data.frame(x, y)
df3_variety <- data.frame(x, y_3)

v_percent <- ggplot(df2_variety, aes(x = reorder(x, y), y = y)) + 
  xlab("") + ylab("%") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle("Varieties 2013-2015") + 
  coord_flip() + theme(plot.title = element_text(size = 20), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"))



######STACKED PLOT BY YEAR 
only_wheat <- samples[samples$Species== 'Wheat',]
wheat_2013 <- data.frame(table(only_wheat[only_wheat$Date=='2013',]$Variety))
wheat_2014 <- data.frame(table(only_wheat[only_wheat$Date=='2014',]$Variety))
wheat_2015 <- data.frame(table(only_wheat[only_wheat$Date=='2015',]$Variety))

l <- data.frame(wheat_2013$Var1, wheat_2013$Freq, wheat_2014$Freq, wheat_2015$Freq, wheat$Freq)
write.csv(l, "var.csv", quote=FALSE)

new <- read.csv("var.csv")

y_total <- new[new$Total > 2,]$Total
x_var <- new[new$Total > 2,]$Var
x_13 <- new[new$Total > 2,]$X2013
x_14<- new[new$Total > 2,]$X2014
x_15 <- new[new$Total > 2,]$X2015

df <- data.frame(x_var, x_13, x_14, x_15)

library(reshape2)
dat.m <- melt(df, id.vars = "x_var") ## just melt(dat) should work
ggplot(dat.m, aes(x = x_var, y =value,fill=variable)) +
  geom_bar(stat='identity')


v_absolute <- ggplot(dat.m, aes(x = reorder(x_var, value), y = value, fill=variable)) +  xlab("") + ylab("Counts") +
  theme_bw() + geom_bar(stat = "identity",  width=.8)  + 
  coord_flip() + theme(axis.text=element_text(size=12), axis.title=element_text(size=16,face="bold"))


colors = brewer.pal(3, "YlGn")
v_absolute + scale_fill_manual(values=colors, 
                                 name="Year",
                                breaks=c("x_13", "x_14", "x_15"),
                               labels=c("2013", "2014", "2015"))


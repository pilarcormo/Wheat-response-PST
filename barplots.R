library(ggplot2)

samples <- read.csv("~/Desktop/samples-2013-2015.csv")

samples <- read.csv("~/Desktop/RNA-seq_samples_1314.csv")

summary(samples)

########COUNTRIES

df <- data.frame(table(samples$Location))
newdata <- df[order(df$Freq),] 
df2 <- data.frame(newdata$Var1, newdata$Freq)

x <- df2$newdata.Var1
y <- df2$newdata.Freq*100/sum(df2$newdata.Freq)
y
df3 <- data.frame(x,y)

#barplot(df2$newdata.Freq, main="Country",ylab="Frequency",las=2)
#labs <- paste(names(df2$newdata.Var1))
#text(cex=1, x=x-.25, y=-1.25, labs, xpd=TRUE, srt=45, pos=2)


p_percent <- ggplot(df3, aes(x = reorder(x, y), y = y)) + 
  xlab("") + ylab("%") +theme_bw()  + geom_bar(stat = "identity",  width=.8) + ggtitle("Countries 2013-2015") +
  coord_flip()  + theme(plot.title = element_text(size = 20), axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"))

p_absolute <- ggplot(df2, aes(x = reorder(newdata.Var1, newdata.Freq), y = newdata.Freq)) + 
  xlab("") + ylab("Frequency") +theme_bw() + geom_bar(stat = "identity",  width=.8) + ggtitle("Countries 2013-2014") + 
  coord_flip()  + theme(axis.text=element_text(size=12), axis.title=element_text(size=16,face="bold"))


########VARIETIES

#Choose varieties isolated more than once (1) or twice (2)
df_variety <- data.frame(table(samples$Variety))
order <- df_variety[order(df_variety$Freq),] 
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



LIB <- read.csv("~/Identify_varieties/final_scores/final_scores_LIB4362.csv")
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colSort <- function(data, ...) sapply(data, sort, ...)
varieties <- colSums(LIB[,-1])
df <- data.frame(varieties)
colMax(df)
csv <- write.csv(df, "try.csv", quote=FALSE)
caca <- read.csv("~/try.csv")
df2 <- data.frame(caca$X, caca$varieties)
df2



caracol <- ggplot(df2, aes(x = reorder(caca.X, caca.varieties), y = caca.varieties)) +
  geom_bar(width = 1, stat = "identity")  +
  coord_polar() +theme_bw() + ylab(" ") + xlab(" ") 


caracol + scale_fill_discrete("darkseagreen1") + theme(plot.title = element_text(size = 16), axis.text=element_text(size=18), axis.title=element_text(size=16,face="bold"))

colors = brewer.pal(3, "YlGn")

caracol + scale_fill_manual(values=colors)

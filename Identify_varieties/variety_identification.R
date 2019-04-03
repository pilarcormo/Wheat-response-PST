
LIB <- read.csv("final_scores_test.csv")
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colSort <- function(data, ...) sapply(data, sort, ...)
varieties <- colSums(LIB[,-1])
df <- data.frame(varieties)
write.csv(df, "test.csv")
colMax(df)
varieties


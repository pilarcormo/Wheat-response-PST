LIB <- read.csv("final_scores/final_scores.csv")
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colSort <- function(data, ...) sapply(data, sort, ...)

varieties <- colSums(LIB[,-1])

df <- data.frame(varieties)
colMax(df)

varieties
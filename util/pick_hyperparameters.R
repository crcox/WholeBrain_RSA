library('dplyr')
library('ggplot2')
library('jsonlite')
d <- read.csv('./tune.csv', header=TRUE)
d$subject <- as.factor(d$subject)
d <- filter(d, cvholdout != 10, finalholdout != 10)

tmp <- na.omit(select(d, subject,finalholdout,cvholdout,lambda,lambda1,err1,nzv,nvox)) %>%
  group_by(subject,finalholdout,lambda,lambda1) %>%
  summarize(err = mean(err1), pnzv = mean(nzv/nvox)) %>%
  group_by(subject,finalholdout) %>%
  summarize(ind = which.min(err)[1], lambda = lambda[ind], lambda1 = lambda1[ind], err = err[ind], pnzv = pnzv[ind]) %>%
  ungroup()

write.csv(x = tmp, file = "best-parameters.csv")
jsonlite::write_json(x = tmp, path = "best-parameters.json", dataframe = 'columns', pretty = TRUE)

subjects <- levels(d$subject)
holdouts <- unique(d$finalholdout)
P <- rep(list(), length(subjects))
 <- rep(NA, length(subjects))
for (i in 1:length(subjects)) {
  s <- subjects[i]
  lambda <- 
}

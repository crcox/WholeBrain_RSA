library('dplyr')
library('ggplot2')
library('jsonlite')
tab <- function(indent, tabwidth = 4) {
    return(indent * tabwidth)
}
d <- read.csv('./tune.csv', header=TRUE)
d$subject <- as.factor(d$subject)
d <- filter(d, cvholdout != 10, finalholdout != 10)

tmp <- na.omit(select(d, subject,finalholdout,cvholdout,lambda,lambda1,err1,nzv,nvox)) %>%
  group_by(subject,finalholdout,lambda,lambda1) %>%
  summarize(err = mean(err1), pnzv = mean(nzv/nvox)) %>%
  group_by(subject,finalholdout) %>%
  summarize(ind = which.min(err)[1], lambda = lambda[ind], lambda1 = lambda1[ind], err = err[ind], pnzv = pnzv[ind]) %>%
  ungroup()

write.csv(x = tmp, file = "best-parameters-bysubject.csv")
jsonlite::write_json(x = tmp, path = "best-parameters-bysubject.json", dataframe = 'columns', pretty = TRUE)

subjects <- levels(d$subject)
holdouts <- unique(d$finalholdout)
P <- rep(list(list(lambda = rep(0, length(holdouts)),lambda1 = rep(0, length(holdouts)))), length(subjects))
for (i in 1:length(subjects)) {
  s <- subjects[i]
  P[[i]]$lambda <- filter(tmp, subject == s)$lambda
  P[[i]]$lambda1 <- filter(tmp, subject == s)$lambda1
}

fid <- file(description = './best-parameters-bysubject.yaml', open = 'w')
cat('lambda:\n', file = fid)
for (i in 1:length(subjects)) cat(paste('    - [',paste(P[[i]]$lambda, collapse=','),']\n', sep=''), file = fid)

cat('lambda1:\n', file = fid)
for (i in 1:length(subjects)) cat(paste('    - [',paste(P[[i]]$lambda1, collapse=','),']\n', sep=''), file = fid)
close(fid)

lambda <- unique(d$lambda)
lambda1 <- unique(d$lambda1)
fid <- file(description = './README.best-parameters-bysubject', open = 'w')
msg <- c(
    strwrap("These parameters were selected by searching the space defined by crossing:", width = 80, indent = tab(0)),
    strwrap(sprintf("lambda: [%s]\n", paste(lambda, collapse = ',')), width = 80, indent = tab(1), exdent = tab(2)),
    strwrap(sprintf("lambda1: [%s]\n", paste(lambda1, collapse = ',')), width = 80, indent = tab(1), exdent = tab(2)),
    "",
    strwrap(paste(
        sprintf("This space was searched for each of %d subjects and %d holdout sets.", length(subjects), length(holdouts)),
        sprintf("The decision was made based on a performance measure obtained by %d-fold cross validation", length(holdouts) - 1),
        "(within each subject, holdoutset, and combination of parameters)."), width = 80, indent = tab(0)),
    "",
    strwrap(paste(
        "The selected parameters are represented as csv, json, and (perhaps most immediately useful) yaml.",
        "The yaml formatted parameter values have one row per subject, and the elements of each row correspond to a holdout set."), width = 80, indent = tab(0)),
    "",
    strwrap("Copy the contents of 'best-parameters-bysubject.yaml' into your stub.yaml file (replacing the relevant parameters).", width = 80, indent = tab(0)),
    "",
    strwrap(paste(
        "Under the EXPAND field (in the stub.yaml), you want to specify that each row of each hyper-parameter field",
        "is associated with a different dataset, and that each column is associated with a different holdout set.",
        "This can be achieved with something like the following syntax:"), width = 80, indent = tab(0)),
    strwrap("EXPAND:", width = 80, indent = tab(1)),
    strwrap("- [data,lambda,lambda1]", width = 80, indent = tab(2)),
    strwrap("- [cvholdout,lambda_,lambda1_]", width = 80, indent = tab(2)),
    "",
    strwrap(paste(
        "Sublists under EXPAND specify that the elements of the referenced fields should be iterated over with the",
        "same index. So, if there are 23 subjects and 9 cvholdouts, then this specification amounts to:"), width = 80, indent = tab(0)),
    strwrap("for i = 1:23", width = 80, indent = tab(1)),
    strwrap("for j = 1:9", width = 80, indent = tab(2)),
    strwrap("d  = data(i);", width = 80, indent = tab(3)),
    strwrap("c  = cvholdout(j);", width = 80, indent = tab(3)),
    strwrap("L  = lambda(i,j);", width = 80, indent = tab(3)),
    strwrap("L1 = lambda1(i,j);", width = 80, indent = tab(3)),
    strwrap("end", width = 80, indent = tab(2)),
    strwrap("end", width = 80, indent = tab(1))
)
writeLines(msg, con = fid)
close(fid)

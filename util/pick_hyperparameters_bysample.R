library('dplyr')
library('ggplot2')
library('jsonlite')
tab <- function(indent, tabwidth = 4) {
    return(indent * tabwidth)
}
d <- read.csv('./tune.csv', header=TRUE)
d$subject <- as.factor(d$subject)
d <- filter(d, cvholdout != 10, finalholdout != 10)

tmp <- na.omit(select(d, finalholdout,cvholdout,lambda,lambda1,err1,nzv,nvox)) %>%
  group_by(finalholdout,lambda,lambda1) %>%
  summarize(err = mean(err1), pnzv = mean(nzv/nvox)) %>%
  group_by(finalholdout) %>%
  summarize(ind = which.min(err)[1], lambda = lambda[ind], lambda1 = lambda1[ind], err = err[ind], pnzv = pnzv[ind]) %>%
  ungroup()

write.csv(x = tmp, file = "best-parameters-bysample.csv")
jsonlite::write_json(x = tmp, path = "best-parameters-bysample.json", dataframe = 'columns', pretty = TRUE)

fid <- file(description = './best-parameters-bysample.yaml', open = 'w')
msg <- c(
    paste('lambda: [',paste(tmp$lambda, collapse=','),']', sep=''),
    paste('lambda1: [',paste(tmp$lambda1, collapse=','),']', sep='')
)
writeLines(msg, con = fid)
close(fid)

subjects <- unique(d$subject)
holdouts <- unique(d$finalholdout)
lambda <- unique(d$lambda)
lambda1 <- unique(d$lambda1)
fid <- file(description = './README.best-parameters-bysample', open = 'w')
msg <- c(
    strwrap("These parameters were selected by searching the space defined by crossing:", width = 80, indent = tab(0)),
    strwrap(sprintf("lambda: [%s]\n", paste(lambda, collapse = ',')), width = 80, indent = tab(1), exdent = tab(2)),
    strwrap(sprintf("lambda1: [%s]\n", paste(lambda1, collapse = ',')), width = 80, indent = tab(1), exdent = tab(2)),
    "",
    strwrap(paste(
        sprintf("This space was searched for each of %d holdout sets, aggregating over subjects.", length(holdouts)),
        sprintf("The decision was made based on a performance measure obtained by %d-fold cross validation", length(holdouts) - 1),
        "(within each subject, holdoutset, and combination of parameters)."), width = 80, indent = tab(0)),
    "",
    strwrap(paste(
        "The selected parameters are represented as csv, json, and (perhaps most immediately useful) yaml.",
        "The yaml formatted parameter values correspond to the best performing values for each holdout set."), width = 80, indent = tab(0)),
    "",
    strwrap("Copy the contents of 'best-parameters-bysample.yaml' into your stub.yaml file (replacing the relevant parameters).", width = 80, indent = tab(0)),
    "",
    strwrap(paste(
        "Under the EXPAND field (in the stub.yaml), you want to specify that each hyper-parameter value",
        "is associated with a different holdout set, and that all subjects will be fit with the same hyperparameters.",
        "This can be achieved with something like the following syntax:"), width = 80, indent = tab(0)),
    strwrap("EXPAND:", width = 80, indent = tab(1)),
    strwrap("- [cvholdout,lambda,lambda1]", width = 80, indent = tab(2)),
    strwrap("- data", width = 80, indent = tab(2)),
    "",
    strwrap(paste(
        "Sublists under EXPAND specify that the elements of the referenced fields should be iterated over with the",
        "same index. So, if there are 23 subjects and 9 cvholdouts, then this specification amounts to:"), width = 80, indent = tab(0)),
    strwrap("for i = 1:23", width = 80, indent = tab(1)),
    strwrap("for j = 1:9", width = 80, indent = tab(2)),
    strwrap("d  = data(i);", width = 80, indent = tab(3)),
    strwrap("c  = cvholdout(j);", width = 80, indent = tab(3)),
    strwrap("L  = lambda(j);", width = 80, indent = tab(3)),
    strwrap("L1 = lambda1(j);", width = 80, indent = tab(3)),
    strwrap("end", width = 80, indent = tab(2)),
    strwrap("end", width = 80, indent = tab(1))
)
writeLines(msg, con = fid)
close(fid)

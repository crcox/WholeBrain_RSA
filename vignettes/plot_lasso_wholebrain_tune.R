library(dplyr)
library(ggplot2)

datadir <- file.path('results','lasso','tune')
stopifnot(dir.exists(datadir))

datafile <- file.path(datadir, 'wholebrain_tune.csv')
stopifnot(file.exists(datafile))

d <- read.csv(datafile, header=T)

# Average over cvholdouts
# NB. lambda values may differ by subject. The GLMNET default is to attempt 100 lambda at a time.
nlambda <- 100
nsubject <- 10
nfinalholdout <- 10
ncvholdout <- nfinalholdout - 1
ntarget <- 3
d$lamind <- rep(1:nlambda,nsubject*nfinalholdout*ncvholdout*ntarget)

d <- d %>%
  group_by(subject,finalholdout,lamind,target) %>%
  summarize(nt1=mean(nt1), nt2=mean(nt2), nd1=mean(nd1), nd2=mean(nd2),
            h1=mean(h1), h2=mean(h2), f1=mean(f1), f2=mean(f2),
            Wnz=mean(Wnz), nvox=mean(nvox))

# NB. Lambda are in *descending* order. So the sparsest solution is on the *left*
# Plot sensitivity by lambda, with lines for final holdout sets and panel for subjects
ggplot(filter(d, target=="TrueFaces"), aes(x=lamind, y=(h1/nt1)-(f1/nd1), color=as.factor(finalholdout))) + geom_line() + facet_wrap("subject") + theme_bw(base_size=18) + ggtitle("Lasso - Face - Whole Brain") + scale_fill_discrete(name="Final\nHoldout")
ggplot(filter(d, target=="TruePlaces"), aes(x=lamind, y=(h1/nt1)-(f1/nd1), color=as.factor(finalholdout))) + geom_line() + facet_wrap("subject") + theme_bw(base_size=18) + ggtitle("Lasso - Place - Whole Brain")
ggplot(filter(d, target=="TrueObjects"), aes(x=lamind, y=(h1/nt1)-(f1/nd1), color=as.factor(finalholdout))) + geom_line() + facet_wrap("subject") + theme_bw(base_size=18) + ggtitle("Lasso - Object - Whole Brain")

# Plot number of voxels selected by lambda, with lines for final holdout sets and panel for subjects
ggplot(filter(d, target=="TrueFaces"), aes(x=lamind, y=Wnz, color=as.factor(finalholdout))) + geom_line() + facet_wrap("subject") + theme_bw(base_size=18) + ggtitle("Lasso - Face - Whole Brain")
ggplot(filter(d, target=="TruePlaces"), aes(x=lamind, y=Wnz, color=as.factor(finalholdout))) + geom_line() + facet_wrap("subject") + theme_bw(base_size=18) + ggtitle("Lasso - Place - Whole Brain")
ggplot(filter(d, target=="TrueObjects"), aes(x=lamind, y=Wnz, color=as.factor(finalholdout))) + geom_line() + facet_wrap("subject") + theme_bw(base_size=18) + ggtitle("Lasso - Object - Whole Brain")

# Select best lambda
d <- read.csv(datafile, header=T)
d <- d %>%
  group_by(subject,finalholdout,cvholdout,target) %>%
  summarize(
    lamind=which.max((h1/nt1)-(f1/nd1)),lambda=lambda[lamind],
    nt1=nt1[lamind], nt2=nt2[lamind], nd1=nd1[lamind], nd2=nd2[lamind],
    h1=h1[lamind], h2=h2[lamind], f1=f1[lamind], f2=f2[lamind],
    Wnz=Wnz[lamind], nvox=nvox[lamind])

ggplot(filter(d,target=="TrueFaces"), aes(x=as.factor(finalholdout), y=(h1/nt1)-(f1/nd1))) + geom_point() + facet_wrap("subject")
ggplot(filter(d,target=="TrueFaces"), aes(x=as.factor(finalholdout), y=lambda)) + geom_point() + facet_wrap("subject")

best_lambda <- d %>%
  group_by(subject,target) %>%
  summarize(lambda=mean(lambda))

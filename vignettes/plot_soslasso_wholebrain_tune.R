library(dplyr)
library(ggplot2)

datadir <- file.path('results','soslasso','tune')
stopifnot(dir.exists(datadir))

datafile <- file.path(datadir, 'soslasso_face_wholebrain.csv')
stopifnot(file.exists(datafile))
face <- read.csv(datafile, header=T)
face$target <- factor(1, levels=1:3, labels=c('face','place','object'))

datafile <- file.path(datadir, 'soslasso_place_wholebrain.csv')
stopifnot(file.exists(datafile))
place <- read.csv(datafile, header=T)
place$target <- factor(2, levels=1:3, labels=c('face','place','object'))

datafile <- file.path(datadir, 'soslasso_object_wholebrain.csv')
stopifnot(file.exists(datafile))
object <- read.csv(datafile, header=T)
object$target <- factor(3, levels=1:3, labels=c('face','place','object'))

d <- rbind(face,place,object)
rm(face,place,object)

# Average over cvholdouts
# NB. lambda values may differ by subject. The GLMNET default is to attempt 100 lambda at a time.
nlambda <- n_distinct(d$lambda)
nalpha <- n_distinct(d$alpha)
nsubject <- n_distinct(d$subject)
nfinalholdout <- n_distinct(d$finalholdout)
ncvholdout <- nfinalholdout - 1
ntarget <- n_distinct(d$target)

d <- d %>%
  group_by(subject,lambda,alpha,target) %>%
  summarize(nt1=mean(nt1), nt2=mean(nt2), nd1=mean(nd1), nd2=mean(nd2),
            h1=mean(h1), h2=mean(h2), f1=mean(f1), f2=mean(f2),
            Wnz=mean(Wnz))

# NB. Lambda are in *descending* order. So the sparsest solution is on the *left*
# Plot sensitivity by lambda, with lines for final holdout sets and panel for subjects
ggplot(filter(d, target=="face"), aes(x=lambda, y=(h1/nt1)-(f1/nd1), color=as.factor(subject))) + geom_line() + facet_wrap("alpha") + theme_bw(base_size=18) + ggtitle("Lasso - Face - Whole Brain") + scale_fill_discrete(name="Subject")
ggplot(filter(d, target=="TruePlaces"), aes(x=lamind, y=(h1/nt1)-(f1/nd1), color=as.factor(finalholdout))) + geom_line() + facet_wrap("subject") + theme_bw(base_size=18) + ggtitle("Lasso - Place - Whole Brain")
ggplot(filter(d, target=="TrueObjects"), aes(x=lamind, y=(h1/nt1)-(f1/nd1), color=as.factor(finalholdout))) + geom_line() + facet_wrap("subject") + theme_bw(base_size=18) + ggtitle("Lasso - Object - Whole Brain")

# Plot number of voxels selected by lambda, with lines for final holdout sets and panel for subjects
ggplot(filter(d, target=="TrueFaces"), aes(x=lamind, y=Wnz, color=as.factor(finalholdout))) + geom_line() + facet_wrap("subject") + theme_bw(base_size=18) + ggtitle("Lasso - Face - Whole Brain")
ggplot(filter(d, target=="TruePlaces"), aes(x=lamind, y=Wnz, color=as.factor(finalholdout))) + geom_line() + facet_wrap("subject") + theme_bw(base_size=18) + ggtitle("Lasso - Place - Whole Brain")
ggplot(filter(d, target=="TrueObjects"), aes(x=lamind, y=Wnz, color=as.factor(finalholdout))) + geom_line() + facet_wrap("subject") + theme_bw(base_size=18) + ggtitle("Lasso - Object - Whole Brain")

# Select best lambda
#d <- read.csv(datafile, header=T)
load(file.path(datadir,"soslasso_wholebrain.Rdat"))
d <- d %>%
  group_by(subject,finalholdout,cvholdout,target) %>%
  summarize(
    ind=which.max((h1/nt1)-(f1/nd1)),lambda=lambda[ind],alpha=alpha[ind],
    nt1=nt1[ind], nt2=nt2[ind], nd1=nd1[ind], nd2=nd2[ind],
    h1=h1[ind], h2=h2[ind], f1=f1[ind], f2=f2[ind],
    Wnz=Wnz[ind])

ggplot(filter(d,target=="TrueFaces"), aes(x=as.factor(finalholdout), y=(h1/nt1)-(f1/nd1))) + geom_point() + facet_wrap("subject")
ggplot(filter(d,target=="TrueFaces"), aes(x=as.factor(finalholdout), y=lambda)) + geom_point() + facet_wrap("subject")

best_lambda <- d %>%
  group_by(target) %>%
  summarize(lambda=mean(lambda),alpha=mean(alpha))

library("dplyr")
library("ggplot2")

datadir <- file.path("Manchester","results","WholeBrain_RSA","semantic","audio","grOWL2","tune")
stopifnot(file.exists(datadir))

datapath <- file.path(datadir, "GrOWL2_tune_audio_semantic.csv")
stopifnot(file.exists(datapath))

d <- read.csv(datapath)
d <- d[!is.na(d$lambda),]
str(d)

d.avg <- d %>%
  # filter(finalholdout==3) %>%
  group_by(subject,lambda,lambda1) %>%
  summarize(err=mean(err1,na.rm=TRUE)) %>% ungroup()

p <- ggplot(d.avg, aes(x=lambda,y=err,color=as.factor(lambda1))) + geom_line() + facet_wrap("subject")
png("audio_errorByLambda.png", width=800, height=700)
print(p)
dev.off()

ggplot(d.avg, aes(x=as.factor(lambda),y=as.factor(lambda1))) + geom_tile(aes(fill=err)) + facet_wrap("subject")
d.tune <- d %>%
  group_by(subject,finalholdout,lambda1,lambda) %>%
  summarize(err=mean(err1, na.rm=TRUE)) %>%
  group_by(subject,finalholdout) %>%
  summarize(lambda1=lambda1[which.min(err)],lambda=lambda[which.min(err)])
summary(d)

paste(d.tune$subject,collapse = ",")
paste(d.tune$finalholdout,collapse = ",")
paste(d.tune$lambda1,collapse = ",")
paste(d.tune$lambda,collapse = ",")

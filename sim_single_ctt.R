library(irtimv)
source("00funs.R")

##########################################################
##one example calculation
b0<-0
b1<-1
b2<-2
a<-1
bs<-c(b0, b1,b2)
set.seed(123)
th<-rnorm(10000)

x <- generate_grm_data(th, a, b=bs)
x2 <- generate_grm_data(th, a, b=bs)

grmpars<-estimate_grm(x,th)
pcmpars<-estimate_pcm(x,th)
##  plot(th,grmpr[,1])
##  for (i in 1:4) points(th,grmpr[,i])
##  for (i in 1:4) points(th,pcmpr[,i],col='red')

grmpr<-predict_grm(th,grmpars$a,grmpars$b)
pcmpr<-predict_pcm(th,pcmpars$a,pcmpars$d)

## Calculate frequency table
pctt.tab <- table(factor(x, levels = 0:length(bs)))
pctt.tab <- as.numeric(pctt.tab)/length(x)

y <- data.frame(
    resp = x2
  )
for (i in 1:ncol(pcmpr)) y[[paste('p1',i-1,sep='')]]<-pcmpr[,i]
for (i in 1:ncol(pcmpr)) y[[paste('p2',i-1,sep='')]]<-grmpr[,i]


omega_c <- imv_c(y, pctt.tab, p1 = "p1", p2 = "p2")
omega_t <- imv_t(y, pctt.tab, p1 = "p1", p2 = "p2")

##########################################################
##many items

           
simfun<-function(K,gen) {
    a<-rlnorm(1, .2, .3)
    th<-rnorm(5000)
    ##generate thresholds
    if (gen=='grm') {
        bs <- sort(rnorm(K-1)) ##ordered
    } 
    if (gen=='pcm') {
        bs <- sort(rnorm(K-1))
    }
    ##
    if (gen=='grm') {
        x <- generate_grm_data(th, a, bs)
        x2 <- generate_grm_data(th, a, bs)
    } else {
        x <- generate_pcm_data(th, a, bs)
        x2 <- generate_pcm_data(th, a, bs)
    }
    grmpars<-estimate_grm(x,th)
    pcmpars<-estimate_pcm(x,th)
    ##confirm that these look great https://www.desmos.com/calculator/xev7dsgdqx
    grmpr<-predict_grm(th,grmpars$a,grmpars$b)
    pcmpr<-predict_pcm(th,pcmpars$a,pcmpars$d)
    if (gen=='grm') pr<-grmpr else alt<-grmpr
    if (gen=='pcm') pr<-pcmpr else alt<-pcmpr
    ## Calculate frequency table
    pctt.tab <- table(factor(x, levels = 0:length(bs)))
    pctt.tab <- as.numeric(pctt.tab)/length(x)
    ##
    y <- data.frame(
        resp = x2
    )
    for (i in 1:ncol(pr)) y[[paste('p1',i-1,sep='')]]<-mean(x==(i-1))
    for (i in 1:ncol(pr)) y[[paste('p2',i-1,sep='')]]<-pr[,i]
    omega_c0 <- imv_c(y, pctt.tab, p1 = "p1", p2 = "p2")
    omega_t0 <- imv_t(y, pctt.tab, p1 = "p1", p2 = "p2")
    ##
    y <- data.frame(
        resp = x2
    )
    for (i in 1:ncol(alt)) y[[paste('p1',i-1,sep='')]]<-alt[,i]
    for (i in 1:ncol(pr)) y[[paste('p2',i-1,sep='')]]<-pr[,i]
    omega_c <- imv_c(y, pctt.tab, p1 = "p1", p2 = "p2")
    omega_t <- imv_t(y, pctt.tab, p1 = "p1", p2 = "p2")
    ##
    c(omega_c0,omega_t0,omega_c,omega_t)
}

library(parallel)
nsim<-100
out<-list()
for (ncat in c(3,5,7)) {
    L1<-mclapply(rep(ncat,nsim),simfun,gen='grm',mc.cores=10)
    z1<-do.call("rbind",L1)
    z1<-data.frame(z1)
    z1$K<-ncat
    L2<-mclapply(rep(ncat,nsim),simfun,gen='pcm',mc.cores=10)
    z2<-do.call("rbind",L2)
    z2<-data.frame(z2)
    z2$K<-ncat
    out[[as.character(ncat)]]<-list(grm=colMeans(z1),pcm=colMeans(z2))
}

x<-lapply(out,function(x) x$grm)
x<-do.call("rbind",x)
grm<-data.frame(x)
names(grm)<-c("ctt.c","ctt.t","alt.c","alt.t","K")
x<-lapply(out,function(x) x$pcm)
x<-do.call("rbind",x)
pcm<-data.frame(x)
names(pcm)<-c("ctt.c","ctt.t","alt.c","alt.t","K")

    
par(mfrow=c(2,2),mgp=c(2,1,0),mar=c(3,3,1,1),oma=rep(.5,4))
##
plot(NULL,xlim=c(3,7),xaxt='n',ylab='IMV(CTT,Model)',ylim=c(0,.5),pch=19,xlab='N categories')
abline(h=0)
axis(side=1,c(3,5,7),c('3','5','7'))
legend("bottomleft",bty='n',c("Cat","Thr"),fill=c("black","red"))
mtext(side=3,"DGM: GRM",line=0)
lines(grm$K,grm$ctt.c)
lines(grm$K,grm$ctt.t,col='red')
##
plot(NULL,xlim=c(3,7),xaxt='n',ylab='IMV(CTT,Model)',ylim=c(0,.5),pch=19,xlab='N categories')
abline(h=0)
axis(side=1,c(3,5,7),c('3','5','7'))
legend("bottomleft",bty='n',c("Cat","Thr"),fill=c("black","red"))
mtext(side=3,"DGM: PCM",line=0)
lines(pcm$K,pcm$ctt.c)
lines(pcm$K,pcm$ctt.t,col='red')
##
plot(NULL,xlim=c(3,7),xaxt='n',ylab='IMV(Alt,Model)',ylim=c(-.01,.01),pch=19,xlab='N categories')
abline(h=0)
axis(side=1,c(3,5,7),c('3','5','7'))
legend("bottomleft",bty='n',c("Cat","Thr"),fill=c("black","red"))
mtext(side=3,"DGM: PCM",line=0)
lines(grm$K,grm$alt.c)
lines(grm$K,grm$alt.t,col='red')
##
plot(NULL,xlim=c(3,7),xaxt='n',ylab='IMV(Alt,Model)',ylim=c(-.01,.01),pch=19,xlab='N categories')
abline(h=0)
axis(side=1,c(3,5,7),c('3','5','7'))
legend("bottomleft",bty='n',c("Cat","Thr"),fill=c("black","red"))
mtext(side=3,"DGM: PCM",line=0)
lines(pcm$K,pcm$alt.c)
lines(pcm$K,pcm$alt.t,col='red')






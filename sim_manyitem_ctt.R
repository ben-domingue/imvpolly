source("00funs.R")
library(irtimv)
ntimes<-4
lf<-list.files()

#################################################################
##sim

library(mirt)
ni<-20
np<-2500
out<-list()

for (ncat in c(3,5,7)) {
    for (nsim in 1:25) {
        a <- matrix(rlnorm(ni,.2,.3))
        for (mod in c("graded","gpcm")) {
            if (mod=='graded') {
                                        # for the graded model, ensure that there is enough space between the intercepts,
                                        # otherwise closer categories will not be selected often (minimum distance of 0.3 here)
                diffs <- t(apply(matrix(runif(ni*(ncat-1), .3, 1), ni), 1, cumsum))
                diffs <- -(diffs - rowMeans(diffs))
                d <- diffs + rnorm(ni)
            } 
            if (mod=='gpcm') {
                diffs <- t(apply(matrix(runif(ni*(ncat-1), -0.2, 1), ni), 1, cumsum)) ##allowing for cross-overs
                diffs <- -(diffs - rowMeans(diffs))
                d <- diffs + rnorm(ni)
                d<-cbind(0,d)
            }
            dat <- simdata(a, d, np, itemtype = mod)
            id<-1:nrow(dat)
            L<-list()
            for (i in 1:ncol(dat)) L[[i]]<-data.frame(id=id,item=colnames(dat)[i],resp=dat[,i])
            df<-data.frame(do.call("rbind",L))
            df$item<-paste("item_",df$item,sep='')
            print(table(df$resp))
            out[[paste(mod,ncat,nsim)]]<-df
        }
    }
}

f<-function(df) {
    ni<-length(unique(df$item))
    ncat<-length(unique(df$resp))
    out[[paste(mod,ncat,nsim)]]<-c(ni,ncat,getom(df)) 
}
library(parallel)
out<-mclapply(out,f,mc.cores=10)

z<-data.frame(do.call("rbind",out))
z$mod<-ifelse(grepl("graded",rownames(z)),'graded','gpcm')
l<-split(z,paste(z$mod,z[,2]))
s<-data.frame(t(sapply(l,function(x) colMeans(x[,-ncol(x)]))))
pcm<-s[1:3,]
grm<-s[4:6,]

##
par(mfrow=c(2,2),mgp=c(2,1,0),mar=c(3,3,1,1),oma=rep(.5,4))
plot(NULL,xlim=c(1,3),xaxt='n',ylab='IMV(CTT,Model)',ylim=c(0,.5),pch=19,xlab='N categories')
abline(h=0)
axis(side=1,1:3,c('3','5','7'))
legend("bottomleft",bty='n',c("Cat","Thr"),fill=c("black","red"))
legend("topleft",bty='n',lty=c(1,2),c("GRM","GPCM"))
lines(1:3,grm$ctt.grm)
lines(1:3,grm$ctt.gpcm,lty=2)
lines(1:3,grm$ctt.grm.1,col='red')
lines(1:3,grm$ctt.gpcm.1,lty=2,col='red')
mtext(side=3,"DGM: GRM",line=0)
##
plot(NULL,xlim=c(1,3),xaxt='n',ylab='IMV(CTT,Model)',ylim=c(0,.5),pch=19,xlab='N categories')
abline(h=0)
axis(side=1,1:3,c('3','5','7'))
legend("bottomleft",bty='n',c("Cat","Thr"),fill=c("black","red"))
legend("topleft",bty='n',lty=c(1,2),c("GRM","GPCM"))
lines(1:3,pcm$ctt.grm)
lines(1:3,pcm$ctt.gpcm,lty=2)
lines(1:3,pcm$ctt.grm.1,col='red')
lines(1:3,pcm$ctt.gpcm.1,lty=2,col='red')
mtext(side=3,"DGM: GPCM",line=0)
##################
plot(NULL,xlim=c(1,3),xaxt='n',ylab='IMV(GPCM,GRM)',ylim=c(0,.01),pch=19,xlab='N categories')
abline(h=0)
axis(side=1,1:3,c('3','5','7'))
lines(1:3,grm$gpcm.grm)
lines(1:3,grm$gpcm.grm.1,col='red')
plot(NULL,xlim=c(1,3),xaxt='n',ylab='IMV(GRM,GPCM)',ylim=c(0,.01),pch=19,xlab='N categories')
abline(h=0)
axis(side=1,1:3,c('3','5','7'))
lines(1:3,pcm$grm.gpcm)
lines(1:3,pcm$grm.gpcm.1,col='red')

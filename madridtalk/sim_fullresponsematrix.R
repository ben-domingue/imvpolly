##load functions from irw.R



grm<-function(x,modtype="graded") {
    ##
    resp<-irw::irw_long2resp(x[x$oos==0,])
    library(mirt)
    id<-resp$id
    ni<-ncol(resp)-1
    m<-mirt(resp[,-1],1,modtype)
    th<-fscores(m)
    ##
    id<-resp[,-1]
    stud<-data.frame(id=resp$id,th=th[,1])
    x<-merge(x[x$oos==1,],stud)
    ##
    resp<-resp[,-1]
    maxcat<-max(apply(resp,1,max,na.rm=TRUE))+1
    ##
    L<-split(x,x$item)
    for (j in 1:length(L)) {
        y<-L[[j]]
        y$ncat<-max(y$resp)+1
        iii<-match(unique(y$item),names(coef(m)))
        extr <- extract.item(m, iii)
        y$p <- expected.item(extr, y$th) #min() of first item
        pcat<-probtrace(extr,y$th)
        ##
        for (k in 1:maxcat) {
            if (k>ncol(pcat)) {
                y[[paste("p",k-1,sep='')]]<-NA
            } else {
                y[[paste("p",k-1,sep='')]]<-pcat[,k]
            }
         }
        L[[j]]<-y
    }
    x<-data.frame(do.call("rbind",L))
    x
}

getom<-function(x) {
    x<-x[!is.na(x$resp),]
    ##
    id<-unique(x$id)
    if (length(id)>10000) x<-x[x$id %in% sample(id,10000),]
    ##get rid of lightly used categories
    ll<-split(x,x$item)
    f<-function(y) {
        tab<-table(y$resp)
        nms<-names(tab)[tab<250]
        if (length(nms)>0) {
            for (nm in as.numeric(nms)) {
                if (nm==0) {
                    y$resp<-ifelse(y$resp==nm,1,y$resp)
                } else {
                    y$resp<-ifelse(y$resp==nm,nm-1,y$resp)
                }
            }
            mm<-max(y$resp,na.rm=TRUE)
            vals<-sort(unique(y$resp))
            index<-match(y$resp,vals)
            y$resp<-(0:mm)[index]
        }
        y
    }
    ll<-lapply(ll,f)
    ##lower all categories to 0
    for (i in 1:length(ll)) {
        y<-ll[[i]]
        y$resp<-y$resp-min(y$resp)
        ll[[i]]<-y
    }
    x<-data.frame(do.call("rbind",ll))
    ##
    np<-length(unique(x$id))
    ni<-length(unique(x$item))
    ##
    if (length(unique(by(x$resp,x$item,function(x) max(x))))==1) samecat<-TRUE else samecat<-FALSE
    ##
    x$gr<-sample(1:ntimes,nrow(x),replace=TRUE)
    x.hold<-x
    omega<-list()
    for (i in 1:ntimes) {
        x<-x.hold
        x$oos<-ifelse(x$gr==i,1,0)
        #test<-runif(nrow(x))
        #x$resp<-ifelse(x$resp==4 & test>.5,5,x$resp)
        ##
        xr<-grm(x)
        ##
        x2<-grm(x,modtype="gpcm")
        x2<-x2[,paste("p",0:max(x2$resp,na.rm=TRUE),sep='')]
        names(x2)<-paste("gpcm",0:(ncol(x2)-1),sep='')
        xr<-data.frame(cbind(xr,x2))
        ##
        x2<-grm(x,modtype="nominal")
        x2<-x2[,paste("p",0:max(x2$resp,na.rm=TRUE),sep='')]
        names(x2)<-paste("nom",0:(ncol(x2)-1),sep='')
        xr<-data.frame(cbind(xr,x2))
        ##
        if (samecat) {
            x2<-grm(x,modtype="rsm")
            x2<-x2[,paste("p",0:max(x2$resp,na.rm=TRUE),sep='')]
            names(x2)<-paste("rsm",0:(ncol(x2)-1),sep='')
            xr<-data.frame(cbind(xr,x2))
        }
        ##experimenting
        x0<-x[x$oos==0,]
        L<-split(x0,x0$item)
        pctt.tab<-lapply(L,function(x) table(x$resp)/nrow(x))
        L<-split(xr,xr$item)
        f<-function(y,pctt.tab) {
            nn<-unique(y$ncat)
            for (i in 0:(nn-1)) y[[paste("p0",i,sep='')]]<-pctt.tab[i+1]
            ##
            om0<-imv_c(y,p1="p0",p2="p",pctt.tab) #grm
            om0.gpcm<-imv_c(y,p1="p0",p2="gpcm",pctt.tab) #gpcm
            om<-imv_c(y,p1="p",p2="gpcm",pctt.tab)
            om.alt<-imv_c(y,p1="gpcm",p2="p",pctt.tab)
            omcat<-c(ctt.grm=om0,
                     ctt.gpcm=om0.gpcm,
                     grm.gpcm=om,
                     gpcm.grm=om.alt)
            ##
            om0<-imv_t(y,p1="p0",p2="p",pctt.tab) #grm
            om0.gpcm<-imv_t(y,p1="p0",p2="gpcm",pctt.tab) #gpcm
            om<-imv_t(y,p1="p",p2="gpcm",pctt.tab)
            om.alt<-imv_t(y,p1="gpcm",p2="p",pctt.tab)
            omthr<-c(ctt.grm=om0,
                     ctt.gpcm=om0.gpcm,
                     grm.gpcm=om,
                     gpcm.grm=om.alt)
            ##
            c(omcat,omthr)
        }
        om<-list()
        for (ii in 1:length(L)) om[[ii]]<-f(L[[ii]],pctt.tab[[ii]])
        omega[[i]]<-colMeans(do.call("rbind",om),na.rm=TRUE)
    }
    c(ni,np,colMeans(do.call("rbind",omega)))
}


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

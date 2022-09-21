##looking at the variation in discriminations

parfun<-function(args) {
    library(irtimv)
    sim<-function(np,a=NULL,J=15,K=5,KK) {
        th<-rnorm(np)
        if (is.null(a)) a<-rep(1,J)
        L<-list()
        for (j in 1:J) { #over items
            bb<-seq(-2.5,2.5,length.out=KK-1)
            b<-sort(sample(bb,K-1))
            psi<-list()
            psi[[1]]<-rep(1,length(th))
            for (k in 1:(K-1)) {
                kern<-k*th-sum(b[1:k])
                psi[[k+1]]<-exp(a[j]*kern)
            }
            psi<-do.call("cbind",psi)
            den<-rowSums(psi)
            p<-psi/den
            colnames(p)<-paste("truep",0:(ncol(p)-1),sep='')
            resp<-numeric()
            for (i in 1:np) resp[i]<-which(rmultinom(1,1,p[i,])[,1]>0)-1
            L[[j]]<-data.frame(item=paste0("item_",10+j),id=1:np,resp=resp,p)
        }
        data.frame(do.call("rbind",L))
    }
    ##
    J<-args[[3]]
    np<-args[[2]]
    K<-args[[1]]
    KK<-args[[4]]
    ##
    library(mirt)
    irt<-function(x,mod) { ###re beta prior on guessing param https://groups.google.com/g/mirt-package/c/8Usx53BoXyw
        ##
        resp<-makeresponse(x)
        index<-grep("id",names(resp))
        ni<-ncol(resp)-1
        m<-mirt(resp[,-index],1,mod)
        list(resp,m)
    }
    ##
    x<-sim(np,K=K,J=J,KK=KK)
    ##
    z<-list()
    for (me in c("gpcm","rsm")) {
        m<-irt(x,mod=me)
        if (me=="gpcm") {
            resp<-m[[1]]
            index<-grep("id",colnames(resp))
            id<-resp[,index]
            resp<-resp[,-index]
        }
        th.est<-fscores(m[[2]],response.pattern=resp)
        index<-grep("^F",colnames(th.est))
        th<-data.frame(id=1:nrow(th.est),th=th.est[,index])
        ##
        tmp<-list()
        for (j in unique(x$item)) {
            item<-extract.item(m[[2]],j)
            pr<-probtrace(item,th$th)
            pr<-data.frame(pr)
            names(pr)<-paste(me,0:(K-1),sep='')
            tmp[[j]]<-data.frame(th,pr,item=j)
        }
        tmp<-data.frame(do.call("rbind",tmp))
        tmp$th<-NULL
        z[[me]]<-tmp
    }
    x<-merge(x,z[[1]])
    x<-merge(x,z[[2]])
    ##
    L<-split(x,x$item)
    pctt.tab<-lapply(L,function(x) table(x$resp)/nrow(x))
    ##
    x$resp0<-x$resp
    pr<-x[,paste("truep",0:(K-1),sep='')]
    resp<-numeric()
    for (i in 1:nrow(pr)) resp[i]<-which(rmultinom(1,1,pr[i,])[,1]>0)-1
    x$resp<-resp
    ########################
    ## imv category
    L<-split(x,x$item)
    f<-function(y,pctt.tab) {
        nn<-length(unique(y$resp))
        for (i in 0:(nn-1)) y[[paste("p0",i,sep='')]]<-pctt.tab[i+1]
        om<-imv_c(y,p1="rsm",p2="gpcm",pctt.tab)
        om
    }
    om<-list()
    for (i in 1:length(L)) om[[i]]<-f(L[[i]],pctt.tab[[i]])
    om1<-mean(unlist(om))
    ############################
    ##cumulative imv
    L<-split(x,x$item)
    f<-function(y,pctt.tab) {
        nn<-length(unique(y$resp))
        for (i in 0:(nn-1)) y[[paste("p0",i,sep='')]]<-pctt.tab[i+1]
        om<-imv_t(y,p1="rsm",p2="gpcm",pctt.tab)
        om
    }
    om<-list()
    for (i in 1:length(L)) om[[i]]<-f(L[[i]],pctt.tab[[i]])
    om2<-mean(unlist(om))
    ##
    c(unlist(args),om1,om2)
}


K<-c(3,6)
J<-c(20)
np<-c(5000)
z<-expand.grid(1:25,K,np,J)
argvals<-list()
for (i in 1:nrow(z)) {
    for (KK in c(1,2,4)) argvals[[paste(i,KK)]]<-list(z[i,2],z[i,3],z[i,4],KK*z[i,2])
}

set.seed(867.5309^2)
library(parallel)
tab1<-mclapply(argvals,parfun,mc.cores=2)
tab1<-do.call("rbind",tab1)

tab<-tab1
z<-data.frame(tab)
names(z)<-c("K","n","J","KK","omc","omt")
    
save(z,file="poly_rsm.Rdata")

load("poly_rsm.Rdata")


zz<-split(z,paste(z$n,z$K,z$J,z$KK))
zz<-lapply(zz,colMeans)
z<-data.frame(do.call("rbind",zz))
z$id<-paste(z$n,z$K,z$J,z$KK)
z<-z[order(z$K,z$KK),]


pdf("/home/bd/Dropbox/Apps/Overleaf/IMV_IRT/poly_rating.pdf",width=4,height=3)
par(mgp=c(2,1,0),mar=c(5.5,3,.5,.5),bty='n',oma=rep(.3,4))
plot(NULL,xlim=c(1,nrow(z)+.3),ylim=c(0,.04),xlab="",xaxt="n",ylab="IMV")
for (i in seq(-.025,.025,by=.025)) abline(h=i,lwd=1,col='gray')
axis(side=1,at=1:nrow(z),labels=rep("",nrow(z)))
mtext(side=1,at=1:nrow(z),z$J,line=.75,cex=.8)
mtext(side=1,line=.75,at=0,"J",cex=.8)
mtext(side=1,at=1:nrow(z),z$n,line=1.75,cex=.8)
mtext(side=1,line=1.75,at=0,"N",cex=.8)
mtext(side=1,at=1:nrow(z),z$K,line=2.75,cex=.8)
mtext(side=1,line=2.75,at=0,"K",cex=.8)
mtext(side=1,at=1:nrow(z),z$KK,line=3.75,cex=.8)
mtext(side=1,line=3.75,at=0,"KK",cex=.8)
##
lines(1:nrow(z),z$omc,type='b',col='blue',lwd=2,pch=19)
text(nrow(z),z$omc[nrow(z)],expression(omega[c]),col='blue',pos=4)
lines(1:nrow(z),z$omt,type='b',col='red',lwd=2,pch=19)
text(nrow(z),z$omt[nrow(z)],expression(omega[t]),col='red',pos=4)
dev.off()


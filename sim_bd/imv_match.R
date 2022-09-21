parfun<-function(args) {
    library(irtimv)
    K<-3
    sim<-function(np,a=NULL,J=15,K=5,b.delta) {
        th<-rnorm(np)
        if (is.null(a)) a<-rep(1,J)
        L<-list()
        for (j in 1:J) { #over items
            b<-runif(1,min=-1,max=1)
            b<-c(b,b+b.delta)
            ##
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
    J<-args[[2]]
    np<-args[[1]]
    b.delta<-args[[3]]
    xxx<-sim(np,J=J,b.delta=b.delta,K=K)
    per1<-sum(xxx$resp==1)/nrow(xxx)
    per0<-sum(xxx$resp==0)/nrow(xxx)
    per2<-sum(xxx$resp==2)/nrow(xxx)
    ##     #############################################
    ##dichotomous
    x<-xxx
    x$resp<-ifelse(x$resp>0,1,0)
    #test.p<-runif(nrow(x))
    #x$resp<-ifelse(x$resp==2 | (x$resp==1 & test.p>.5),2,0)
    #x$resp<-ifelse(x$resp==2,1,0)
    ##
    library(mirt)
    irt.dich<-function(x,mod) { ###re beta prior on guessing param https://groups.google.com/g/mirt-package/c/8Usx53BoXyw
        ##
        resp<-makeresponse(x)
        library(mirt)
        index<-grep("id",names(resp))
        ni<-ncol(resp)-1
        if (mod=="Rasch") {
            m<-mirt(resp[,-index],1,"Rasch")
        }
        if (mod=="2PL") {
            s<-paste("F=1-",ni,"\nPRIOR = (1-",ni,", a1, lnorm, 0.2, 0.2)",
                     sep="") 
            model<-mirt.model(s)
            test<-try(m<-mirt(resp[,-index],model,itemtype=rep("2PL",ni),method="EM",technical=list(NCYCLES=10000)))
        }
        if (mod=="3PL") {
            s<-paste("F=1-",ni,"\nPRIOR = (1-",ni,", a1, lnorm, 0.2, 0.2),(1-",ni,", g, expbeta, 2, 17)",
                     sep="") 
            model<-mirt.model(s)
            test<-try(m<-mirt(resp[,-index],model,itemtype=rep(mod,ni),method="EM",technical=list(NCYCLES=10000)))
        }
        list(resp,m)
    }
    z<-list()
    for (me in c("Rasch")) {
        m<-irt.dich(x,mod=me)
        if (me=="Rasch") {
            resp<-m[[1]]
            index<-grep("id",colnames(resp))
            id<-resp[,index]
            resp<-resp[,-index]
        }
        m<-m[[2]]
        co<-coef(m)
        nms<-names(co)
        co<-do.call("rbind",co[-length(co)])
        item<-data.frame(item=nms[-length(nms)],easy=co[,2],load=co[,1],gues=co[,3])
        th.est<-fscores(m,response.pattern=resp)
        index<-grep("^F",colnames(th.est))
        stud<-data.frame(id=id,th=th.est[,index])
        oos<-x
        y<-merge(oos,stud)
        y<-merge(y,item)
        ##
        kk<-y$load*y$th+y$easy
        kk<-exp(kk)
        y$p<-y$gues+(1-y$gues)*kk/(1+kk)
        hold<-y$resp
        truep<-y$truep
        y<-y[,c("item","id","p")]
        if (me=="Rasch") {
            y$resp<-hold
            y$truep<-truep
        }
        names(y)[3]<-me
        z[[me]]<-y                
    }
    x<-merge(x,z[[1]])
    ##
    x$resp0<-x$resp
    x$p<-mean(x$resp0)
    x$truep<-1-x$truep0
    x$resp<-rbinom(nrow(x),1,x$truep)
    ##
    omctt<-imv(x,p1='p',p2='Rasch')
    ## #############################################################
    ##polytomous
    x<-xxx
    irt<-function(x,mod) { ###re beta prior on guessing param https://groups.google.com/g/mirt-package/c/8Usx53BoXyw
        ##
        resp<-makeresponse(x)
        index<-grep("id",names(resp))
        ni<-ncol(resp)-1
        m<-mirt(resp[,-index],1,mod)
        list(resp,m)
    }
    ##
    ##
    z<-list()
    for (me in c("gpcm")) {
        m<-irt(x,mod=me)
        resp<-m[[1]]
        index<-grep("id",colnames(resp))
        id<-resp[,index]
        resp<-resp[,-index]
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
    #x<-merge(x,z[[2]])
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
        om0<-imv_c(y,p1='p0',p2='gpcm',pctt.tab)
        c(om0)
    }
    om<-list()
    for (i in 1:length(L)) om[[i]]<-f(L[[i]],pctt.tab[[i]])
    om<-colMeans(do.call("rbind",om))
    om1<-c(unlist(args),om)
############################
    ##cumulative imv
    L<-split(x,x$item)
    f<-function(y,pctt.tab) {
        nn<-length(unique(y$resp))
        for (i in 0:(nn-1)) y[[paste("p0",i,sep='')]]<-pctt.tab[i+1]
        om0<-imv_t(y,p1='p0',p2='gpcm',pctt.tab)
        c(om0)
    }
    om<-list()
    for (i in 1:length(L)) om[[i]]<-f(L[[i]],pctt.tab[[i]])
    om<-colMeans(do.call("rbind",om))
    om2<-c(unlist(args),om)
    ##
    list(omctt=omctt,
         om.cat=om1,
         om.cum=om2,
         per0=per0,per1=per1,per2=per2
         )
}


#K<-c(3,5,7)
#J<-c(10,20,40)
J<-20
np<-c(5000)
b.delta<- seq(-3,2,by=.2)
##
z<-expand.grid(1:1,np,J,b.delta)
argvals<-list()
for (i in 1:nrow(z)) argvals[[i]]<-list(z[i,2],z[i,3],z[i,4])

library(parallel)
tab1<-mclapply(argvals,parfun,mc.cores=2)

f<-function(x) c(x$om.cum,x$omctt,x$om.cat[4],x$per0,x$per1,x$per2)
z<-do.call("rbind",lapply(tab1,f))
par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(2,1,0))
plot(z[,3],z[,7],type='l',xlim=c(min(b.delta),max(b.delta)*1.1),xlab="b.delta",ylab="% responses",ylim=c(0,1))
text(z[,3],z[,7],'0',pos=4)
lines(z[,3],z[,8],type='l')
text(z[,3],z[,8],'1',pos=4)
lines(z[,3],z[,9],type='l')
text(z[,3],z[,9],'2',pos=4)
##
plot(z[,3],z[,5],type='l',ylim=c(0,.5),xlim=c(min(b.delta),max(b.delta)*1.1),xlab="b.delta",ylab="imv")
nn<-nrow(z)
text(z[nn,3],z[nn,5],'dich',pos=4)
points(z[,3],z[,4],col='red',pch=19)
text(z[nn,3],z[nn,4],'poly-t',pos=4,col='red')
points(z[,3],z[,6],col='blue',pch=19)
text(z[nn,3],z[nn,6],'poly-c',pos=2,col='blue')
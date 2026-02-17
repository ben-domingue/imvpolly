tabs<-irw::irw_filter(n_categories=c(3,6),n_participants=c(1000,25000),density=c(0.8,1))



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
            om<-imv_c(y,p1="p",p2="gpcm",pctt.tab)
            om2<-imv_c(y,p1="gpcm",p2="nom",pctt.tab)
            if (samecat) om3<-imv_c(y,p1="rsm",p2="gpcm",pctt.tab) else om3<-NA
            omcat<-c(ctt_grm_cat=om0,grm_gpcm_cat=om,gpcm_nom_cat=om2,rsm_gpcm_cat=om3)
            ##
            om0<-imv_t(y,p1="p0",p2="p",pctt.tab) #grm
            om<-imv_t(y,p1="p",p2="gpcm",pctt.tab)
            #om2<-imv_t(y,p1="gpcm",p2="nom",pctt.tab)
            if (samecat) om3<-imv_t(y,p1="rsm",p2="gpcm",pctt.tab) else om3<-NA
            omthr<-c(grmthr=om0,gpcmthr=om,
                     #nomthr=om2,  #doesn't make sense?
                     rsmthr=om3)
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
##empirical
##154
out<-list()

for (tab in tabs[154:154]) {
    print(tab)
    fn<-paste(tab,'.Rdata',sep='')
    if (fn %in% lf) {#load(fn) else  df<-irw::irw_fetch(tab)
        load(fn)
        df$item<-paste("item_",df$item,sep='')
        out[[tab]]<-getom(df)
    }
}

tab<-do.call("rbind",out)
write.csv(tab,'')

x<-read.csv("results.csv")
par(mfrow=c(2,1),mgp=c(2,1,0),mar=c(3,3,1,1),oma=rep(.5,4))
plot(x$ctt_grm_cat,x$grmthr,
     xlab='IMV-c(CTT,GRM)',ylab='IMV-t(CTT,GRM)',
     pch=19)
abline(h=0); abline(v=0)
xx<-x[x$grm_gpcm_cat<.05,]
plot(xx$grm_gpcm_cat,xx$gpcmthr,
     xlab='IMV-c(GRM,GPCM)',ylab='IMV-t(GRM,GPCM)',
     pch=19)
abline(h=0); abline(v=0)



z<-x[!is.na(x$rsm_gpcm_cat),]
par(mfrow=c(1,2),mgp=c(2,1,0),mar=c(3,3,1,1),oma=rep(.5,4))
plot(z[,2],z$rsmthr,xlab='N items',ylab="IMV-t",pch=19); abline(h=0)
plot(z[,3],z$rsmthr,xlab='N ppl',ylab="IMV-t",pch=19); abline(h=0)


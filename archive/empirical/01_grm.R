

filenames <-
    list(
        `State Assessment 8` = "cr__c1_2007_8_responses.Rdata",
        `State Assessment 9` = "cr__c1_2007_9_responses.Rdata",
        `State Assessment 7` = "cr__c1_2007_7_responses.Rdata", 
        `State Assessment 6` = "cr__c1_2007_6_responses.Rdata",
        `State Assessment 10` = "cr__c1_2007_10_responses.Rdata", 
        `State Assessment 5` = "cr__c1_2007_5_responses.Rdata",
        `State Assessment 3` = "cr__c1_2007_3_responses.Rdata", 
        `State Assessment 4` = "cr__c1_2007_4_responses.Rdata",
        Openness="ffm_OPN.Rdata",
        Conscientousness="ffm_CSN.Rdata",
        Agreeableness="ffm_AGR.Rdata",
        Stability="ffm_EST.Rdata",
        Extroversion="ffm_EXT.Rdata",
        Grit="grit.Rdata"
         )

library(irtimv)
ntimes<-4



grm<-function(x,modtype="graded") {
    ##
    resp<-makeresponse(x[x$oos==0,])
    library(mirt)
    index<-grep("id",names(resp))
    ni<-ncol(resp)-1
    m<-mirt(resp[,-index],1,modtype)
    th<-fscores(m)
    ##
    id<-resp[,index]
    stud<-data.frame(id=resp$id,th=th[,1])
    x<-merge(x[x$oos==1,],stud)
    ##
    resp<-resp[,-index]
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

out<-list()
ni<-np<-list()
for (fn in filenames) {
    print(fn)
    load(fn)
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
    np[[fn]]<-length(unique(x$id))
    ni[[fn]]<-length(unique(x$item))
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
            omcat<-c(om0,om,om2,om3)
            ##
            om0<-imv_t(y,p1="p0",p2="p",pctt.tab) #grm
            om<-imv_t(y,p1="p",p2="gpcm",pctt.tab)
            om2<-imv_t(y,p1="gpcm",p2="nom",pctt.tab)
            if (samecat) om3<-imv_t(y,p1="rsm",p2="gpcm",pctt.tab) else om3<-NA
            omthr<-c(om0,om,om2,om3)
            ##
            c(omcat,omthr)
        }
        om<-list()
        for (ii in 1:length(L)) om[[ii]]<-f(L[[ii]],pctt.tab[[ii]])
        omega[[i]]<-colMeans(do.call("rbind",om),na.rm=TRUE)
    }
    out[[fn]]<-colMeans(do.call("rbind",omega))
}
tab<-do.call("rbind",out)

tab<-cbind(unlist(np),unlist(ni),tab)
rownames(tab)<-names(filenames)
tab<-rbind(tab,colMeans(tab,na.rm=TRUE))

dump("tab","")

tab <-
structure(c(10000, 10000, 10000, 10000, 10000, 10000, 10000, 
10000, 9950, 9952, 9954, 9953, 9962, 3158, 9494.9285714285706, 
14, 14, 14, 14, 14, 14, 8, 14, 10, 10, 10, 10, 10, 12, 12, 0.14377193041126216, 
0.15251090076319651, 0.15737244851176924, 0.15060393694340976, 
0.14218483579466829, 0.14467853971861983, 0.097704146974378625, 
0.11634318256292792, 0.13070525787387272, 0.18945828706872708, 
0.14639385881503494, 0.27385322138019619, 0.30926206755336177, 
0.22477820079897978, 0.16997291536931464, -0.00025912191485780346, 
0.002768060516226723, 0.00077775521877442163, 0.0020421638043092907, 
0.0019238132607409121, 0.001580473596215885, -0.001618986090883951, 
0.0012030151876041663, -0.010576484454577057, -0.0067158098986547831, 
-0.012981210572007322, -0.0097419260795579658, -0.0096289753798808938, 
-0.0034257909307869091, -0.0031895016955239492, 0.0031696149105267342, 
0.0011958187013096837, 0.0048726260327103463, 0.0032647368281956538, 
0.0029986411222685393, 0.0016002221964944003, 0.004027891312284369, 
0.00072503865818318819, 0.014750724205671752, 0.001817489777030787, 
0.014222027879996898, 0.00039619381502747306, 0.0011685955066624426, 
0.0013841685253512435, 0.0039709849622652504, NaN, NaN, NaN, 
NaN, NaN, NaN, NaN, NaN, 0.0096772535031777168, 0.026441280372597607, 
0.033930571709907917, 0.023088740958376305, 0.0099373206517183003, 
NaN, 0.020615033439155571, 0.15610519068044376, 0.1478819250669349, 
0.1662438355906759, 0.17349307136962852, 0.13385612127658708, 
0.15947004743881565, 0.13167050581962825, 0.11757502374451989, 
0.15677257182647875, 0.14808993354287039, 0.1374357819575468, 
0.18135918178271088, 0.21440462700340568, 0.15282673375600095, 
0.15551318220401766, 0.00011389432303993624, 0.0013232676442330525, 
0.00044475761828878117, 0.00025407059098814868, 0.00055356317003927972, 
0.00088195964623404759, 0.00072713610973471467, 0.00077694568086579985, 
-0.0058619438305907096, -0.00091559921268068742, -0.0048239939634059916, 
-0.00035653032782064214, -0.0011249958504374832, 0.0004713038312321408, 
-0.00053829746930568662, 0.0014395065420344964, 0.00043768504677363346, 
0.0014582739991355323, 0.0014956836571195157, 0.0013221465814800722, 
0.00061031468108967233, 0.0013184173967052267, 0.00074422760841317431, 
0.010010602759017382, 0.0015122734094018407, 0.0067074670323574032, 
0.00059111289728581106, 0.00068744116729063112, -5.1578791340237628e-05, 
0.0020202552847688681, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 
0.0068868099877768037, 0.0068501719625131737, 0.020928584970446956, 
0.014538011753311006, 0.0033761144326694193, NaN, 0.010515938621343472
), .Dim = c(15L, 10L), .Dimnames = list(c("State Assessment 8", 
"State Assessment 9", "State Assessment 7", "State Assessment 6", 
"State Assessment 10", "State Assessment 5", "State Assessment 3", 
"State Assessment 4", "Openness", "Conscientousness", "Agreeableness", 
"Stability", "Extroversion", "Grit", ""), NULL))


#write.csv(tab,'')
library(xtable)
print(xtable(tab,digits=1,display=c("f","d","d",
                                    "e","e","e","e",
                                    "e","e","e","e")
             ))





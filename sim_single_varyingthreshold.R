library(irtimv)
source("00funs.R")

##########################################################
##one example calculation
b0<-0
b1<-1
a<-1
set.seed(123)
th<-rnorm(10000)

x <- generate_grm_data(th, a, b0, b1)
x2 <- generate_grm_data(th, a, b0, b1)

grmpars<-estimate_grm(x,th)
pcmpars<-estimate_pcm(x,th)
##confirm that these look great https://www.desmos.com/calculator/xev7dsgdqx

grmpr<-predict_grm(th,grmpars$a,grmpars$b0,grmpars$b1)
pcmpr<-predict_pcm(th,pcmpars$a,pcmpars$d1,pcmpars$d2)

## Calculate frequency table
pctt.tab <- table(factor(x, levels = 0:2))
pctt.tab <- as.numeric(pctt.tab)/length(x)

y <- data.frame(
    resp = x2,
    p10 = pcmpr[,1],
    p11 = pcmpr[,2],
    p12 = pcmpr[,3],
    p20 = grmpr[,1],
    p21 = grmpr[,2],
    p22 = grmpr[,3]
  )

omega_c <- imv_c(y, pctt.tab, p1 = "p1", p2 = "p2")
omega_t <- imv_t(y, pctt.tab, p1 = "p1", p2 = "p2")

##########################################################
##many items
simfun<-function(b1,gen) {
    b0<-0
    a<-1
    set.seed(123)
    th<-rnorm(5000)
    if (gen=='grm') {
        x <- generate_grm_data(th, a, b0, b1)
        x2 <- generate_grm_data(th, a, b0, b1)
    } else {
        x <- generate_pcm_data(th, a, b0, b1)
        x2 <- generate_pcm_data(th, a, b0, b1)
    }
    grmpars<-estimate_grm(x,th)
    pcmpars<-estimate_pcm(x,th)
    ##confirm that these look great https://www.desmos.com/calculator/xev7dsgdqx
    grmpr<-predict_grm(th,grmpars$a,grmpars$b0,grmpars$b1)
    pcmpr<-predict_pcm(th,pcmpars$a,pcmpars$d1,pcmpars$d2)
    if (gen=='grm') pr<-grmpr else alt<-grmpr
    if (gen=='pcm') pr<-pcmpr else alt<-pcmpr
    ## Calculate frequency table
    pctt.tab <- table(factor(x, levels = 0:2))
    pctt.tab <- as.numeric(pctt.tab)
    y <- data.frame(
        resp = x2,
        p10 = mean(x==0),
        p11 = mean(x==1),
        p12 = mean(x==2),
        p20 = pr[,1],
        p21 = pr[,2],
        p22 = pr[,3]
    )
    omega_c0 <- imv_c(y, pctt.tab, p1 = "p1", p2 = "p2")
    omega_t0 <- imv_t(y, pctt.tab, p1 = "p1", p2 = "p2")
    ##
    y <- data.frame(
        resp = x2,
        p10 = alt[,1],
        p11 = alt[,2],
        p12 = alt[,3],
        p20 = pr[,1],
        p21 = pr[,2],
        p22 = pr[,3]
    )
    omega_c <- imv_c(y, pctt.tab, p1 = "p1", p2 = "p2")
    omega_t <- imv_t(y, pctt.tab, p1 = "p1", p2 = "p2")
    ##
    c(omega_c0,omega_t0,omega_c,omega_t)
}

b1<-sort(runif(100,0,1.5))
L1<-lapply(b1,simfun,gen='grm')
z1<-do.call("rbind",L1)
z1<-data.frame(z1)
z1$b1<-b1
b1<-sort(runif(100,-1,1.5))
L2<-lapply(b1,simfun,gen='pcm')
z2<-do.call("rbind",L2)
z2<-data.frame(z2)
z2$b1<-b1


pf<-function(x,y,...) {
    m<-loess(y~x)
    lines(x,predict(m),...)
}
par(mfrow=c(2,2),mgp=c(2,1,0),mar=c(3,3,2,1),oma=rep(.5,4))
L<-list(grm=z1,pcm=z2)
for (i in 1:length(L)) {
    z<-L[[i]]
    b1<-z$b1
    plot(NULL,xlim=range(b1),ylim=c(-.005,.5),xlab='b1',ylab='imv(mean,model)')
    mtext(side=3,line=0,names(L)[i])
    abline(h=0)
    pf(b1,z[,1])
    pf(b1,z[,2],col='red')
    plot(NULL,xlim=range(b1),ylim=c(-.005,.005),xlab='b1',ylab='imv(alt,model)')
    abline(h=0)
    pf(b1,z[,3])
    pf(b1,z[,4],col='red')
}
legend("topright",bty='n',fill=c("black","red"),c("category","threshold"))



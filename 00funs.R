
# ============================================================================
# FUNCTION 1: Generate data from Graded Response Model
# ============================================================================
generate_grm_data <- function(th, a, b0, b1) {
  # th: vector of abilities
  # a: discrimination parameter
  # b0, b1: threshold parameters (b0 < b1)
  
  n <- length(th)
  
  # Calculate cumulative probabilities
  # P(X >= 1 | th)
  p1_star <- 1 / (1 + exp(-a * (th - b0)))
  # P(X >= 2 | th)
  p2_star <- 1 / (1 + exp(-a * (th - b1)))
  
  # Calculate category probabilities
  # P(X = 0)
  p0 <- 1 - p1_star
  # P(X = 1)
  p1 <- p1_star - p2_star
  # P(X = 2)
  p2 <- p2_star
  
  # Generate responses
  x <- numeric(n)
  for(i in 1:n) {
    x[i] <- sample(0:2, size = 1, prob = c(p0[i], p1[i], p2[i]))
  }
  
  return(x)
}

# ============================================================================
# FUNCTION 2: Estimate GRM parameters given theta
# ============================================================================
estimate_grm <- function(x, th) {
  # x: response vector (0, 1, 2)
  # th: known ability vector
  
  # Negative log-likelihood for GRM
  nll_grm <- function(par) {
    a <- exp(par[1])  # Use exp to ensure a > 0
    b0 <- par[2]
    b1 <- par[3]
    
    if(b1 <= b0) return(1e10)  # Ensure b1 > b0
    
    # Calculate cumulative probabilities
    p1_star <- 1 / (1 + exp(-a * (th - b0)))
    p2_star <- 1 / (1 + exp(-a * (th - b1)))
    
    # Category probabilities
    p0 <- 1 - p1_star
    p1 <- p1_star - p2_star
    p2 <- p2_star
    
    # Ensure probabilities are valid
    p0 <- pmax(p0, 1e-10)
    p1 <- pmax(p1, 1e-10)
    p2 <- pmax(p2, 1e-10)
    
    # Log-likelihood
    ll <- sum((x == 0) * log(p0) + (x == 1) * log(p1) + (x == 2) * log(p2))
    
    return(-ll)
  }
  
  # Optimize
  fit <- optim(par = c(0, -0.5, 0.5), fn = nll_grm, method = "BFGS")
  
  a_hat <- exp(fit$par[1])
  b0_hat <- fit$par[2]
  b1_hat <- fit$par[3]
  
  return(list(a = a_hat, b0 = b0_hat, b1 = b1_hat))
}

# ============================================================================
# FUNCTION 3: Estimate 2PL parameters given theta
# ============================================================================
estimate_2pl <- function(x_binary, th) {
  # x_binary: binary response vector (0, 1)
  # th: known ability vector
  
  # Negative log-likelihood for 2PL
  nll_2pl <- function(par) {
    a <- exp(par[1])  # Use exp to ensure a > 0
    b <- par[2]
    
    # Probability of correct response
    p <- 1 / (1 + exp(-a * (th - b)))
    p <- pmax(pmin(p, 1 - 1e-10), 1e-10)  # Bound probabilities
    
    # Log-likelihood
    ll <- sum(x_binary * log(p) + (1 - x_binary) * log(1 - p))
    
    return(-ll)
  }
  
  # Optimize
  fit <- optim(par = c(0, 0), fn = nll_2pl, method = "BFGS")
  
  a_hat <- exp(fit$par[1])
  b_hat <- fit$par[2]
  
  return(list(a = a_hat, b = b_hat))
}

# ============================================================================
# FUNCTION 4: Predict probabilities from GRM
# ============================================================================
predict_grm <- function(th, a, b0, b1) {
  # Returns matrix with columns: p0, p1, p2
  
  p1_star <- 1 / (1 + exp(-a * (th - b0)))
  p2_star <- 1 / (1 + exp(-a * (th - b1)))
  
  p0 <- 1 - p1_star
  p1 <- p1_star - p2_star
  p2 <- p2_star
  
  return(data.frame(p0 = p0, p1 = p1, p2 = p2))
}

# ============================================================================
# FUNCTION 5: Predict probabilities from 2PL
# ============================================================================
predict_2pl <- function(th, a, b) {
  # Returns probability of response = 1
  p <- 1 / (1 + exp(-a * (th - b)))
  return(p)
}



generate_pcm_data <- function(th, a, d1, d2) {
  # th: vector of abilities
  # a: discrimination parameter
  # d1, d2: step difficulty parameters (d0 = 0 for identification)
  
  n <- length(th)
  
  # Calculate unnormalized probabilities
  p0_num <- rep(1, n)
  p1_num <- exp(a * (th - d1))
  p2_num <- exp(a * (2 * th - d1 - d2))
  
  # Denominator
  denom <- p0_num + p1_num + p2_num
  
  # Category probabilities
  p0 <- p0_num / denom
  p1 <- p1_num / denom
  p2 <- p2_num / denom
  
  # Generate responses
  x <- numeric(n)
  for(i in 1:n) {
    x[i] <- sample(0:2, size = 1, prob = c(p0[i], p1[i], p2[i]))
  }
  
  return(x)
}
estimate_pcm <- function(x, th) {
  # x: response vector (0, 1, 2)
  # th: known ability vector
  
  # Negative log-likelihood for PCM
  nll_pcm <- function(par) {
    a <- exp(par[1])  # Use exp to ensure a > 0
    d0 <- 0           # First step difficulty (set to 0 for identification)
    d1 <- par[2]      # Second step difficulty
    d2 <- par[3]      # Third step difficulty
    
    # Calculate unnormalized probabilities (numerators)
    # P*(X = 0) = exp(0)
    p0_num <- rep(1, length(th))
    
    # P*(X = 1) = exp(a * (th - d1))
    p1_num <- exp(a * (th - d1))
    
    # P*(X = 2) = exp(a * (2*th - d1 - d2))
    p2_num <- exp(a * (2 * th - d1 - d2))
    
    # Denominator (partition function)
    denom <- p0_num + p1_num + p2_num
    
    # Category probabilities
    p0 <- p0_num / denom
    p1 <- p1_num / denom
    p2 <- p2_num / denom
    
    # Ensure probabilities are valid
    p0 <- pmax(p0, 1e-10)
    p1 <- pmax(p1, 1e-10)
    p2 <- pmax(p2, 1e-10)
    
    # Log-likelihood
    ll <- sum((x == 0) * log(p0) + (x == 1) * log(p1) + (x == 2) * log(p2))
    
    return(-ll)
  }
  
  # Optimize
  fit <- optim(par = c(0, 0, 0.5), fn = nll_pcm, method = "BFGS")
  
  a_hat <- exp(fit$par[1])
  d0_hat <- 0  # Fixed for identification
  d1_hat <- fit$par[2]
  d2_hat <- fit$par[3]
  
  return(list(a = a_hat, d0 = d0_hat, d1 = d1_hat, d2 = d2_hat))
}

# ============================================================================
# FUNCTION: Predict probabilities from PCM
# ============================================================================
predict_pcm <- function(th, a, d1, d2) {
  # Returns dataframe with columns: p0, p1, p2
  
  # Calculate unnormalized probabilities
  p0_num <- rep(1, length(th))
  p1_num <- exp(a * (th - d1))
  p2_num <- exp(a * (2 * th - d1 - d2))
  
  # Denominator
  denom <- p0_num + p1_num + p2_num
  
  # Category probabilities
  p0 <- p0_num / denom
  p1 <- p1_num / denom
  p2 <- p2_num / denom
  
  return(data.frame(p0 = p0, p1 = p1, p2 = p2))
}



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

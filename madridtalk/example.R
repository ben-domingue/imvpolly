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
# FUNCTION: Estimate PCM parameters given theta
# ============================================================================
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
pctt.tab <- as.numeric(pctt.tab)

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



## ##deprecated
## simfun<-function(b1,gen) {
##     b0<-0
##     a<-1
##     set.seed(123)
##     th<-rnorm(5000)
##     if (gen=='grm') {
##         x <- generate_grm_data(th, a, b0, b1)
##         x2 <- generate_grm_data(th, a, b0, b1)
##     } else {
##         x <- generate_pcm_data(th, a, b0, b1)
##         x2 <- generate_pcm_data(th, a, b0, b1)
##     }
##     grmpars<-estimate_grm(x,th)
##     pcmpars<-estimate_pcm(x,th)
##     ##confirm that these look great https://www.desmos.com/calculator/xev7dsgdqx
##     grmpr<-predict_grm(th,grmpars$a,grmpars$b0,grmpars$b1)
##     pcmpr<-predict_pcm(th,pcmpars$a,pcmpars$d1,pcmpars$d2)
##     if (gen=='grm') pr<-grmpr else alt<-grmpr
##     if (gen=='pcm') pr<-pcmpr else alt<-pcmpr
##     ## Calculate frequency table
##     pctt.tab <- table(factor(x, levels = 0:2))
##     pctt.tab <- as.numeric(pctt.tab)
##     y <- data.frame(
##         resp = x2,
##         p10 = mean(x==0),
##         p11 = mean(x==1),
##         p12 = mean(x==2),
##         p20 = pr[,1],
##         p21 = pr[,2],
##         p22 = pr[,3]
##     )
##     omega_c0 <- imv_c(y, pctt.tab, p1 = "p1", p2 = "p2")
##     omega_t0 <- imv_t(y, pctt.tab, p1 = "p1", p2 = "p2")
##     ##
##     y <- data.frame(
##         resp = x2,
##         p10 = alt[,1],
##         p11 = alt[,2],
##         p12 = alt[,3],
##         p20 = pr[,1],
##         p21 = pr[,2],
##         p22 = pr[,3]
##     )
##     omega_c <- imv_c(y, pctt.tab, p1 = "p1", p2 = "p2")
##     omega_t <- imv_t(y, pctt.tab, p1 = "p1", p2 = "p2")
##     ##
##     c(omega_c0,omega_t0,omega_c,omega_t)
## }

## b1<-sort(runif(1000,0,1.5))
## L1<-lapply(b1,simfun,gen='grm')
## z1<-do.call("rbind",L1)
## z1<-data.frame(z1)
## z1$b1<-b1
## b1<-sort(runif(1000,-1,1.5))
## L2<-lapply(b1,simfun,gen='pcm')
## z2<-do.call("rbind",L2)
## z2<-data.frame(z2)
## z2$b1<-b1


## pf<-function(x,y,...) {
##     m<-loess(y~x)
##     lines(x,predict(m),...)
## }
## par(mfrow=c(2,2),mgp=c(2,1,0),mar=c(3,3,2,1),oma=rep(.5,4))
## L<-list(grm=z1,pcm=z2)
## for (i in 1:length(L)) {
##     z<-L[[i]]
##     b1<-z$b1
##     plot(NULL,xlim=range(b1),ylim=c(-.005,.5),xlab='b1',ylab='imv(mean,model)')
##     mtext(side=3,line=0,names(L)[i])
##     abline(h=0)
##     pf(b1,z[,1])
##     pf(b1,z[,2],col='red')
##     plot(NULL,xlim=range(b1),ylim=c(-.005,.005),xlab='b1',ylab='imv(alt,model)')
##     abline(h=0)
##     pf(b1,z[,3])
##     pf(b1,z[,4],col='red')
## }
## legend("topright",bty='n',fill=c("black","red"),c("category","threshold"))





                                        # Load required library
library(mirt)
library(irtimv)

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

# ============================================================================
# FUNCTION 6: Run single simulation iteration
# ============================================================================
run_single_simulation <- function(th, a, b0, b1) {
  
  # Step 1: Generate initial data x
  x <- generate_grm_data(th, a, b0, b1)
  
  # Step 2: Fit models to x
  # Fit GRM
  grm_fit <- estimate_grm(x, th)
  
  # Dichotomize and fit 2PLs
  x_prime <- as.numeric(x > 0)
  x_double <- as.numeric(x > 1)
  
  fit_2pl_prime <- estimate_2pl(x_prime, th)
  fit_2pl_double <- estimate_2pl(x_double, th)
  
  # Step 3: Generate new test data x2
  x2 <- generate_grm_data(th, a, b0, b1)
  
  # Step 4: Get predictions for x2
  # Predictions from GRM (for polytomous response)
  pred_grm <- predict_grm(th, grm_fit$a, grm_fit$b0, grm_fit$b1)
  
  # Step 5: Calculate omega for dichotomizations
  # For x' dichotomization (x > 0)
  x2_prime <- as.numeric(x2 > 0)
  p_prime_2pl <- predict_2pl(th, fit_2pl_prime$a, fit_2pl_prime$b)
  
  z_prime <- data.frame(
    resp = x2_prime,
    p1 = mean(x),
    p2 = p_prime_2pl
  )
  omega_x_prime <- imv(z_prime, p1 = "p1", p2 = "p2")
  
  # For x'' dichotomization (x > 1)
  x2_double <- as.numeric(x2 > 1)
  p_double_2pl <- predict_2pl(th, fit_2pl_double$a, fit_2pl_double$b)
  
  z_double <- data.frame(
    resp = x2_double,
    p1 = mean(x),
    p2 = p_double_2pl
  )
  omega_x_double <- imv(z_double, p1 = "p1", p2 = "p2")
  
  # Step 6: Calculate omega_c and omega_t
  # Reconstruct polytomous predictions from the two 2PLs
  
  # Create full prediction dataframe for imv_c and imv_t
  y <- data.frame(
    resp = x2,
    p10 = mean(x==0),
    p11 = mean(x==1),
    p12 = mean(x==2),
    p20 = pred_grm$p0,
    p21 = pred_grm$p1,
    p22 = pred_grm$p2
  )
  
  # Calculate frequency table
  pctt.tab <- table(factor(x2, levels = 0:2))
  pctt.tab <- as.numeric(pctt.tab)
  
  omega_c <- imv_c(y, pctt.tab, p1 = "p1", p2 = "p2")
  omega_t <- imv_t(y, pctt.tab, p1 = "p1", p2 = "p2")
  
  # Return results
  return(data.frame(
    b1_true = b1,
    omega_x_prime = omega_x_prime,
    omega_x_double = omega_x_double,
    omega_c = omega_c,
    omega_t = omega_t
  ))
}

# ============================================================================
# FUNCTION 7: Main simulation
# ============================================================================
run_simulation <- function(n_theta = 1000, n_sim = 10, seed = 123,b0=0,mu=0) {
  
  set.seed(seed)
  
  # Fixed parameters
  a <- 1
  
  # Generate theta values (abilities)
  th <- rnorm(n_theta,mean=mu)
  
  # Sample b1 values
  b1_values <- sort(runif(n_sim, b0, 1))
  
  # Run simulations
  results <- list()
  
  for(i in 1:n_sim) {
    cat("Running simulation", i, "of", n_sim, "with b1 =", 
        round(b1_values[i], 3), "\n")
    
    results[[i]] <- run_single_simulation(th, a, b0, b1_values[i])
  }
  
  # Combine results
  results_df <- do.call(rbind, results)
  results_df$b0<-b0
  return(results_df)
}

# ============================================================================
# RUN THE SIMULATION
# ============================================================================

# Run with default settings
results<-list()
for (b0 in seq(-1,.5,by=.5)) results[[as.character(b0)]] <- run_simulation(n_theta = 1000, n_sim = 100, seed = 123,b0=b0,mu=b0)
results2<-list()
for (b0 in seq(-1,.5,by=.5)) results2[[as.character(b0)]] <- run_simulation(n_theta = 1000, n_sim = 100, seed = 123,b0=b0,mu=0)


f<-function(results) {
    pf<-function(x,y,...) {
        m<-loess(y~x)
        lines(x,predict(m),...)
    }
    plot(NULL,xlim=c(-1,1),ylim=c(0,.75),xlab='second threshold',ylab='IMV')
                                        #pf(results$b1_true,results$omega_c)
    pf(results$b1_true,results$omega_t,lwd=2)
    pf(results$b1_true,results$omega_x_prime,col='red',lwd=2)
    pf(results$b1_true,results$omega_x_double,col='blue',lwd=2)
    pf(results$b1_true,results$omega_c,lty=2)
    abline(v=unique(results$b0))
    }
par(mfcol=c(4,1),mgp=c(2,1,0),mar=c(3,3,1,1),oma=rep(.5,4))
for (i in 1:length(results)) f(results[[i]])
legend("left",bty='n',fill=c("black","red","blue"),c("grm","2pl_12","2pl_2"))

#for (i in 1:length(results2)) f(results2[[i]])



## a figure
th<-seq(-4,4,length.out=1000)
par(mfrow=c(5,1),mgp=c(2,1,0),mar=c(3,3,1,1),oma=rep(.5,4))
for (b1 in seq(1.5,-1+.05,length.out=5)) {
    z<-predict_grm(th,a=1,b=-1,b1=b1)
    plot(NULL,xlim=range(th),ylim=c(0,1),xlab=expression(theta),ylab="CRF")
    cols<-c("black","red","black")
    for (i in 1:ncol(z)) lines(th,z[,i],col=cols[i])
    abline(v=-1)
    abline(v=b1)
}

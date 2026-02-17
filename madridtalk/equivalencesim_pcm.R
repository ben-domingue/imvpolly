##i don't actually think this makes sense as pcm never collapses




## # Load required library
## library(mirt)
## library(irtimv)


## generate_pcm_data <- function(th, a, d1, d2) {
##   # th: vector of abilities
##   # a: discrimination parameter
##   # d1, d2: step difficulty parameters (d0 = 0 for identification)
  
##   n <- length(th)
  
##   # Calculate unnormalized probabilities
##   p0_num <- rep(1, n)
##   p1_num <- exp(a * (th - d1))
##   p2_num <- exp(a * (2 * th - d1 - d2))
  
##   # Denominator
##   denom <- p0_num + p1_num + p2_num
  
##   # Category probabilities
##   p0 <- p0_num / denom
##   p1 <- p1_num / denom
##   p2 <- p2_num / denom
  
##   # Generate responses
##   x <- numeric(n)
##   for(i in 1:n) {
##     x[i] <- sample(0:2, size = 1, prob = c(p0[i], p1[i], p2[i]))
##   }
  
##   return(x)
## }
## estimate_pcm <- function(x, th) {
##   # x: response vector (0, 1, 2)
##   # th: known ability vector
  
##   # Negative log-likelihood for PCM
##   nll_pcm <- function(par) {
##     a <- exp(par[1])  # Use exp to ensure a > 0
##     d0 <- 0           # First step difficulty (set to 0 for identification)
##     d1 <- par[2]      # Second step difficulty
##     d2 <- par[3]      # Third step difficulty
    
##     # Calculate unnormalized probabilities (numerators)
##     # P*(X = 0) = exp(0)
##     p0_num <- rep(1, length(th))
    
##     # P*(X = 1) = exp(a * (th - d1))
##     p1_num <- exp(a * (th - d1))
    
##     # P*(X = 2) = exp(a * (2*th - d1 - d2))
##     p2_num <- exp(a * (2 * th - d1 - d2))
    
##     # Denominator (partition function)
##     denom <- p0_num + p1_num + p2_num
    
##     # Category probabilities
##     p0 <- p0_num / denom
##     p1 <- p1_num / denom
##     p2 <- p2_num / denom
    
##     # Ensure probabilities are valid
##     p0 <- pmax(p0, 1e-10)
##     p1 <- pmax(p1, 1e-10)
##     p2 <- pmax(p2, 1e-10)
    
##     # Log-likelihood
##     ll <- sum((x == 0) * log(p0) + (x == 1) * log(p1) + (x == 2) * log(p2))
    
##     return(-ll)
##   }
  
##   # Optimize
##   fit <- optim(par = c(0, 0, 0.5), fn = nll_pcm, method = "BFGS")
  
##   a_hat <- exp(fit$par[1])
##   d0_hat <- 0  # Fixed for identification
##   d1_hat <- fit$par[2]
##   d2_hat <- fit$par[3]
  
##   return(list(a = a_hat, d0 = d0_hat, d1 = d1_hat, d2 = d2_hat))
## }

## # ============================================================================
## # FUNCTION: Predict probabilities from PCM
## # ============================================================================
## predict_pcm <- function(th, a, d1, d2) {
##   # Returns dataframe with columns: p0, p1, p2
  
##   # Calculate unnormalized probabilities
##   p0_num <- rep(1, length(th))
##   p1_num <- exp(a * (th - d1))
##   p2_num <- exp(a * (2 * th - d1 - d2))
  
##   # Denominator
##   denom <- p0_num + p1_num + p2_num
  
##   # Category probabilities
##   p0 <- p0_num / denom
##   p1 <- p1_num / denom
##   p2 <- p2_num / denom
  
##   return(data.frame(p0 = p0, p1 = p1, p2 = p2))
## }


## # ============================================================================
## # FUNCTION 3: Estimate 2PL parameters given theta
## # ============================================================================
## estimate_2pl <- function(x_binary, th) {
##   # x_binary: binary response vector (0, 1)
##   # th: known ability vector
  
##   # Negative log-likelihood for 2PL
##   nll_2pl <- function(par) {
##     a <- exp(par[1])  # Use exp to ensure a > 0
##     b <- par[2]
    
##     # Probability of correct response
##     p <- 1 / (1 + exp(-a * (th - b)))
##     p <- pmax(pmin(p, 1 - 1e-10), 1e-10)  # Bound probabilities
    
##     # Log-likelihood
##     ll <- sum(x_binary * log(p) + (1 - x_binary) * log(1 - p))
    
##     return(-ll)
##   }
  
##   # Optimize
##   fit <- optim(par = c(0, 0), fn = nll_2pl, method = "BFGS")
  
##   a_hat <- exp(fit$par[1])
##   b_hat <- fit$par[2]
  
##   return(list(a = a_hat, b = b_hat))
## }

## # ============================================================================
## # FUNCTION 5: Predict probabilities from 2PL
## # ============================================================================
## predict_2pl <- function(th, a, b) {
##   # Returns probability of response = 1
##   p <- 1 / (1 + exp(-a * (th - b)))
##   return(p)
## }

## # ============================================================================
## # FUNCTION 6: Run single simulation iteration
## # ============================================================================
## run_single_simulation <- function(th, a, b0, b1) {
  
##   # Step 1: Generate initial data x
##   x <- generate_pcm_data(th, a, b0, b1)
  
##   # Step 2: Fit models to x
##   # Fit PCM
##   pcm_fit <- estimate_pcm(x, th)
  
##   # Dichotomize and fit 2PLs
##   x_prime <- as.numeric(x > 0)
##   x_double <- as.numeric(x > 1)
  
##   fit_2pl_prime <- estimate_2pl(x_prime, th)
##   fit_2pl_double <- estimate_2pl(x_double, th)
  
##   # Step 3: Generate new test data x2
##   x2 <- generate_pcm_data(th, a, b0, b1)
  
##   # Step 4: Get predictions for x2
##   # Predictions from PCM (for polytomous response)
##   pred_pcm <- predict_pcm(th, pcm_fit$a, pcm_fit$d1, pcm_fit$d2)
  
##   # Step 5: Calculate omega for dichotomizations
##   # For x' dichotomization (x > 0)
##   x2_prime <- as.numeric(x2 > 0)
##   p_prime_2pl <- predict_2pl(th, fit_2pl_prime$a, fit_2pl_prime$b)
  
##   z_prime <- data.frame(
##     resp = x2_prime,
##     p1 = mean(x),
##     p2 = p_prime_2pl
##   )
##   omega_x_prime <- imv(z_prime, p1 = "p1", p2 = "p2")
  
##   # For x'' dichotomization (x > 1)
##   x2_double <- as.numeric(x2 > 1)
##   p_double_2pl <- predict_2pl(th, fit_2pl_double$a, fit_2pl_double$b)
  
##   z_double <- data.frame(
##     resp = x2_double,
##     p1 = mean(x),
##     p2 = p_double_2pl
##   )
##   omega_x_double <- imv(z_double, p1 = "p1", p2 = "p2")
  
##   # Step 6: Calculate omega_c and omega_t
##   # Reconstruct polytomous predictions from the two 2PLs
  
##   # Create full prediction dataframe for imv_c and imv_t
##   y <- data.frame(
##     resp = x2,
##     p10 = mean(x==0),
##     p11 = mean(x==1),
##     p12 = mean(x==2),
##     p20 = pred_pcm$p0,
##     p21 = pred_pcm$p1,
##     p22 = pred_pcm$p2
##   )
  
##   # Calculate frequency table
##   pctt.tab <- table(factor(x2, levels = 0:2))
##   pctt.tab <- as.numeric(pctt.tab)
  
##   omega_c <- imv_c(y, pctt.tab, p1 = "p1", p2 = "p2")
##   omega_t <- imv_t(y, pctt.tab, p1 = "p1", p2 = "p2")
  
##   # Return results
##   return(data.frame(
##     b1_true = b1,
##     omega_x_prime = omega_x_prime,
##     omega_x_double = omega_x_double,
##     omega_c = omega_c,
##     omega_t = omega_t
##   ))
## }

## # ============================================================================
## # FUNCTION 7: Main simulation
## # ============================================================================
## run_simulation <- function(n_theta = 1000, n_sim = 10, seed = 123,b0=0) {
  
##   set.seed(seed)
  
##   # Fixed parameters
##   a <- 1
  
##   # Generate theta values (abilities)
##   th <- rnorm(n_theta)
  
##   # Sample b1 values
##   b1_values <- sort(runif(n_sim, -1, 1))
  
##   # Run simulations
##   results <- list()
  
##   for(i in 1:n_sim) {
##     cat("Running simulation", i, "of", n_sim, "with b1 =", 
##         round(b1_values[i], 3), "\n")
    
##     results[[i]] <- run_single_simulation(th, a, b0, b1_values[i])
##   }
  
##   # Combine results
##   results_df <- do.call(rbind, results)
##   results_df$b0<-b0
##   return(results_df)
## }

## # ============================================================================
## # RUN THE SIMULATION
## # ============================================================================

## # Run with default settings
## results<-list()
## for (b0 in seq(-1,.5,by=.5)) results[[as.character(b0)]] <- run_simulation(n_theta = 1000, n_sim = 100, seed = 123)


## f<-function(results) {
##     pf<-function(x,y,...) {
##         m<-loess(y~x)
##         lines(x,predict(m),...)
##     }
##     plot(NULL,xlim=c(-1,1),ylim=c(0,.75))
##                                         #pf(results$b1_true,results$omega_c)
##     pf(results$b1_true,results$omega_t,lwd=2)
##     pf(results$b1_true,results$omega_x_prime,col='red',lwd=2)
##     pf(results$b1_true,results$omega_x_double,col='blue',lwd=2)
##     pf(results$b1_true,results$omega_c,lty=2)
##     abline(v=unique(results$b0))
##     }
## par(mfrow=c(4,1),mgp=c(2,1,0),mar=c(3,3,1,1),oma=rep(.1,4))
## for (i in 1:length(results)) f(results[[i]])
## legend("left",bty='n',fill=c("black","red","blue"),c("pcm","2pl_12","2pl_2"))

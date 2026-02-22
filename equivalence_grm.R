##will need to fix the way `y` is defined in run_single_singulation.R


library(mirt)
library(irtimv)
source("00funs.R")

# ============================================================================
# FUNCTION 6: Run single simulation iteration
# ============================================================================
run_single_simulation <- function(th, a, b0, b1) {
  
  # Step 1: Generate initial data x
  x <- generate_grm_data(th, a, b=c(b0, b1))
  
  # Step 2: Fit models to x
  # Fit GRM
  grm_fit <- estimate_grm(x, th)
  
  # Dichotomize and fit 2PLs
  x_prime <- as.numeric(x > 0)
  x_double <- as.numeric(x > 1)
  
  fit_2pl_prime <- estimate_2pl(x_prime, th)
  fit_2pl_double <- estimate_2pl(x_double, th)
  
  # Step 3: Generate new test data x2
  x2 <- generate_grm_data(th, a, b=c(b0, b1))
  
  # Step 4: Get predictions for x2
  # Predictions from GRM (for polytomous response)
  pred_grm <- predict_grm(th, grm_fit$a, grm_fit$b)
  
  # Step 5: Calculate omega for dichotomizations
  # For x' dichotomization (x > 0)
  x2_prime <- as.numeric(x2 > 0)
  p_prime_2pl <- predict_2pl(th, fit_2pl_prime$a, fit_2pl_prime$b)
  
  z_prime <- data.frame(
    resp = x2_prime,
    p1 = mean(x_prime),
    p2 = p_prime_2pl
  )
  omega_x_prime <- imv(z_prime, p1 = "p1", p2 = "p2")
  
  # For x'' dichotomization (x > 1)
  x2_double <- as.numeric(x2 > 1)
  p_double_2pl <- predict_2pl(th, fit_2pl_double$a, fit_2pl_double$b)
  
  z_double <- data.frame(
    resp = x2_double,
    p1 = mean(x_double),
    p2 = p_double_2pl
  )
  omega_x_double <- imv(z_double, p1 = "p1", p2 = "p2")
  
  # Step 6: Calculate omega_c and omega_t
  # Reconstruct polytomous predictions from the two 2PLs
  
  # Create full prediction dataframe for imv_c and imv_t
    y <- data.frame(
        resp = x2
    )
    for (i in 1:ncol(pred_grm)) y[[paste('p1',i-1,sep='')]]<-mean(x==i-1)
    for (i in 1:ncol(pred_grm)) y[[paste('p2',i-1,sep='')]]<-pred_grm[,i]

    
  # Calculate frequency table
  pctt.tab <- table(factor(x2, levels = 0:2))
  pctt.tab <- as.numeric(pctt.tab)/length(x2)
  
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
for (b0 in seq(-1,.5,by=.5)) results[[as.character(b0)]] <- run_simulation(n_theta = 1000, n_sim = 50, seed = 123,b0=b0,mu=b0)
results2<-list()
for (b0 in seq(-1,.5,by=.5)) results2[[as.character(b0)]] <- run_simulation(n_theta = 1000, n_sim = 50, seed = 123,b0=b0,mu=0)


f<-function(results) {
    pf<-function(x,y,...) {
        m<-loess(y~x)
        lines(x,predict(m),...)
    }
    plot(NULL,xlim=c(-1,1),ylim=c(0,.5),xlab='second threshold',ylab='IMV')
                                        #pf(results$b1_true,results$omega_c)
    pf(results$b1_true,results$omega_t,lwd=2)
    pf(results$b1_true,results$omega_x_prime,col='red',lwd=2)
    pf(results$b1_true,results$omega_x_double,col='blue',lwd=2)
    pf(results$b1_true,results$omega_c,lty=2)
    abline(v=unique(results$b0))
    }
par(mfcol=c(4,2),mgp=c(2,1,0),mar=c(3,3,1,1),oma=rep(.5,4))
for (i in 1:length(results)) f(results[[i]])
legend("left",bty='n',fill=c("black","red","blue"),c("grm","2pl_12","2pl_2"),title='vert line: first threshold')
for (i in 1:length(results2)) f(results2[[i]])
legend("left",bty='n',fill=c("black","red","blue"),c("grm","2pl_12","2pl_2"),title='vert line: first threshold')



## a figure showing the collapse
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



z<-run_simulation(n_theta = 1000, n_sim = 100, seed = 123,b0=-1,mu=-1)

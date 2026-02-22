
# ============================================================================
# FUNCTION 1: Generate data from Graded Response Model
# ============================================================================
## generate_grm_data <- function(th, a, b0, b1) {
##   # th: vector of abilities
##   # a: discrimination parameter
##   # b0, b1: threshold parameters (b0 < b1)
  
##   n <- length(th)
  
##   # Calculate cumulative probabilities
##   # P(X >= 1 | th)
##   p1_star <- 1 / (1 + exp(-a * (th - b0)))
##   # P(X >= 2 | th)
##   p2_star <- 1 / (1 + exp(-a * (th - b1)))
  
##   # Calculate category probabilities
##   # P(X = 0)
##   p0 <- 1 - p1_star
##   # P(X = 1)
##   p1 <- p1_star - p2_star
##   # P(X = 2)
##   p2 <- p2_star
  
##   # Generate responses
##   x <- numeric(n)
##   for(i in 1:n) {
##     x[i] <- sample(0:2, size = 1, prob = c(p0[i], p1[i], p2[i]))
##   }
  
##   return(x)
## }

#' Generate data from the Graded Response Model (GRM)
#'
#' @param th A numeric vector of latent abilities (theta).
#' @param a A single numeric value for the discrimination parameter.
#' @param b A numeric vector of K-1 threshold (or boundary) parameters.
#'   The function will produce responses in K categories (0, 1, ..., K-1).
#'   The thresholds in `b` should be sorted in increasing order.
#'
#' @return A numeric vector of the same length as `th` containing the
#'   generated categorical responses (0 to K-1).
#'
generate_grm_data <- function(th, a, b) {
  # --- 1. Input Validation and Setup ---
  # Ensure b is sorted, which is a requirement for the GRM
  if (is.unsorted(b)) {
    warning("Threshold parameters 'b' should be sorted for the GRM. Sorting them now.")
    b <- sort(b)
  }
  
  n <- length(th)
  num_categories <- length(b) + 1
  
  # --- 2. Calculate Cumulative Probabilities ---
  # P*(k) = P(Response >= k | th), for k = 1, ..., K-1
  # This uses the 2-Parameter Logistic (2PL) model formula for each threshold.
  # The `outer` function efficiently calculates this for all combinations of th and b,
  # resulting in an n x (K-1) matrix.
  p_star <- outer(th, b, function(theta, b_k) {
    1 / (1 + exp(-a * (theta - b_k)))
  })
  
  # Ensure p_star is a matrix, even if n=1 or length(b)=1
  if (!is.matrix(p_star)) {
    p_star <- matrix(p_star, nrow = n)
  }
  
  # --- 3. Calculate Category Probabilities ---
  # P(Response = k) is calculated from the cumulative probabilities.
  # To do this efficiently, we augment the p_star matrix.
  # P(Response >= 0) is always 1.
  # P(Response >= K) is always 0.
  p_star_augmented <- cbind(1, p_star, 0)
  
  # Category prob P(k) = P*(k) - P*(k+1). This vectorized subtraction
  # calculates all category probabilities at once.
  # The result is an n x K matrix of probabilities.
  prob_matrix <- p_star_augmented[, 1:num_categories] - p_star_augmented[, 2:(num_categories + 1)]
  
  # Due to floating point precision, very small negative numbers can occur.
  # We set them to 0 and re-normalize each row to sum to 1.
  prob_matrix[prob_matrix < 0] <- 0
  prob_matrix <- prob_matrix / rowSums(prob_matrix)

  # --- 4. Generate Responses ---
  # For each person, sample a response from their specific probability distribution.
  responses <- numeric(n)
  categories <- 0:(num_categories - 1)
  for (i in 1:n) {
    responses[i] <- sample(categories, size = 1, prob = prob_matrix[i, ])
  }
  
  return(responses)
}

# ============================================================================
# FUNCTION 2: Estimate GRM parameters given theta
# ============================================================================
## estimate_grm <- function(x, th) {
##   # x: response vector (0, 1, 2)
##   # th: known ability vector
  
##   # Negative log-likelihood for GRM
##   nll_grm <- function(par) {
##     a <- exp(par[1])  # Use exp to ensure a > 0
##     b0 <- par[2]
##     b1 <- par[3]
    
##     if(b1 <= b0) return(1e10)  # Ensure b1 > b0
    
##     # Calculate cumulative probabilities
##     p1_star <- 1 / (1 + exp(-a * (th - b0)))
##     p2_star <- 1 / (1 + exp(-a * (th - b1)))
    
##     # Category probabilities
##     p0 <- 1 - p1_star
##     p1 <- p1_star - p2_star
##     p2 <- p2_star
    
##     # Ensure probabilities are valid
##     p0 <- pmax(p0, 1e-10)
##     p1 <- pmax(p1, 1e-10)
##     p2 <- pmax(p2, 1e-10)
    
##     # Log-likelihood
##     ll <- sum((x == 0) * log(p0) + (x == 1) * log(p1) + (x == 2) * log(p2))
    
##     return(-ll)
##   }
  
##   # Optimize
##   fit <- optim(par = c(0, -0.5, 0.5), fn = nll_grm, method = "BFGS")
  
##   a_hat <- exp(fit$par[1])
##   b0_hat <- fit$par[2]
##   b1_hat <- fit$par[3]
  
##   return(list(a = a_hat, b0 = b0_hat, b1 = b1_hat))
## }

#' Estimate Graded Response Model Parameters for Varying K
#'
#' This function estimates the discrimination (a) and threshold (b) parameters
#' for a single item under the Graded Response Model (GRM), given known
#' person abilities (theta). It generalizes the model for any number of 
#' response categories K (where responses are 0, 1, ..., K-1).
#'
#' @param x A numeric vector of observed responses for a single item. 
#'          Responses must be integers starting from 0 (e.g., 0, 1, 2, ...).
#' @param th A numeric vector of known person ability values (theta).
#'           Must be the same length as x.
#'
#' @return A list containing the estimated discrimination parameter 'a' and a
#'         named vector of threshold parameters 'b'.
#'
#' @examples
#' # 1. SETUP: Simulate data for a K=4 item (responses 0, 1, 2, 3)
#' set.seed(123)
#' n_persons <- 500
#' true_a <- 1.5
#' true_b <- c(-1.0, 0.5, 1.5) # K-1 = 3 thresholds
#' K <- length(true_b) + 1
#'
#' # Generate person abilities
#' theta <- rnorm(n_persons, 0, 1)
#'
#' # 2. SIMULATE RESPONSES based on the GRM
#' # Calculate cumulative probabilities (P*)
#' p_star <- matrix(NA, nrow = n_persons, ncol = K - 1)
#' for (k in 1:(K - 1)) {
#'   p_star[, k] <- 1 / (1 + exp(-true_a * (theta - true_b[k])))
#' }
#' 
#' # Calculate category probabilities (P)
#' p_cat <- matrix(NA, nrow = n_persons, ncol = K)
#' p_cat[, 1] <- 1 - p_star[, 1] # P(X=0)
#' for (k in 2:(K - 1)) {
#'   p_cat[, k] <- p_star[, k - 1] - p_star[, k] # P(X=k-1)
#' }
#' p_cat[, K] <- p_star[, K - 1] # P(X=K-1)
#'
#' # Generate observed responses by sampling
#' responses <- apply(p_cat, 1, function(p) sample(0:(K-1), 1, prob = p))
#'
#' # 3. ESTIMATE PARAMETERS using the generalized function
#' estimated_params <- estimate_grm_generalized(responses, theta)
#'
#' # 4. COMPARE true vs. estimated parameters
#' print("True Parameters:")
#' print(list(a = true_a, b = true_b))
#' 
#' print("Estimated Parameters:")
#' print(estimated_params)

estimate_grm <- function(x, th) {
  # Determine the number of categories K from the data
  # Assumes responses are 0, 1, ..., K-1
  K <- max(x, na.rm = TRUE) + 1
  
  if (K < 2) {
    stop("GRM requires at least 2 response categories.")
  }

  # Negative log-likelihood function for a generalized GRM
  nll_grm <- function(par) {
    # par[1] is log(a), par[2:K] are the b parameters
    a <- exp(par[1])
    b <- par[-1] # All parameters except the first

    # Constraint: b parameters must be ordered
    # diff(b) computes b[i] - b[i-1]. We need all these to be > 0.
    if (length(b) > 1 && any(diff(b) <= 0)) {
      return(1e10) # Return a large value to penalize invalid parameter sets
    }

    # --- Vectorized Calculation of Probabilities ---

    # 1. Calculate cumulative probabilities (P_star)
    # This matrix will have one row per person and one column per threshold (K-1)
    # outer() creates a matrix of all (th_i - b_j) combinations
    logits <- -a * outer(th, b, FUN = "-")
    p_star_matrix <- 1 / (1 + exp(logits))

    # 2. Calculate category probabilities (P)
    # We create a matrix with K columns for P(X=0), P(X=1), ..., P(X=K-1)
    cat_probs <- matrix(NA, nrow = length(th), ncol = K)

    # Prob for the first category (0)
    cat_probs[, 1] <- 1 - p_star_matrix[, 1]

    # Probs for middle categories (1 to K-2)
    if (K > 2) {
      for (k in 2:(K - 1)) {
        cat_probs[, k] <- p_star_matrix[, k - 1] - p_star_matrix[, k]
      }
    }
    
    # Prob for the last category (K-1)
    cat_probs[, K] <- p_star_matrix[, K - 1]

    # Ensure probabilities are not exactly zero to avoid log(0)
    cat_probs <- pmax(cat_probs, 1e-10)

    # 3. Calculate log-likelihood using matrix indexing
    # This is an efficient way to get the probability for each person's actual response
    # It avoids loops or complex sums.
    # We create a 2-column matrix to index `cat_probs`: 
    #   - 1st col: row number (person index)
    #   - 2nd col: column number (response category + 1, for 1-based indexing)
    idx <- cbind(seq_along(x), x + 1)
    
    # Extract the probabilities corresponding to the observed responses
    observed_probs <- cat_probs[idx]
    
    # Sum the log-probabilities and return the negative
    ll <- sum(log(observed_probs))
    
    return(-ll)
  }

  # Generate reasonable starting values for optimization
  # log(a)=0 => a=1. b's are spread out between -1 and 1.
  # Total parameters = 1 (for a) + (K-1) (for b's) = K
  if (K == 2) {
    # For a 2PL model, there's only one b
    start_vals <- c(0, 0) 
  } else {
    start_vals <- c(0, seq(-1, 1, length.out = K - 1))
  }
  
  # Optimize to find the parameters that minimize the negative log-likelihood
  fit <- optim(par = start_vals, fn = nll_grm, method = "BFGS")

  # Extract and name the final parameters
  a_hat <- exp(fit$par[1])
  b_hats <- fit$par[-1]
  names(b_hats) <- paste0("b", 0:(K - 2))

  return(list(a = a_hat, b = b_hats))
}

# ============================================================================
# FUNCTION 3: Estimate 2PL parameters given theta
# ============================================================================
estimate_2pl <- function(x_binary, th) {
  # x_binary: binary response vector (0, 1)
  # th: known ability vector
  
  # Negative log-likelihood for 2PL
  nll_2pl <- function(par) {
    a <- par[1]
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
  
  a_hat <- fit$par[1]
  b_hat <- fit$par[2]
  
  return(list(a = a_hat, b = b_hat))
}

# ============================================================================
# FUNCTION 4: Predict probabilities from GRM
# ============================================================================
## predict_grm <- function(th, a, b0, b1) {
##   # Returns matrix with columns: p0, p1, p2
  
##   p1_star <- 1 / (1 + exp(-a * (th - b0)))
##   p2_star <- 1 / (1 + exp(-a * (th - b1)))
  
##   p0 <- 1 - p1_star
##   p1 <- p1_star - p2_star
##   p2 <- p2_star
  
##   return(data.frame(p0 = p0, p1 = p1, p2 = p2))
## }

#' Predict Probabilities for the Graded Response Model
#'
#' Calculates the probability of responding in each category for a given set of
#' person abilities (theta), a discrimination parameter (a), and a vector of
#' threshold parameters (b).
#'
#' @param th A numeric vector of person ability values (theta).
#' @param a The discrimination parameter (a scalar, must be positive).
#' @param b A numeric vector of K-1 threshold parameters, which must be
#'          in strictly increasing order (b_1 < b_2 < ...).
#'
#' @return A data frame where each row corresponds to a theta value and each
#'         column corresponds to a response category probability (p0, p1, ...).
#'
#' @examples
#' # Define parameters for a K=4 item (responses 0, 1, 2, 3)
#' a_param <- 1.5
#' b_params <- c(-1.0, 0.5, 1.5) # K-1 = 3 thresholds
#'
#' # Define a range of theta values to predict probabilities for
#' theta_values <- seq(-3, 3, by = 0.5)
#'
#' # Get the predicted probabilities
#' predicted_probs <- predict_grm_generalized(theta_values, a_param, b_params)
#'
#' # Print the results
#' print(round(predicted_probs, 3))
#'
#' # --- Sanity Check ---
#' # The probabilities for each person (row) should sum to 1
#' print("Row sums should all be 1:")
#' print(rowSums(predicted_probs))
#'
#' # --- Plotting the results (optional but very useful) ---
#' # library(ggplot2)
#' # library(tidyr)
#' #
#' # predicted_probs$theta <- theta_values
#' # predicted_probs %>%
#' #   pivot_longer(-theta, names_to = "Category", values_to = "Probability") %>%
#' #   ggplot(aes(x = theta, y = Probability, color = Category)) +
#' #   geom_line(linewidth = 1.2) +
#' #   labs(title = "Category Response Curves (CRC) for a GRM Item",
#' #        x = "Ability (Theta)",
#' #        y = "Probability of Response") +
#' #   theme_minimal()

predict_grm <- function(th, a, b) {
  # --- Input Validation ---
  if (a <= 0) {
    stop("Discrimination parameter 'a' must be positive.")
  }
  if (length(b) > 1 && any(diff(b) <= 0)) {
    stop("Threshold parameters 'b' must be in strictly increasing order.")
  }

  # Number of categories K is determined by the number of thresholds
  K <- length(b) + 1
  
  # --- Vectorized Calculation ---

  # 1. Calculate cumulative probabilities (P*)
  # Creates a matrix where element [i, j] is P*(theta_i, b_j)
  logits <- -a * outer(th, b, FUN = "-")
  p_star_matrix <- 1 / (1 + exp(logits))

  # 2. Calculate category probabilities (P)
  # Initialize a matrix to hold the results
  cat_probs <- matrix(NA, nrow = length(th), ncol = K)

  # Probability for the first category (0)
  # P(X=0) = 1 - P*(b_1)
  cat_probs[, 1] <- 1 - p_star_matrix[, 1]

  # Probabilities for middle categories (if they exist)
  # P(X=k) = P*(b_k) - P*(b_{k+1})
  if (K > 2) {
    for (k in 2:(K - 1)) {
      # The k-th column of cat_probs is for response k-1.
      # It uses the (k-1)th and k-th columns of p_star_matrix.
      cat_probs[, k] <- p_star_matrix[, k - 1] - p_star_matrix[, k]
    }
  }

  # Probability for the last category (K-1)
  # P(X=K-1) = P*(b_{K-1})
  cat_probs[, K] <- p_star_matrix[, K - 1]
  
  # Name the columns for clarity
  colnames(cat_probs) <- paste0("p", 0:(K - 1))
  
  return(as.data.frame(cat_probs))
}

# ============================================================================
# FUNCTION : Predict probabilities from 2PL
# ============================================================================
predict_2pl <- function(th, a, b) {
  # Returns probability of response = 1
  p <- 1 / (1 + exp(-a * (th - b)))
  return(p)
}


# ============================================================================
# FUNCTION : Predict data from PCM
# ============================================================================

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

#' Generate data from the Partial Credit Model (PCM)
#'
#' @param th A numeric vector of latent abilities (theta).
#' @param a A single numeric value for the discrimination parameter (often fixed to 1
#'   in the original PCM, but included here for the Generalized PCM).
#' @param d A numeric vector of K-1 step difficulty parameters (d_1, d_2, ...).
#'   The function will produce responses in K categories (0, 1, ..., K-1).
#'
#' @return A numeric vector of the same length as `th` containing the
#'   generated categorical responses (0 to K-1).
#'
#' @details The probability of a response in category k is given by:
#'   P(X=k) = exp(a * sum_{j=1 to k}(th - d_j)) / (sum_{m=0 to K-1} exp(a * sum_{j=1 to m}(th - d_j)))
#'   where the sum from j=1 to 0 is defined as 0.
#'
generate_pcm_data <- function(th, a, d) {
  # --- 1. Setup ---
  n <- length(th)
  num_categories <- length(d) + 1
  categories <- 0:(num_categories - 1)
  
  # --- 2. Calculate Numerators of the Probabilities ---
  # The exponent for category k > 0 is: a * (k*th - sum(d_1, ..., d_k))
  # We will calculate these for all persons and categories in a vectorized way.
  
  # Create a matrix to hold the unnormalized log-probabilities (the exponents)
  # It has n rows (persons) and K columns (categories)
  log_p_numerators <- matrix(0, nrow = n, ncol = num_categories)
  
  # The cumulative sum of d parameters is needed for the formula.
  # d_cumsum will be [d_1, d_1+d_2, d_1+d_2+d_3, ...]
  d_cumsum <- cumsum(d)
  
  # For categories k = 1, 2, ..., K-1
  for (k in 1:(num_categories - 1)) {
    # The exponent is a * (k*theta - sum of d's up to k)
    # The column index is k+1 because R is 1-based and our categories are 0-based
    log_p_numerators[, k + 1] <- a * (k * th - d_cumsum[k])
  }
  
  # The exponents are now calculated. Exponentiate to get the numerators.
  p_numerators <- exp(log_p_numerators)
  
  # --- 3. Normalize to get Final Probabilities ---
  
  # The denominator for each person is the sum of their numerators across categories
  denominator <- rowSums(p_numerators)
  
  # Divide each person's numerators by their denominator to get the probability matrix.
  # R's recycling handles this row-wise division correctly.
  prob_matrix <- p_numerators / denominator
  
  # --- 4. Generate Responses ---
  # For each person, sample one response from their specific probability distribution.
  responses <- numeric(n)
  for (i in 1:n) {
    responses[i] <- sample(categories, size = 1, prob = prob_matrix[i, ])
  }
  
  return(responses)
}

# ============================================================================
# FUNCTION : Estimate the PCM
# ============================================================================
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

#' Estimate (Generalized) Partial Credit Model Parameters
#'
#' Estimates the discrimination (a) and step difficulty (d) parameters
#' for a single item under the (G)PCM, given known person abilities (theta).
#' The function generalizes for any number of response categories K.
#'
#' @param x A numeric vector of observed responses for a single item. 
#'          Responses must be integers starting from 0 (e.g., 0, 1, 2, ...).
#' @param th A numeric vector of known person ability values (theta).
#'           Must be the same length as x.
#'
#' @return A list containing the estimated discrimination parameter 'a' and a
#'         named vector of step difficulty parameters 'd'. Note: The first step
#'         difficulty, d_0, is implicitly fixed to 0 for model identification.
#'         The function estimates d_1, d_2, ..., d_{K-1}.
#'
#' @examples
#' # 1. SETUP: Simulate data for a K=4 item (responses 0, 1, 2, 3)
#' set.seed(456)
#' n_persons <- 1000
#' true_a <- 1.2
#' true_d <- c(-1.0, 0.2, 0.8) # K-1 = 3 step difficulties
#' 
#' # Generate person abilities
#' theta <- rnorm(n_persons, 0, 1)
#'
#' # 2. SIMULATE RESPONSES using the generalized data generation function
#' responses <- generate_pcm_data(th = theta, a = true_a, d = true_d)
#'
#' # 3. ESTIMATE PARAMETERS using the generalized estimation function
#' estimated_params <- estimate_pcm_generalized(responses, theta)
#'
#' # 4. COMPARE true vs. estimated parameters
#' print("True Parameters:")
#' print(list(a = true_a, d = true_d))
#' 
#' print("Estimated Parameters:")
#' print(estimated_params)

estimate_pcm <- function(x, th) {
  # Determine the number of categories K from the data
  K <- max(x, na.rm = TRUE) + 1
  n <- length(th)
  
  if (K < 2) {
    stop("PCM requires at least 2 response categories.")
  }

  # Negative log-likelihood function for a generalized PCM
  nll_pcm <- function(par) {
    # par[1] is a, par[2:K] are the d_1, ..., d_{K-1} parameters
    a <- par[1]
    d <- par[-1]

    # --- Vectorized Calculation of Probabilities ---
    
    # 1. Calculate the exponents for the numerators
    # The exponent for category k is: a * (k*th - sum(d_1, ..., d_k))
    log_p_numerators <- matrix(0, nrow = n, ncol = K)
    d_cumsum <- cumsum(d)
    
    for (k in 1:(K - 1)) {
      # Column k+1 corresponds to response category k
      log_p_numerators[, k + 1] <- a * (k * th - d_cumsum[k])
    }

    # 2. Calculate the full numerators and the denominator
    p_numerators <- exp(log_p_numerators)
    denominator <- rowSums(p_numerators)
    
    # 3. Calculate the final category probabilities
    # Use pmax to avoid division by zero if all numerators are zero (highly unlikely)
    prob_matrix <- p_numerators / pmax(denominator, 1e-10)

    # 4. Calculate log-likelihood using efficient matrix indexing
    # This extracts the probability corresponding to each person's observed response
    idx <- cbind(seq_along(x), x + 1)
    observed_probs <- prob_matrix[idx]
    
    # Use pmax to avoid log(0) for stability
    ll <- sum(log(pmax(observed_probs, 1e-10)))
    
    return(-ll)
  }

  # --- Optimization ---
  
  # Generate sensible starting values for optimization
  # log(a)=0 => a=1. d's are spread out.
  # Total parameters = 1 (for a) + (K-1) (for d's) = K
  if (K == 2) {
    # For a 2PL model equivalent, there's only one d
    start_vals <- c(0, 0) 
  } else {
    # Spread initial d values around 0
    start_vals <- c(0, seq(-0.5, 0.5, length.out = K - 1))
  }

  # Optimize to find parameters that minimize the negative log-likelihood
  fit <- optim(par = start_vals, fn = nll_pcm, method = "BFGS")

  # --- Extract and Return Final Parameters ---
  a_hat <- fit$par[1]
  d_hats <- fit$par[-1]
  names(d_hats) <- paste0("d", 1:(K - 1))

  return(list(a = a_hat, d = d_hats))
}

# ============================================================================
# FUNCTION: Predict probabilities from PCM
# ============================================================================
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
#' Predict Probabilities for the (Generalized) Partial Credit Model
#'
#' Calculates the probability of responding in each category for a given set of
#' person abilities (theta), a discrimination parameter (a), and a vector of
#' step difficulty parameters (d).
#'
#' @param th A numeric vector of person ability values (theta).
#' @param a The discrimination parameter (a scalar).
#' @param d A numeric vector of K-1 step difficulty parameters (d_1, d_2, ...).
#'
#' @return A data frame where each row corresponds to a theta value and each
#'         column corresponds to a response category probability (p0, p1, ...).
#'
#' @examples
#' # Define parameters for a K=4 item (responses 0, 1, 2, 3)
#' a_param <- 1.2
#' d_params <- c(-1.0, 0.2, 0.8) # K-1 = 3 step difficulties
#'
#' # Define a range of theta values to predict probabilities for
#' theta_values <- seq(-4, 4, by = 0.5)
#'
#' # Get the predicted probabilities
#' predicted_probs <- predict_pcm_generalized(theta_values, a_param, d_params)
#'
#' # Print the results
#' print(round(predicted_probs, 3))
#'
#' # --- Sanity Check ---
#' # The probabilities for each person (row) should sum to 1
#' cat("\nRow sums should all be 1:\n")
#' print(rowSums(predicted_probs))
#'
#' # --- Plotting the results (optional but very useful) ---
#' # library(ggplot2)
#' # library(tidyr)
#' #
#' # predicted_probs$theta <- theta_values
#' # predicted_probs %>%
#' #   pivot_longer(-theta, names_to = "Category", values_to = "Probability") %>%
#' #   ggplot(aes(x = theta, y = Probability, color = Category)) +
#' #   geom_line(linewidth = 1.2) +
#' #   labs(title = "Category Response Curves (CRC) for a PCM Item",
#' #        x = "Ability (Theta)",
#' #        y = "Probability of Response") +
#' #   theme_minimal()

predict_pcm <- function(th, a, d) {
  # --- 1. Setup ---
  n <- length(th)
  num_categories <- length(d) + 1
  
  # --- 2. Calculate Numerators of the Probabilities ---
  # The exponent for category k > 0 is: a * (k*th - sum(d_1, ..., d_k))
  
  # Create a matrix to hold the unnormalized log-probabilities (the exponents)
  log_p_numerators <- matrix(0, nrow = n, ncol = num_categories)
  
  # Pre-calculate the cumulative sum of d parameters for efficiency
  d_cumsum <- cumsum(d)
  
  # Loop over categories k = 1, 2, ..., K-1
  for (k in 1:(num_categories - 1)) {
    # The exponent is a * (k*theta - sum of d's up to k)
    # The column index is k+1 (for response k)
    log_p_numerators[, k + 1] <- a * (k * th - d_cumsum[k])
  }
  
  # Exponentiate to get the actual numerators
  p_numerators <- exp(log_p_numerators)
  
  # --- 3. Normalize to get Final Probabilities ---
  
  # The denominator for each person is the sum of their numerators
  denominator <- rowSums(p_numerators)
  
  # Divide each numerator by the corresponding denominator
  # R's recycling handles this row-wise division correctly
  prob_matrix <- p_numerators / denominator
  
  # --- 4. Format and Return Output ---
  
  # Name the columns for clarity (p0, p1, p2, ...)
  colnames(prob_matrix) <- paste0("p", 0:(num_categories - 1))
  
  return(as.data.frame(prob_matrix))
}

# ============================================================================
# FUNCTION: Predict probabilities from PCM
# ============================================================================


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

# ============================================================================
# FUNCTION: Predict probabilities from PCM
# ============================================================================

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

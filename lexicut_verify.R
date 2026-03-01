################################################################################
# LexiCut Conjecture Verification for n > 20
# Efficient R code (no external packages required)
#
# Conjecture: For n >= 7 and r in (1,2), the optimal cut is always a k-isolated cut.
#
# Strategy:
#   1. Compute thresholds r_k(n) efficiently
#   2. Exhaustively verify near-isolated cuts S*_k = {1,...,k,n} (only known competitors)
#   3. Random sampling for other cuts
################################################################################

#-------------------------------------------------------------------------------
# SECTION 1: Basic Setup and Helper Functions
#-------------------------------------------------------------------------------

# Number of edges in K_n
N_edges <- function(n) choose(n, 2)

# Edge index (i,j) in lexicographic order (1-indexed)
# Edges: (1,2), (1,3), ..., (1,n), (2,3), ..., (n-1,n)
edge_index <- function(i, j, n) {
  if (i >= j) stop("Require i < j")
  # Edges from vertex i: (i, i+1), ..., (i, n) have indices starting at
  # 1 + sum_{m=1}^{i-1}(n-m) = 1 + (i-1)*n - i*(i-1)/2
  start_i <- (i - 1) * n - i * (i - 1) / 2
  return(start_i + (j - i))
}

# Vectorized edge index for pairs
edge_index_vec <- function(i_vec, j_vec, n) {
  start_i <- (i_vec - 1) * n - i_vec * (i_vec - 1) / 2
  return(start_i + (j_vec - i_vec))
}

#-------------------------------------------------------------------------------
# SECTION 2: Precompute Edge Exponents (Key Optimization)
#-------------------------------------------------------------------------------

# Precompute exponents for k-isolated cut
precompute_isolated_exponents <- function(k, n) {
  Nn <- N_edges(n)
  
  # Crossing edges: (i, j) where i <= k < j
  pairs <- expand.grid(i = 1:k, j = (k+1):n)
  indices <- edge_index_vec(pairs$i, pairs$j, n)
  exponents <- Nn - indices
  
  return(sort(exponents, decreasing = TRUE))
}

# Precompute exponents for near-isolated cut S*_k = {1,...,k,n}
precompute_near_isolated_exponents <- function(k, n) {
  Nn <- N_edges(n)
  
  S_star <- c(1:k, n)
  complement <- (k+1):(n-1)
  
  # Crossing edges between S_star and complement
  pairs <- expand.grid(v1 = S_star, v2 = complement)
  pairs$i <- pmin(pairs$v1, pairs$v2)
  pairs$j <- pmax(pairs$v1, pairs$v2)
  
  indices <- edge_index_vec(pairs$i, pairs$j, n)
  exponents <- Nn - indices
  
  return(sort(exponents, decreasing = TRUE))
}

#-------------------------------------------------------------------------------
# SECTION 3: Efficient Weight Computation Using Log-Sum-Exp
#-------------------------------------------------------------------------------

# Compute sum(r^exponents) stably using log-sum-exp trick for large exponents
# Returns log of the weight when use_log = TRUE
weight_from_exponents <- function(exponents, r, use_log = FALSE) {
  if (length(exponents) == 0) return(if(use_log) -Inf else 0)
  
  log_r <- log(r)
  log_terms <- exponents * log_r
  
  # Log-sum-exp: log(sum(exp(x))) = max(x) + log(sum(exp(x - max(x))))
  max_log <- max(log_terms)
  log_sum <- max_log + log(sum(exp(log_terms - max_log)))
  
  if (use_log) {
    return(log_sum)
  } else {
    return(exp(log_sum))
  }
}

# Compare two weights given as exponent vectors
# Returns: positive if first > second, negative if first < second
# Uses log-space for numerical stability
compare_weights_log <- function(exp1, exp2, r) {
  log_w1 <- weight_from_exponents(exp1, r, use_log = TRUE)
  log_w2 <- weight_from_exponents(exp2, r, use_log = TRUE)
  return(log_w1 - log_w2)
}

#-------------------------------------------------------------------------------
# SECTION 4: Threshold Computation
#-------------------------------------------------------------------------------

# Threshold polynomial P^{n,k}(r) = W(C_k) - W(C_{k+1})
threshold_poly <- function(r, exp_k, exp_k1) {
  weight_from_exponents(exp_k, r) - weight_from_exponents(exp_k1, r)
}

# Compute threshold r_k(n) using uniroot with log-space comparison
compute_threshold <- function(k, n, tol = 1e-12) {
  if (k >= floor(n / 2)) return(NA)
  
  exp_k <- precompute_isolated_exponents(k, n)
  exp_k1 <- precompute_isolated_exponents(k + 1, n)
  
  # Use log-space for numerical stability with large exponents
  f <- function(r) {
    log_r <- log(r)
    
    # Compute log-sum-exp for each cut
    log_terms_k <- exp_k * log_r
    log_terms_k1 <- exp_k1 * log_r
    
    max_k <- max(log_terms_k)
    max_k1 <- max(log_terms_k1)
    
    log_w_k <- max_k + log(sum(exp(log_terms_k - max_k)))
    log_w_k1 <- max_k1 + log(sum(exp(log_terms_k1 - max_k1)))
    
    # Return difference of logs (sign preserved)
    return(log_w_k - log_w_k1)
  }
  
  # Check endpoints
  f1 <- tryCatch(f(1.0001), error = function(e) NA)
  f2 <- tryCatch(f(1.9999), error = function(e) NA)
  
  if (is.na(f1) || is.na(f2) || !is.finite(f1) || !is.finite(f2)) {
    return(NA)
  }
  
  if (sign(f1) == sign(f2)) {
    # No root in interval
    return(NA)
  }
  
  result <- tryCatch({
    uniroot(f, interval = c(1.0001, 1.9999), tol = tol)$root
  }, error = function(e) NA)
  
  return(result)
}

# Compute all thresholds for given n
compute_all_thresholds <- function(n, tol = 1e-12) {
  max_k <- floor(n / 2) - 1
  if (max_k < 1) return(numeric(0))
  
  thresholds <- sapply(1:max_k, function(k) compute_threshold(k, n, tol))
  names(thresholds) <- paste0("r_", 1:max_k)
  return(thresholds)
}

#-------------------------------------------------------------------------------
# SECTION 5: Find Optimal k for Given r
#-------------------------------------------------------------------------------

find_optimal_k <- function(r, n, thresholds = NULL) {
  if (is.null(thresholds)) {
    thresholds <- compute_all_thresholds(n)
  }
  
  max_k <- floor(n / 2)
  
  # Filter out NA thresholds
  valid_thresh <- thresholds[!is.na(thresholds)]
  
  if (length(valid_thresh) == 0) {
    # No valid thresholds, use balanced partition
    return(max_k)
  }
  
  if (r >= valid_thresh[1]) return(1)
  
  # Find the right interval
  for (i in 1:(length(valid_thresh) - 1)) {
    if (r > valid_thresh[i + 1] && r <= valid_thresh[i]) {
      return(i + 1)
    }
  }
  
  # r <= smallest threshold
  return(max_k)
}

#-------------------------------------------------------------------------------
# SECTION 6: Main Verification Function
#-------------------------------------------------------------------------------

verify_conjecture <- function(n, n_r_values = 100, n_samples = 10000,
                              seed = 42, verbose = TRUE) {
  
  if (verbose) {
    cat(paste(rep("=", 60), collapse = ""), "\n")
    cat(sprintf("LexiCut Conjecture Verification for n = %d\n", n))
    cat(paste(rep("=", 60), collapse = ""), "\n")
  }
  
  Nn <- N_edges(n)
  if (verbose) cat(sprintf("Number of edges N = %d\n", Nn))
  
  # Generate r values - focus more on the interesting region near 1
  r_values <- c(
    seq(1.001, 1.1, length.out = ceiling(n_r_values / 2)),
    seq(1.1, 1.999, length.out = floor(n_r_values / 2))
  )
  r_values <- unique(sort(r_values))
  
  # Step 1: Compute thresholds
  if (verbose) cat("\n1. Computing thresholds...\n")
  t_start <- Sys.time()
  thresholds <- compute_all_thresholds(n)
  t_thresh <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  
  if (verbose) {
    cat(sprintf("   Time: %.2f seconds\n", t_thresh))
    cat("   Thresholds:\n")
    for (k in seq_along(thresholds)) {
      if (!is.na(thresholds[k])) {
        cat(sprintf("     r_%d(%d) = %.8f\n", k, n, thresholds[k]))
      }
    }
  }
  
  # Step 2: Precompute all exponents
  if (verbose) cat("\n2. Precomputing edge exponents...\n")
  t_start <- Sys.time()
  
  max_k_iso <- floor(n / 2)
  max_k_near <- n - 2
  
  iso_exponents <- lapply(1:max_k_iso, function(k) precompute_isolated_exponents(k, n))
  near_exponents <- lapply(1:max_k_near, function(k) precompute_near_isolated_exponents(k, n))
  
  t_precomp <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (verbose) cat(sprintf("   Time: %.2f seconds\n", t_precomp))
  
  # Step 3: Verify near-isolated cuts (exhaustive)
  if (verbose) cat("\n3. Verifying near-isolated cuts (exhaustive)...\n")
  t_start <- Sys.time()
  
  near_violations <- data.frame(
    r = numeric(), k_near = integer(), k_opt = integer(),
    log_diff = numeric()
  )
  
  for (r in r_values) {
    k_opt <- find_optimal_k(r, n, thresholds)
    log_W_opt <- weight_from_exponents(iso_exponents[[k_opt]], r, use_log = TRUE)
    
    for (k in 1:max_k_near) {
      log_W_near <- weight_from_exponents(near_exponents[[k]], r, use_log = TRUE)
      
      # Compare in log space
      log_diff <- log_W_near - log_W_opt
      if (log_diff > 1e-12) {  # near > opt
        near_violations <- rbind(near_violations, data.frame(
          r = r, k_near = k, k_opt = k_opt, log_diff = log_diff
        ))
      }
    }
  }
  
  t_near <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (verbose) {
    cat(sprintf("   Time: %.2f seconds\n", t_near))
    cat(sprintf("   Violations found: %d\n", nrow(near_violations)))
    if (nrow(near_violations) > 0) {
      cat("   First few violations:\n")
      print(head(near_violations, 5))
    }
  }
  
  # Step 4: Random sampling
  if (verbose) cat(sprintf("\n4. Random sampling (%d samples)...\n", n_samples))
  t_start <- Sys.time()
  set.seed(seed)
  
  random_violations <- data.frame(
    r = numeric(), S = character(), k_opt = integer(), log_diff = numeric()
  )
  
  # Sample random cuts
  progress_interval <- max(1, n_samples %/% 10)
  
  for (i in 1:n_samples) {
    # Random membership for vertices 2 to n
    membership <- c(1, sample(0:1, n - 1, replace = TRUE))
    S <- which(membership == 1)
    
    # Skip trivial cuts
    if (length(S) == 1 || length(S) == n) next
    
    # Skip isolated cuts (already verified optimal among themselves)
    if (all(diff(S) == 1) && S[1] == 1) next
    
    # Compute exponents for this random cut
    complement <- setdiff(1:n, S)
    pairs <- expand.grid(v1 = S, v2 = complement)
    pairs$i <- pmin(pairs$v1, pairs$v2)
    pairs$j <- pmax(pairs$v1, pairs$v2)
    random_exponents <- Nn - edge_index_vec(pairs$i, pairs$j, n)
    
    for (r in r_values) {
      k_opt <- find_optimal_k(r, n, thresholds)
      log_W_opt <- weight_from_exponents(iso_exponents[[k_opt]], r, use_log = TRUE)
      log_W_random <- weight_from_exponents(random_exponents, r, use_log = TRUE)
      
      log_diff <- log_W_random - log_W_opt
      if (log_diff > 1e-12) {
        random_violations <- rbind(random_violations, data.frame(
          r = r, S = paste(S, collapse = ","), k_opt = k_opt, log_diff = log_diff
        ))
      }
    }
    
    if (verbose && i %% progress_interval == 0) {
      cat(sprintf("   Progress: %d/%d samples (%.0f%%)\n", 
                  i, n_samples, 100 * i / n_samples))
    }
  }
  
  t_random <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (verbose) {
    cat(sprintf("   Time: %.2f seconds\n", t_random))
    cat(sprintf("   Violations found: %d\n", nrow(random_violations)))
    if (nrow(random_violations) > 0) {
      cat("   First few violations:\n")
      print(head(random_violations, 5))
    }
  }
  
  # Summary
  total_violations <- nrow(near_violations) + nrow(random_violations)
  total_time <- t_thresh + t_precomp + t_near + t_random
  
  if (verbose) {
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("SUMMARY\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
    cat(sprintf("Total time: %.2f seconds\n", total_time))
    cat(sprintf("Near-isolated violations: %d (exhaustive over %d r-values)\n", 
                nrow(near_violations), length(r_values)))
    cat(sprintf("Random sampling violations: %d (%d samples x %d r-values)\n",
                nrow(random_violations), n_samples, length(r_values)))
    cat("\n")
    
    if (total_violations == 0) {
      cat(sprintf("*** CONJECTURE VERIFIED for n = %d ***\n", n))
    } else {
      cat(sprintf("*** CONJECTURE FAILED for n = %d ***\n", n))
    }
    cat(paste(rep("=", 60), collapse = ""), "\n")
  }
  
  return(list(
    n = n,
    thresholds = thresholds,
    near_violations = near_violations,
    random_violations = random_violations,
    verified = (total_violations == 0),
    timing = list(
      thresholds = t_thresh,
      precompute = t_precomp,
      near_isolated = t_near,
      random = t_random,
      total = total_time
    )
  ))
}

#-------------------------------------------------------------------------------
# SECTION 7: Batch Verification
#-------------------------------------------------------------------------------

verify_batch <- function(n_values, n_r_values = 100, n_samples = 10000,
                         seed = 42, verbose = TRUE) {
  results <- list()
  
  cat("\n", paste(rep("#", 70), collapse = ""), "\n")
  cat("# BATCH VERIFICATION OF LEXICUT CONJECTURE\n")
  cat("# n values:", paste(n_values, collapse = ", "), "\n")
  cat(paste(rep("#", 70), collapse = ""), "\n\n")
  
  for (n in n_values) {
    results[[as.character(n)]] <- verify_conjecture(
      n, n_r_values = n_r_values, n_samples = n_samples,
      seed = seed, verbose = verbose
    )
    cat("\n")
  }
  
  # Summary table
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("BATCH SUMMARY\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("%-6s %-12s %-12s %-12s %-10s\n", 
              "n", "Near Viol.", "Random Viol.", "Time (s)", "Status"))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  for (n in n_values) {
    result <- results[[as.character(n)]]
    status <- if(result$verified) "VERIFIED" else "FAILED"
    cat(sprintf("%-6d %-12d %-12d %-12.2f %-10s\n", 
                n, nrow(result$near_violations), nrow(result$random_violations),
                result$timing$total, status))
  }
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  return(results)
}

#-------------------------------------------------------------------------------
# SECTION 8: Focused Threshold Analysis
#-------------------------------------------------------------------------------

# Verify near threshold boundaries with fine granularity
verify_near_thresholds <- function(n, n_points_per_threshold = 50, verbose = TRUE) {
  if (verbose) cat(sprintf("Focused threshold analysis for n = %d\n", n))
  
  thresholds <- compute_all_thresholds(n)
  max_k_near <- n - 2
  
  # Precompute exponents
  max_k_iso <- floor(n / 2)
  iso_exponents <- lapply(1:max_k_iso, function(k) precompute_isolated_exponents(k, n))
  near_exponents <- lapply(1:max_k_near, function(k) precompute_near_isolated_exponents(k, n))
  
  violations <- data.frame()
  
  for (k in seq_along(thresholds)) {
    thresh <- thresholds[k]
    if (is.na(thresh)) next
    
    # Test points around threshold
    delta <- 0.01
    r_test <- seq(thresh - delta, thresh + delta, length.out = n_points_per_threshold)
    r_test <- r_test[r_test > 1 & r_test < 2]
    
    for (r in r_test) {
      k_opt <- find_optimal_k(r, n, thresholds)
      W_opt <- weight_from_exponents(iso_exponents[[k_opt]], r)
      
      for (k_near in 1:max_k_near) {
        W_near <- weight_from_exponents(near_exponents[[k_near]], r)
        rel_diff <- (W_near - W_opt) / max(abs(W_opt), 1e-10)
        
        if (rel_diff > 1e-9) {
          violations <- rbind(violations, data.frame(
            near_threshold = thresh, r = r, k_near = k_near, k_opt = k_opt,
            diff = W_near - W_opt
          ))
        }
      }
    }
    
    if (verbose) {
      cat(sprintf("  Checked around r_%d = %.6f\n", k, thresh))
    }
  }
  
  if (verbose) {
    if (nrow(violations) == 0) {
      cat("  No violations found near thresholds\n")
    } else {
      cat(sprintf("  Found %d violations near thresholds!\n", nrow(violations)))
    }
  }
  
  return(violations)
}

#-------------------------------------------------------------------------------
# SECTION 9: Display Information
#-------------------------------------------------------------------------------

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("LexiCut Conjecture Verification Code\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Main Conjecture: For n >= 7 and r in (1,2), the optimal cut\n")
cat("is always a k-isolated cut C_k = {1,...,k} | {k+1,...,n}.\n\n")

cat("Available functions:\n")
cat("  verify_conjecture(n, n_r_values, n_samples)\n")
cat("    - Main verification for single n\n\n")
cat("  verify_batch(n_values, n_r_values, n_samples)\n")
cat("    - Batch verification for multiple n values\n\n")
cat("  compute_all_thresholds(n)\n")
cat("    - Compute threshold values r_k(n)\n\n")
cat("  verify_near_thresholds(n)\n")
cat("    - Focused analysis near threshold boundaries\n\n")

cat("Example usage:\n")
cat("  result <- verify_conjecture(25)  # Verify for n=25\n")
cat("  results <- verify_batch(21:30)   # Batch verify n=21 to 30\n\n")
# Matrizen

S_vals <- c(-4,-3,-2,-1,-1/2,-1/3,-1/4, 1/4,1/3,1/2,1,2,3,4)

# Hilfsfunktion A mit gewünschten Kriterien
ok_matrix <- function(A, region = c("inside","outside"),
                       det_range = NULL, smax_max = NULL) {
  region <- match.arg(region)
  ev <- eigen(A)$values
  mags <- Mod(ev)
  if (region == "inside"  && !all(mags < 1)) return(FALSE)
  if (region == "outside" && !all(mags > 1)) return(FALSE)
  if (!is.null(det_range)) {
    d <- det(A)
    if (d < det_range[1] || d > det_range[2]) return(FALSE) ########### kleiner als -2???
  }
  if (!is.null(smax_max)) {
    s <- svd(A, nu = 0, nv = 0)$d
    if (max(s) > smax_max) return(FALSE)
  }
  TRUE
}

# zufällige 2x2-Matrix aus S, mit Spektrum + optionalen Nebenbedingungen
rand_2x2_from_S <- function(region = c("inside","outside"),
                            S = S_vals, max_tries = 100000,
                            det_range = NULL, smax_max = NULL) {
  region <- match.arg(region)
  for (k in seq_len(max_tries)) {
    A <- matrix(sample(S, 4, replace = TRUE), 2, 2)
    if (ok_matrix(A, region, det_range, smax_max)) return(A)
  }
  stop("Kein passendes A gefunden (max_tries erreicht).")
}

# Pools bauen (nebstdem: dieselben Nebenbedingungen weiterreichen)
make_matrix_pool <- function(region = c("inside","outside"),
                             pool_size = 200L, S = S_vals,
                             det_range = NULL, smax_max = NULL) {
  region <- match.arg(region)
  pool <- vector("list", pool_size)
  for (i in seq_len(pool_size)) {
    pool[[i]] <- rand_2x2_from_S(region, S, max_tries = 100000,
                                 det_range = det_range, smax_max = smax_max)
  }
  pool
}

inside_pool  <- make_matrix_pool("inside",  pool_size = 300,
                                 det_range = c(-2, 2), smax_max = 4)
outside_pool <- make_matrix_pool("outside", pool_size = 300,
                                 det_range = c(-2, 2), smax_max = 4)

### ohne null

simulate_2D_both_rect_mat <- function(p, inside_pool, outside_pool,
                                      x0, box_limit = 4, n_trials, tol = 1e-12) {
  success_count <- 0L
  total_steps   <- 0L
  
  for (i in seq_len(n_trials)) {
    x <- x0
    steps <- 0L
    
    while (!(all(abs(x) < tol)) && all(abs(x) <= box_limit)) {
      # Mit Wkt. p eine Outside-Matrix, sonst Inside-Matrix
      if (runif(1) < p) {
        A <- outside_pool[[ sample.int(length(outside_pool), 1) ]]
      } else {
        A <- inside_pool[[ sample.int(length(inside_pool), 1) ]]
      }
      
      # Matrix-Update
      x <- as.vector(A %*% x)
      steps <- steps + 1L
    }
    
    # Erfolg = rausgeflogen
    if (any(abs(x) > box_limit)) {
      success_count <- success_count + 1L
    }
    
    total_steps <- total_steps + steps
  }
  
  list(
    success_prob      = success_count / n_trials,
    expected_duration = as.numeric(total_steps) / n_trials
  )
}

# Parameter (wie bei dir)
x0 <- c(1, 1)
box_limit <- 4
p_grid <- seq(0.01, 0.99, by = 0.025)
n_trials <- 100000  # ggf. anheben

success_probs <- numeric(length(p_grid))
durations     <- numeric(length(p_grid))

system.time({
  for (i in seq_along(p_grid)) {
    p <- p_grid[i]
    
    res <- simulate_2D_both_rect_mat(
      p          = p,
      inside_pool = inside_pool,
      outside_pool = outside_pool,
      x0         = x0,
      box_limit  = box_limit,
      n_trials   = n_trials
    )
    
    success_probs[i] <- res$success_prob
    durations[i]     <- res$expected_duration
    
    print(p)
  }
})

par(mar = c(5, 5, 4, 2) + 0.1) 
plot(p_grid, success_probs, type = "l", lwd = 2,
     xlab = expression(p), ylab = expression(P(S[T] == b)),
     col = "deepskyblue3",bty = "l", cex.axis = 1.5, cex.lab  = 1.5)

plot(p_grid, durations, type = "l", lwd = 2,
     xlab = expression(p), ylab = expression(E(T)),
     col = "deepskyblue3", bty = "l", cex.axis = 1.5, cex.lab  = 1.5)


### Mit null
simulate_2D_both_rect_mat <- function(p, q, inside_pool, outside_pool, zero_pool,
                                      x0, box_limit = 4, n_trials, tol = 1e-12) {
  success_count <- 0L
  total_steps   <- 0L
  
  for (i in seq_len(n_trials)) {
    x <- x0
    steps <- 0L
    
    while (!(all(abs(x) < tol)) && all(abs(x) <= box_limit)) {
      u <- runif(1)
      if (u < p) {
        A <- outside_pool[[ sample.int(length(outside_pool), 1) ]]
      } else if (u < p + q) {
        A <- inside_pool[[ sample.int(length(inside_pool), 1) ]]
      } else {
        A <- zero_pool[[1L]]  # sofortige Nullmatrix
      }
      
      x <- as.vector(A %*% x)
      steps <- steps + 1L
    }
    
    if (any(abs(x) > box_limit)) {
      success_count <- success_count + 1L
    }
    total_steps <- total_steps + steps
  }
  
  list(
    success_prob      = success_count / n_trials,
    expected_duration = as.numeric(total_steps) / n_trials
  )
}

## --- Pools vorbereiten ---
zero_pool <- list(matrix(c(0,0,0,0), 2, 2))   # Nullmatrix

## --- Parameter ---
x0        <- c(1, 1)
box_limit <- 4
p_grid    <- seq(0.00, 1, by = 0.05)
q_grid    <- seq(0.00, 1, by = 0.05)
n_trials  <- 100000  # Beispiel

success_surface  <- matrix(NA_real_, nrow = length(p_grid), ncol = length(q_grid),
                           dimnames = list(p = p_grid, q = q_grid))
duration_surface <- matrix(NA_real_, nrow = length(p_grid), ncol = length(q_grid),
                           dimnames = list(p = p_grid, q = q_grid))

tol_simplex <- 1e-12

## --- Doppelloop über p und q ---
system.time({
  for (ip in seq_along(p_grid)) {
    for (iq in seq_along(q_grid)) {
      p <- p_grid[ip]
      q <- q_grid[iq]
      
      if (p + q <= 1) {
        res <- simulate_2D_both_rect_mat(
          p           = p,
          q           = q,
          inside_pool = inside_pool,
          outside_pool = outside_pool,
          zero_pool   = zero_pool,
          x0          = x0,
          box_limit   = box_limit,
          n_trials    = n_trials
        )
        
        success_surface[ip, iq]  <- res$success_prob
        duration_surface[ip, iq] <- res$expected_duration
      }
    }
    print(p_grid[ip])  # Fortschritt
  }
})

## --- 3D Plots ---
par(mar = c(3, 3, 2, 1) + 0.1)

# Erfolgswahrscheinlichkeit
persp(x = p_grid, y = q_grid, z = success_surface,
      theta = 20, phi = 30, expand = 0.8,
      xlab = "p", ylab = "q", zlab = "P(S_T=b)",
      ticktype = "detailed", col = "deepskyblue3",
      cex.axis = 1.3, cex.lab = 1.3)

# Erwartete Dauer
persp(x = p_grid, y = q_grid, z = duration_surface,
      theta = 20, phi = 20, expand = 0.8,
      xlab = "p", ylab = "q", zlab = expression(E(T)),
      ticktype = "detailed", col = "deepskyblue3",
      cex.axis = 1.3, cex.lab = 1.3)
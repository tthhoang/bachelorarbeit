## --- Deine bestehende Kombi-Funktion (Faktoren, 2D, multiplicativ) ---
simulate_2D_both_rect <- function(p_weights, factors, x0,
                                  box_limit, n_trials,
                                  tol = 1e-12) {
  p_weights <- p_weights / sum(p_weights)  # defensiv normalisieren
  
  success_count <- 0L
  total_steps   <- 0L
  
  for (i in seq_len(n_trials)) {
    x <- x0
    steps <- 0L
    
    # solange NICHT im Ruin (nahe 0) und noch in der Box
    while (!(all(abs(x) < tol)) && all(abs(x) <= box_limit)) {
      a <- sample(factors, size = 2, replace = TRUE, prob = p_weights)
      x <- c(a[1] * x[1], a[2] * x[2])
      steps <- steps + 1L
    }
    
    # Erfolg = mindestens eine Koordinate raus
    if (any(abs(x) > box_limit)) success_count <- success_count + 1L
    total_steps <- total_steps + steps
  }
  
  list(
    success_prob      = success_count / n_trials,
    expected_duration = as.numeric(total_steps) / n_trials
  )
}

# Parameter
x0 <- c(1, 1)
box_limit <- 4
factors <- c(1/3, 0.5, 2, 3)  # Beispiel: manchmal Kontraktion, manchmal Expansion
p_grid <- seq(0.01, 0.99, by = 0.025)
q_grid   <- seq(0.00, 1.00, by = 0.05)
n_trials <- 100000
tol_simplex <- 1e-12

### ohne Null
# Ergebnisse
success_probs <- numeric(length(p_grid))
durations <- numeric(length(p_grid))

system.time({
  for (i in seq_along(p_grid)) {
    p <- p_grid[i]
    probs <- c((1 - p)/2, (1 - p)/2, p/2, p/2)  # Erfolg = Expansion
    res <- simulate_2D_both_rect(probs, factors, x0, box_limit, n_trials)
    success_probs[i] <- res$success_prob
    durations[i]     <- res$expected_duration
    print(p)
  }
})

# Plotten

par(mar = c(5, 5.5, 2, 2) + 0.1) 

plot(p_grid, success_probs, type = "l", lwd = 2, col = "deepskyblue3",
     xlab = expression(p), ylab = expression(P(S[T] == b)), bty = "l", cex.axis = 1.5, cex.lab = 1.4)

plot(p_grid, durations, type = "l", lwd = 2, col = "deepskyblue3",
     xlab = expression(p), ylab = expression(E(T)), bty = "l", cex.axis = 1.5, cex.lab = 1.4)

### mit Null
## Ergebnis-Matrizen (p x q), wie im zweiten Block
Z_success  <- matrix(NA_real_, nrow = length(p_grid), ncol = length(q_grid),
                     dimnames = list(p = p_grid, q = q_grid))
Z_duration <- matrix(NA_real_, nrow = length(p_grid), ncol = length(q_grid),
                     dimnames = list(p = p_grid, q = q_grid))

## --- Doppelloop über (p,q) im Simplex ---
system.time({
  for (ip in seq_along(p_grid)) {
    p <- p_grid[ip]
    for (iq in seq_along(q_grid)) {
      q <- q_grid[iq]
      
      # nur gültige Punkte: 1 - p - q >= 0
      if (p + q > 1 + tol_simplex) next
      
      # Roh-Gewichte passend zur Reihenfolge in 'factors'
      # 0, (zwei <1), (zwei >1)
      probs0 <- c(1 - p - q, q/2, q/2, p/2, p/2)
      
      # numerisch robust: kleine Negative clippen, Summe prüfen & normalisieren
      if (any(!is.finite(probs0))) next
      if (any(probs0 < -tol_simplex)) next
      probs0[probs0 < 0] <- 0
      s <- sum(probs0)
      if (s <= 0) next
      probs <- probs0 / s
      
      # Simulation
      res <- simulate_2D_both_rect(
        p_weights = probs,
        factors   = factors,
        x0        = x0,
        box_limit = box_limit,
        n_trials  = n_trials
      )
      
      Z_success[ip, iq]  <- res$success_prob
      Z_duration[ip, iq] <- res$expected_duration
    }
    print(p)  # Fortschritt
  }
})

## --- 3D-Plots (Base R: persp) ---
par(mar = c(3, 3, 2, 1) + 0.1)

# Erfolgswahrscheinlichkeit
persp(x = p_grid, y = q_grid, z = Z_success,
      theta = 20, phi = 25, expand = 0.7,
      xlab = "p", ylab = "q", zlab = "P(S[T]=b)",
      ticktype = "detailed", col = "deepskyblue3",
      cex.axis = 1.5, cex.lab = 1.4)

# Erwartete Laufzeit
persp(x = p_grid, y = q_grid, z = Z_duration,
      theta = 15, phi = 25, expand = 0.7,
      xlab = "p", ylab = "q", zlab = "E(T)",
      ticktype = "detailed", col = "deepskyblue3",
      cex.axis = 1.5, cex.lab = 1.4)

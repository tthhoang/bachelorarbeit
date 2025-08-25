simulate_both <- function(probs, steps, a, b, n_trials, tol = 1e-12) {
  probs <- probs / sum(probs)              # sicherheitshalber normalisieren
  log_tol <- log(tol)   #untere Grenze
  success_count <- 0L
  total_steps   <- 0L
  
  for (i in seq_len(n_trials)) {
    S_n <- S_0
    steps_count <- 0L
    while (S_n > log_tol && S_n < b) {
      x <- sample(steps, size = 1, prob = probs)
      S_n <- S_n + x
      steps_count <- steps_count + 1L
    }
    if (S_n >= b) success_count <- success_count + 1L
    total_steps <- total_steps + steps_count
  }
  list(
    success_prob      = success_count / n_trials,
    expected_duration = as.numeric(total_steps) / n_trials
  )
}

S_0 <- 0    # Start
b <- log(4)   # obere Grenze
p_grid <- seq(0.00, 1.00, by = 0.025)   # p
q_grid <- seq(0.00, 1.00, by = 0.025)   # für den Fall mit Null
n_trials <- 100000
tol = 1e-12

### Schritte ohne Null
steps <- c(log(0.5), log(2))
# steps <- c(log(1/3),log(0.5), log(2),log(3))
# steps <- c(log(1/5),log(1/4),log(1/3),log(0.5), log(2),log(3))
# steps <- c(log(1/5), log(2),log(3), log(4))

# Simulationen für verschiedene p
success_probs <- numeric(length(p_grid))
durations <- numeric(length(p_grid))

system.time({
  for (i in seq_along(p_grid)) {
    p <- p_grid[i]
    probs <- c((1 - p)/4,(1 - p)/4,(1 - p)/4,(1 - p)/4, p/2, p/2)
    
    res <- simulate_both(probs, steps, a, b, n_trials)
    success_probs[i] <- res$success_prob
    durations[i]     <- res$expected_duration
    print(p)
  }
})

# Plots
par(mar = c(5, 5.5, 2, 2) + 0.1) 
# Erfolgswahrscheinlichkeit
plot(p_grid, success_probs, type = "l", lwd = 2, col = "deepskyblue3",
     xlab = expression(p), ylab = expression(P(S[T] == b)),
     bty = "l",cex.axis = 1.5, cex.lab = 1.5)

# Erwartete Dauer
plot(p_grid, durations, type = "l", lwd = 2, col = "deepskyblue3",
     xlab = expression(p), ylab = expression(E(T)),
     bty = "l",cex.axis = 1.5, cex.lab = 1.4)



### Schritte mit null (am besten immer ein davon auskommentieren)
steps <- c(log(0), log(0.5), log(2))
# steps <- c(log(0),log(1/5),log(1/4),log(1/3),log(0.5), log(2),log(3))

# Ergebnis-Matrizen für 3D-Plot
Z_success  <- matrix(NA_real_, nrow = length(q_grid), ncol = length(p_grid),
                     dimnames = list(q = q_grid, p = p_grid))
Z_duration <- matrix(NA_real_, nrow = length(q_grid), ncol = length(p_grid),
                     dimnames = list(q = q_grid, p = p_grid))

system.time({
  for (ip in seq_along(p_grid)) {
    p <- p_grid[ip]
    for (iq in seq_along(q_grid)) {
      q <- q_grid[iq]
      
      # 1) Nur gültige Punkte im Simplex (mit Toleranz)
      if (p + q > 1) {
        next
      }
      
      # 2) Wahrscheinlichkeiten
      probs0 <- c(1 - p - q, q, p)
      # probs0 <- c(1 - p - q, q/4,q/4,q/4,q/4, p/2,p/2)  
      
      # 3) Winzige negative Rundungsreste auf 0 clampen, harte Negativen skippen
      if (any(probs0 < -tol) || !is.finite(sum(probs0))) {
        next
      }
      probs0[probs0 < 0] <- 0
      
      # 4) Summe prüfen und normalisieren
      s <- sum(probs0)
      if (s <= 0) next
      probs <- probs0 / s
      
      # 5) Simulation
      res <- simulate_both(probs, steps, a, b, n_trials)
      Z_success[ip, iq]  <- res$success_prob
      Z_duration[ip, iq] <- res$expected_duration
    } 
    print(p)}
})

# 3D-Plots
par(mar = c(3, 3, 2, 1) + 0.1)

# 3D: Erfolgswahrscheinlichkeit über (p,q)
res <- persp(x = p_grid, y = q_grid, z = Z_success,
             theta = 20, phi = 25, expand = 0.7,
             xlab = "p", ylab = "q", zlab = "P(S_T=b)",
             ticktype = "detailed", col = "deepskyblue3",
             cex.axis = 1.5, cex.lab = 1.4)


# 3D: Erwartete Laufzeit über (p,q)
persp(x = p_grid, y = q_grid, z = Z_duration,
      theta = 15, phi = 25, expand = 0.7,
      xlab = "p", ylab = "q", zlab = "E(T)",
      ticktype = "detailed", col = "deepskyblue3",cex.axis = 1.5, cex.lab = 1.4)



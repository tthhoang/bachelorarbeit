x0 <- c(1, 1)            # Startpunkt
tau <- 4                 # äußerer Absorptionsrand
a <- 1e-12               # 0
n_trials <- 100000       # Anzahl der Simulationen

# Faktoren und p für den Fall ohne 0
schritte <- c(1/3, 1/2, 2, 3)  
p_grid <- seq(0.01, 0.99, by = 0.025)

# Faktoren, p und q für den Fall mit 0
schritte_0 <- c(0, 1/3, 1/2, 2, 3)  
p_grid_0 <- seq(0.00, 1.00, by = 0.05)
q_grid_0 <- seq(0.00, 1.00, by = 0.05)

# Funktion der Simulation
simulation <- function(schritte, wkt, n_trials, x0, a, tau){
  erfolg_count <- 0
  total_steps <- 0
  
  for(i in seq_len(n_trials)){
    x <- x0
    steps_count <- 0
    while(!(all(abs(x) < a)) && all(abs(x) <= tau)){
      m <- sample(schritte, size = 2, replace = TRUE, pro = wkt)
      x <- c(m[1]*x[1], m[2]*x[2])
      steps_count <- steps_count + 1
    }
    if(any(abs(x) > tau)){
      erfolg_count <- erfolg_count + 1
    }
    total_steps <- total_steps + steps_count
  }
  list(
    erfolg_wkt = erfolg_count/n_trials,
    erwartete_absorptionszeit = total_steps/n_trials
  )
}

### Fall ohne 0 ###
absorptionswkt <- numeric(length(p_grid))
absorptionszeit <- numeric(length(p_grid))

# Die Simulation
system.time({
  for (i in seq_along(p_grid)) {
    p <- p_grid[i]
    wkt <- c((1 - p)/2, (1 - p)/2, p/2, p/2)  # Erfolg = Expansion
    res <- simulation(schritte, wkt, n_trials, x0, a, tau)
    absorptionswkt[i] <- res$erfolg_wkt
    absorptionszeit[i] <- res$erwartete_absorptionszeit
    print(p)
  }
})

## Plots
par(mar = c(5, 5.5, 2, 2) + 0.1) 

# Absorptionswahrscheinlichkeit
plot(p_grid, absorptionswkt, type = "l", lwd = 2, col = "deepskyblue3",
     xlab = expression(p), ylab = expression(P(S[T] == b)), bty = "l", cex.axis = 1.5, cex.lab = 1.4)

# Erwartete Absorptionszeit
plot(p_grid, absorptionszeit, type = "l", lwd = 2, col = "deepskyblue3",
     xlab = expression(p), ylab = expression(E(T)), bty = "l", cex.axis = 1.5, cex.lab = 1.4)




### Fall mit 0 ###
absorptionswkt_mit_0  <- matrix(NA_real_, nrow = length(q_grid_0), ncol = length(p_grid_0),
                                dimnames = list(q = q_grid_0, p = p_grid_0))                                  
absorptionszeit_mit_0 <- matrix(NA_real_, nrow = length(q_grid_0), ncol = length(p_grid_0),
                                dimnames = list(q = q_grid_0, p = p_grid_0))

# Die Simulation
system.time({
  for(ip in seq_along(p_grid_0)){
    p <- p_grid_0[ip]
    for(iq in seq_along(q_grid_0)){
      q <- q_grid_0[iq]
      if (p+q > 1){next}
      wkt_0 <- c(1-p-q, q/2, q/2, p/2, p/2)
      if (any(wkt_0 < -1e-12)) {next}
      wkt_0[wkt_0 < 0] <- 0
      wkt_0 <- wkt_0/sum(wkt_0)
      res <- simulation(schritte_0, wkt_0, n_trials, x0, a, tau)
      absorptionswkt_mit_0[ip,iq] <- res$erfolg_wkt
      absorptionszeit_mit_0[ip,iq] <- res$erwartete_absorptionszeit
    }
    print(p)
  }
})


## 3D-Plots
par(mar = c(3, 3, 2, 1) + 0.1)

# Absorptionswahrscheinlichkeit
persp(x = p_grid_0, y = q_grid_0, z = absorptionswkt_mit_0,
      theta = 20, phi = 25, expand = 0.7,
      xlab = "p", ylab = "q", zlab = "P(S_T=b)",
      ticktype = "detailed", col = "deepskyblue3",
      cex.axis = 1.5, cex.lab = 1.4)

# Erwartete Absorptionszeit
persp(x = p_grid_0, y = q_grid_0, z = absorptionszeit_mit_0,
      theta = 15, phi = 25, expand = 0.7,
      xlab = "p", ylab = "q", zlab = "E(T)",
      ticktype = "detailed", col = "deepskyblue3",cex.axis = 1.5, cex.lab = 1.4)


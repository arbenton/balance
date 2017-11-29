library(forecast)
library(Rcpp)

sourceCpp("balancer.cpp")

saa_periodic_review <- function(x, model, h, p, n, iter=25, lower=0, upper=200) {
  costs <- matrix(NA, nrow = iter, ncol = (upper - lower + 1))
  for (i in 1:iter) {
    samples <- simulate(model, nsim = n)
    costs[i, ] <- pr_cost(lower:upper, x, samples, h[1:n], p[1:n])
  }
  approx <- colMeans(costs)
  return(which.min(approx))
}

n_periods <- 100
n_samples <- 50
n_burn_in <- 10
demand <- cumsum(rnorm(n=n_periods, mean=.1, sd=1)) + 100
holding <- rep(.5, n_periods)
ordering <- rep(3, n_periods)
penalty <- rep(1, n_periods)
inventory <- rep(0, n_periods)
inventory[1:n_burn_in] <- 0
orders <- rep(0, n_periods)

pos <- function(x) {
  ifelse(x > 0, x, 0)
}

for (t in n_burn_in:n_periods) {
  model <- Arima(demand[1:(t-1)], order=c(0, 1, 0), include.drift=T)
  orders[t] <- saa_periodic_review(inventory[t-1], model,
      holding[t:n_periods], penalty[t:n_periods], (n_periods - t + 1))
  inventory[t] <- pos(inventory[t-1] + orders[t] - demand[t])
}

inventory_history <- data.frame(
  order_placed = orders[n_burn_in:n_periods],
  demand_observed = demand[n_burn_in:n_periods],

)
inventory_history$demand_met <- pmin(inventory_history$demand_observed,
  inventory_history$starting_inventory + inventory_history$order_placed)
inventory_history$units_backordered <- inventory_history$demand_observed - inventory_history$demand_met
inventory_history$units_held <- inventory_histol

print(summary(inventory_history))

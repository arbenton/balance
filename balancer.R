library(forecast)
library(Rcpp)

sourceCpp("balancer.cpp")

optimize_dual_cost <- function(x, model, holding, penalty, horizon) {
  samples <- simulate(model, nsim = horizon)
  max_cost <- function(q) {
    balancer_cost(q, x, samples, holding[1:horizon], penalty[1:horizon])
  }
  opt <- optimize(max_cost, c(50, 200))$minimum
  return(opt)
}

set.seed(1)

n_periods <- 150
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

for (t in n_burn_in:(n_periods - 1)) {
  model <- Arima(demand[1:t], order=c(0, 1, 0), include.drift=T)
  order_quantity <- c()
  for (i in 1:n_samples) {
    opt <- optimize_dual_cost(inventory[t-1], model,
      holding[t:n_periods], penalty[t:n_periods], (n_periods - t + 1))
    order_quantity <- c(order_quantity, opt)
  }
  orders[t] <- mean(order_quantity)
  inventory[t] <- pos(inventory[t-1] + orders[t] - demand[t])
}

inventory_summary <- data.frame(
  starting_inventory = inventory[n_burn_in:n_periods],
  order_placed = orders[n_burn_in:n_periods],
  demand_observed = demand[n_burn_in:n_periods]
)

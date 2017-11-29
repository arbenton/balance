import numpy as np
from scipy import optimize, stats, integrate

def holding_cost(quantity, inventory, demand, holding_per):
    cost = cum_demand = 0
    for d, h in zip(demand, holding_per):
        cum_demand += d
        n_left = max(cum_demand - inventory, 0)
        if quantity < passive_hold:
            break
        n_held = quantity - n_left
        cost += n_held*h
    return cost

def backorder_cost(quantity, inventory, demand, penalty_per):
    n_back = max(demand[0] - (inventory + quantity), 0)
    cost = penalty_per*n_back
    return cost

def dual_cost(quantity, inventory, demand, holding_per, penalty_per):
    return max(holding_cost(quantity, inventory, demand, holding_per),
               backorder_cost(quantity, inventory, demand, penalty_per))

def optimal_periodic_review(quantity, inventory, demand, holding_per, penalty_per)

if __name__ == "__main__":

import numpy as np
from scipy import stats, optimize

def backorder_loss(q, x, d_mean, d_std, p_cost):
    dist = stats.norm(d_mean[0], d_std[0])
    penalty = dist.expect(lambda u: p_cost*(u - x - q), lb=(x+q))
    return penalty

def der_backorder_loss(q, x, d_mean, d_std, p_cost):
    dist = stats.norm(d_mean[0], d_std[0])
    der_penalty = p_cost*(1 - dist.cdf(x + q))
    return der_penalty

def holding_loss(q, x, d_mean, d_std, h_cost, tol=1e-6):
    holding = 0
    for mean, var in zip(np.cumsum(d_mean), np.cumsum(np.square(d_std))):
        dist = stats.norm(mean, np.sqrt(var))
        held = dist.expect(lambda u: h_cost*(x + q - u), lb=x, ub=(x+q))
        if held < tol:
            break
        holding += held
    return holding

def der_holding_loss(q, x, d_mean, d_std, h_cost, tol=1e-6):
    der_holding = 0
    for mean, var in zip(np.cumsum(d_mean), np.cumsum(np.square(d_std))):
        dist = stats.norm(mean, np.sqrt(var))
        held = dist.cdf(x + q) - dist.cdf(x)
        der_holding += h_cost*held
    return der_holding

def dual_root(q, x, d_mean, d_std, h, p):
    return holding_loss(q, x, d_mean, d_std, h) \
        - backorder_loss(q, x, d_mean, d_std, p)

def dual_root_der(q, x, d_mean, d_std, h, p):
    return der_holding_loss(q, x, d_mean, d_std, h) \
        + der_backorder_loss(q, x, d_mean, d_std, p)

def dual_max(q, x, d_mean, d_std, h, p):
    return max(holding_loss(q, x, d_mean, d_std, h),
               backorder_loss(q, x, d_mean, d_std, p))

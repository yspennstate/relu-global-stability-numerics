# -*- coding: utf-8 -*-
# HISTORICAL COMPARISON: this file evaluates the earlier unoptimized h bound.
# Current publication tables and finite certificates use finite_budget.py (optimized g).
"""The exponent budget of the proof, evaluated exactly at finite depth.

The master inequality gives, for every width n, depth L >= 2 and theta in (0, 1/4),

    P(K_{n,L}(alpha) > 1)  <=  2 exp( n * B_L(alpha, theta) ),
    B_L(alpha, theta) = L H(1/L) + log 5 + theta log 4 + L h_alpha(theta),

with H the binary entropy and h_alpha the finite-width scalar moment bound (paper (47)).  The threshold
the proof itself delivers at depth L is therefore

    alpha*_L = sup{ alpha : min_theta B_L(alpha, theta) < 0 },

a number one can compute, and the honest finite-depth row of any table.  The closed form
2 exp(-sqrt(10 (log L + A)/L)) is the large-L asymptotics of the same budget with the o(1) terms dropped;
at L = 100 it reads 0.83 where the budget delivers 0.68.  (The closed form is the
one that tends to get quoted as if it were the finite-depth number; the tables use the budget.)
"""
import math
import numpy as np

_THETAS = np.linspace(1e-4, 0.2499, 6000)


def H(u):
    if u <= 0.0 or u >= 1.0:
        return 0.0
    return -u * math.log(u) - (1.0 - u) * math.log(1.0 - u)


def h(alpha, theta):
    return theta * math.log(alpha / (2.0 * math.e)) + math.log((1.0 + (1.0 - 4.0 * theta) ** -0.5) / 2.0)


def budget(L, alpha, theta):
    return L * H(1.0 / L) + math.log(5.0) + theta * math.log(4.0) + L * h(alpha, theta)


def rate(L, alpha):
    """(c, theta*): the largest exponential width rate the budget gives, and the theta that gives it."""
    th = _THETAS
    vals = -(L * H(1.0 / L) + math.log(5.0) + th * math.log(4.0)
             + L * (th * math.log(alpha / (2.0 * math.e)) + np.log((1.0 + (1.0 - 4.0 * th) ** -0.5) / 2.0)))
    k = int(np.argmax(vals))
    return float(vals[k]), float(th[k])


def alpha_star(L, it=50):
    """The largest alpha at which the budget is negative for some theta: the proof's own threshold at depth L."""
    lo, hi = 1e-4, 1.999999
    if rate(L, lo)[0] <= 0.0:
        return 0.0
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if rate(L, mid)[0] > 0.0:
            lo = mid
        else:
            hi = mid
    return lo


def depth_needed(alpha, lo=2, hi=10 ** 8):
    """The smallest depth at which the budget turns negative for this alpha (the L0 of the theorem)."""
    if rate(hi, alpha)[0] <= 0.0:
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if rate(mid, alpha)[0] > 0.0:
            hi = mid
        else:
            lo = mid
    return hi


if __name__ == "__main__":
    for L in (10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 10):
        print("L=%-12d alpha*_L=%.4f   closed form 2exp(-sqrt(10 log L/L))=%.4f" % (L, alpha_star(L), 2 * math.exp(-math.sqrt(10 * math.log(L) / L))))
    for a in (1.0, 1.4, 1.7, 1.9):
        print("alpha=%.2f  L0=%s" % (a, depth_needed(a)))

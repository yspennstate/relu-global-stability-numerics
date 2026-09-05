# -*- coding: utf-8 -*-
# HISTORICAL COMPARISON: this file evaluates the earlier unoptimized h bound.
# Current publication tables and finite certificates use finite_budget.py (optimized g).
"""The numbers the finite-width paper quotes, computed here and nowhere else.

    h_alpha(theta) = theta log(alpha/(2e)) + log((1 + (1-4 theta)^{-1/2})/2)          (paper (47))
    B_L(alpha, theta) = L H(1/L) + log 5 + theta log 4 + L h_alpha(theta)             the exponent budget
    P(K_{n,L}(alpha) > 1) <= 2 exp(n B_L(alpha, theta))     for every n >= 1, L >= 2, theta in (0, 1/4)
    c_L(alpha) = max_theta [-B_L(alpha, theta)]            the width rate; alpha*_L = sup{alpha : c_L(alpha) > 0}
"""
import math, sys
import numpy as np
from math import comb, log, exp, sqrt

def H(u):
    if u <= 0 or u >= 1:
        return 0.0
    return -u * log(u) - (1 - u) * log(1 - u)

def h(alpha, th):
    return th * log(alpha / (2 * math.e)) + log((1 + (1 - 4 * th) ** -0.5) / 2)

def B(L, alpha, th):
    return L * H(1.0 / L) + log(5) + th * log(4) + L * h(alpha, th)

THETAS = np.linspace(1e-4, 0.2499, 40000)

def c_rate(L, alpha):
    vals = np.array([-B(L, alpha, t) for t in THETAS])
    k = int(np.argmax(vals))
    return float(vals[k]), float(THETAS[k])

def alpha_star(L, lo=1e-3, hi=1.999999, it=60):
    if c_rate(L, lo)[0] <= 0:
        return 0.0
    for _ in range(it):
        mid = (lo + hi) / 2
        if c_rate(L, mid)[0] > 0:
            lo = mid
        else:
            hi = mid
    return lo

def C_exact(n, L):
    return 2 * sum(comb(n * L - 1, j) for j in range(n))

def main():
    out = []
    P = out.append
    P("%% finite-depth computable thresholds alpha*_L and the width rate at a few alphas")
    P("%% L & alpha*_L & 2exp(-sqrt(10 log L / L)) & 2exp(-sqrt(10(log L + A)/L)), A = 1 + log 5 + 1")
    A = 1 + log(5) + 1.0
    for L in (2, 5, 10, 30, 100, 300, 1000, 3000, 10000, 100000, 1000000):
        a = alpha_star(L)
        q1 = 2 * exp(-sqrt(10 * log(L) / L))
        q2 = 2 * exp(-sqrt(10 * (log(L) + A) / L))
        P("%d & %.4f & %.4f & %.4f \\\\" % (L, a, q1, q2))
    P("")
    P("%% the width rate c_L(alpha) and the width n_0.01 at which 2 exp(-c n) <= 0.01, i.e. n >= log(200)/c")
    for alpha in (0.5, 1.0, 1.5, 1.8, 1.9, 1.95):
        row = []
        for L in (100, 300, 1000, 3000, 10000, 100000):
            c, th = c_rate(L, alpha)
            if c > 0:
                row.append("%.3f (%d)" % (c, math.ceil(log(200) / c)))
            else:
                row.append("--")
        P("%.2f & " % alpha + " & ".join(row) + " \\\\")
    P("")
    P("%% explicit failure bounds 2exp(-c n) at alpha = 1, a few depths and widths")
    for L in (300, 1000, 10000):
        c, th = c_rate(L, 1.0)
        P("L=%d theta*=%.4f c=%.4f: " % (L, th, c) + ", ".join("n=%d: %.2e" % (n, 2 * exp(-c * n)) for n in (10, 20, 50, 100, 1000)))
    P("")
    P("%% the alpha > 2 side: 5 alpha^2 L / (n (alpha-2)^2)")
    for alpha in (2.5, 3.0, 4.0):
        P("alpha=%.1f: " % alpha + ", ".join("L=%d,n=%d: %.3g" % (L, n, 5 * alpha ** 2 * L / (n * (alpha - 2) ** 2)) for L in (10, 100) for n in (100, 1000, 10000)))
    P("")
    P("%% C_{n,L} exact against 2 exp(nL H(1/L)) and against 2^{nL}")
    for n, L in ((2, 2), (3, 3), (5, 3), (10, 5), (20, 10), (50, 10)):
        C = C_exact(n, L)
        P("n=%d L=%d: C=%s  log C/n=%.4f  L H(1/L)=%.4f  L log2=%.4f" % (n, L, ("%.3e" % C) if C > 1e6 else str(C), log(C) / n, L * H(1 / L), L * log(2)))
    P("")
    P("%% all-active realizability P(R_D != empty) = 2^{1-nL} sum_{j<n} C(nL-1, j)")
    for n, L in ((2, 2), (3, 2), (3, 3), (5, 3), (10, 3), (20, 3), (10, 10)):
        p = C_exact(n, L) / 2.0 ** (n * L)
        P("n=%d L=%d: %.6g" % (n, L, p))
    P("")
    P("%% the finite h against the exact exponent g at alpha=2: h_2(theta) and (5/2) theta^2")
    for th in (0.001, 0.01, 0.05, 0.1, 0.2):
        P("theta=%.3f: h_2=%.6f  (5/2)theta^2=%.6f  h_1=%.6f" % (th, h(2.0, th), 2.5 * th * th, h(1.0, th)))
    print("\n".join(out))

if __name__ == "__main__":
    main()

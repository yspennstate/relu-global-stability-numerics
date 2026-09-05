# -*- coding: utf-8 -*-
"""Numbers for the finite-horizon section and for the way the two results meet.

1. The shared-matrix ReLU dynamics x_{t+1} = (W x_t)_+, W iid N(0, 1/N): the worst deviation
   max_{t<=T} |2^t ||x_t||^2 - 1| over T = 8 steps, at several widths, 200 draws each - and the same
   statistic for the feedforward pass with a FRESH matrix at every step (the feedforward network of section 1 at alpha = 1).
2. The Chernoff rate I(c) = sup_{s>0} [ -s c - log((1 + (1+2s)^{-1/2})/2) ] for P(S_n <= n c) with
   S_n = sum of n copies of B Z^2 = (xi_+)^2: the exponential form of the alpha > 2 obstruction.
"""
import math
import numpy as np

rng = np.random.default_rng(20260905)


def worst_dev(N, T, shared, draws=200):
    out = []
    for _ in range(draws):
        x = np.zeros(N); x[0] = 1.0
        W = rng.normal(0.0, 1.0 / math.sqrt(N), size=(N, N))
        worst = 0.0
        for t in range(1, T + 1):
            if not shared:
                W = rng.normal(0.0, 1.0 / math.sqrt(N), size=(N, N))
            x = np.maximum(W @ x, 0.0)
            worst = max(worst, abs(2.0 ** t * float(x @ x) - 1.0))
        out.append(worst)
    return np.array(out)


def chernoff_rate(c):
    s = np.linspace(1e-4, 60.0, 600000)
    vals = -s * c - np.log((1.0 + (1.0 + 2.0 * s) ** -0.5) / 2.0)
    k = int(np.argmax(vals))
    return float(vals[k]), float(s[k])


def main():
    print("%% worst deviation max_t |2^t ||x_t||^2 - 1|, T = 8, 200 draws: median / 90th pct / fraction > 0.5")
    print("%% N & shared W & fresh W per step")
    for N in (20, 50, 100, 500, 2000):
        a = worst_dev(N, 8, True)
        b = worst_dev(N, 8, False)
        print("%d & %.3f / %.3f / %.3f & %.3f / %.3f / %.3f \\\\" % (
            N, np.median(a), np.quantile(a, 0.9), np.mean(a > 0.5), np.median(b), np.quantile(b, 0.9), np.mean(b > 0.5)))
    print()
    print("%% Chernoff rate I(1/alpha) for P(K <= 1) <= L exp(-n I(1/alpha)), and the Chebyshev bound 5 alpha^2 L/(n(alpha-2)^2) at L=10")
    for alpha in (2.2, 2.5, 3.0, 4.0, 6.0):
        I, s = chernoff_rate(1.0 / alpha)
        row = ["alpha=%.1f: I=%.5f (s*=%.2f)" % (alpha, I, s)]
        for n in (100, 1000, 10000):
            row.append("n=%d: chernoff %.2e, chebyshev %.2e" % (n, 10 * math.exp(-n * I), 5 * alpha ** 2 * 10 / (n * (alpha - 2) ** 2)))
        print("   ".join(row))


if __name__ == "__main__":
    main()

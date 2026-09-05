# -*- coding: utf-8 -*-
"""Chapter six's measured numbers: the shared-matrix ReLU dynamics x_{t+1} = (W x_t)_+ with W iid N(0, 1/N),
against the same statistic for a fresh matrix at every step (chapter one's forward pass at alpha = 1).

    worst_t |2^t ||x_t||^2 - 1| over T = 8 steps, 200 draws per width: median, 90th percentile, fraction > 1/2
    the Chernoff rate I(1/alpha) of the sharpened obstruction (section 'together' of the paper)

Writes sim/data/s8_horizon.json.  Seed fixed; every number in the lecture and the paper comes from here.
"""
import io, json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
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
    rows = []
    for N in (20, 50, 100, 500, 2000):
        a = worst_dev(N, 8, True)
        b = worst_dev(N, 8, False)
        rows.append({"N": N, "T": 8, "draws": 200,
                     "shared": {"median": float(np.median(a)), "p90": float(np.quantile(a, 0.9)), "frac_gt_half": float(np.mean(a > 0.5))},
                     "fresh": {"median": float(np.median(b)), "p90": float(np.quantile(b, 0.9)), "frac_gt_half": float(np.mean(b > 0.5))}})
    cher = []
    for alpha in (2.2, 2.5, 3.0, 4.0, 6.0):
        I, s = chernoff_rate(1.0 / alpha)
        cher.append({"alpha": alpha, "I": I, "s_star": s,
                     "chernoff_L10_n1000": 10 * math.exp(-1000 * I), "chebyshev_L10_n1000": 5 * alpha ** 2 * 10 / (1000 * (alpha - 2) ** 2)})
    out = {"note": "measured; numpy default_rng(20260905); T = 8; worst deviation max_t |2^t||x_t||^2 - 1|",
           "worst_dev": rows, "chernoff": cher}
    json.dump(out, io.open(os.path.join(HERE, "data", "s8_horizon.json"), "w", encoding="utf-8"), indent=1)
    for r in rows:
        print(r["N"], "shared %.3f/%.3f/%.3f" % tuple(r["shared"].values()), "fresh %.3f/%.3f/%.3f" % tuple(r["fresh"].values()))
    for c in cher:
        print(c)


if __name__ == "__main__":
    main()

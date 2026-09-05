# -*- coding: utf-8 -*-
"""Simulation 1 - the forward radial law.

For a FIXED x, the coordinates of Wx are independent N(0, alpha||x||^2/n), so
    ||phi(Wx)||^2 / ||x||^2  =d  (alpha/n) sum_i B_i Z_i^2  =  X_n^(alpha),
independently of the direction of x.  That identity is what we sample, and it is exactly the
scalar variable of Section 10.4.  Everything printed below is measured.
"""
import io, json, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20260904)


def X(n, alpha, trials):
    """trials independent copies of X_n^(alpha) = (alpha/n) sum B_i Z_i^2."""
    Z = rng.normal(size=(trials, n))
    B = rng.integers(0, 2, size=(trials, n))
    return (alpha / n) * np.sum(B * Z * Z, axis=1)


def forward_profile(n, L, alpha, trials):
    prof = np.empty((trials, L + 1))
    s = np.sqrt(alpha / n)
    for t in range(trials):
        x = np.zeros(n); x[0] = 1.0
        prof[t, 0] = 1.0
        for l in range(L):
            x = np.maximum(rng.normal(0.0, s, size=(n, n)) @ x, 0.0)
            prof[t, l + 1] = float(x @ x)
    return prof


res = {"note": "measured, numpy default_rng(20260904)"}

res["one_layer"] = {}
for alpha in (1.0, 2.0, 3.0):
    row = {}
    for n in (2, 5, 20, 100, 1000, 10000):
        r = X(n, alpha, 40000)
        row[str(n)] = {"mean": float(r.mean()), "sd": float(r.std()),
                       "p05": float(np.quantile(r, .05)), "p95": float(np.quantile(r, .95)),
                       "p_above_1": float((r > 1).mean()), "trials": int(r.size)}
    res["one_layer"][str(alpha)] = {"target_alpha_over_2": alpha / 2.0, "by_n": row}

for nn in (2, 20, 200):
    h = X(nn, 2.0, 200000)
    lo, hi = (0.0, 4.0) if nn == 2 else ((0.2, 2.0) if nn == 20 else (0.6, 1.5))
    cnt, edges = np.histogram(h, bins=60, range=(lo, hi))
    res["hist_n%d_a2" % nn] = {"counts": cnt.tolist(), "edges": [round(float(e), 5) for e in edges],
                               "mean": float(h.mean()), "sd": float(h.std()), "trials": int(h.size)}

res["profiles"] = {}
for alpha in (1.6, 2.0, 2.6):
    P = forward_profile(200, 8, alpha, 300)
    res["profiles"][str(alpha)] = {
        "median_sq_norm": [float(v) for v in np.median(P, axis=0)],
        "predicted": [float((alpha / 2.0) ** l) for l in range(9)],
        "q10": [float(v) for v in np.quantile(P, .10, axis=0)],
        "q90": [float(v) for v in np.quantile(P, .90, axis=0)],
        "n": 200, "trials": 300}

B = rng.integers(0, 2, size=4_000_000); Z = rng.normal(size=4_000_000)
Y = B * Z ** 2
res["summand_Y"] = {"mean": float(Y.mean()), "var": float(Y.var()), "EY2": float((Y ** 2).mean()),
                    "exact_mean": 0.5, "exact_var": 1.25, "exact_EY2": 1.5,
                    "var_over_2mean2": float(Y.var() / (2 * Y.mean() ** 2)), "exact_ratio": 2.5,
                    "samples": int(Y.size)}
res["gaussian_moments"] = {"EZ2": float((Z ** 2).mean()), "EZ4": float((Z ** 4).mean()),
                           "exact_EZ2": 1.0, "exact_EZ4": 3.0}

json.dump(res, io.open(os.path.join(OUT, "s1_forward.json"), "w", encoding="utf-8"), indent=1)
o = res["one_layer"]["2.0"]["by_n"]
print("alpha=2:  n=2 mean %.4f sd %.4f | n=200.. n=10000 mean %.4f sd %.4f  (target 1)"
      % (o["2"]["mean"], o["2"]["sd"], o["10000"]["mean"], o["10000"]["sd"]))
print("Y=BZ^2 : mean %.5f (0.5)  var %.5f (1.25)  var/2mean^2 %.5f (2.5)"
      % (res["summand_Y"]["mean"], res["summand_Y"]["var"], res["summand_Y"]["var_over_2mean2"]))
print("EZ2 %.4f  EZ4 %.4f" % (res["gaussian_moments"]["EZ2"], res["gaussian_moments"]["EZ4"]))
print("profile alpha=2.6 median:", [round(v, 3) for v in res["profiles"]["2.6"]["median_sq_norm"]])
print("        predicted      :", [round(v, 3) for v in res["profiles"]["2.6"]["predicted"]])
print("wrote s1_forward.json")

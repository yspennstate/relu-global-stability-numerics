# -*- coding: utf-8 -*-
"""Simulation 4 - the scalar high moment, its two rate functions, and the exponent budget.

  h_alpha(theta)  = theta*log(alpha/(2e)) + log((1+(1-4theta)^{-1/2})/2)        (Lemma 12.2, finite n bound)
  Phi_{a,th}(u)   = H(u) - log2 + th*log(2a) + (u/2+th)log(u/2+th) - (u/2)log(u/2) - th    (eq. 51)
  g_alpha(theta)  = max_u Phi                                                    (Prop 13.1, exact rate)

Checks: g'(0+) = h'(0+) = log(alpha/2); at alpha=2 both are (5/2)theta^2 + O(theta^3); the
empirical (1/n)log E[X^{theta n}] converges to g; the maximiser u_theta = 1/2 + theta/2 + O(theta^2).
"""
import io, json, math, os
import numpy as np
from scipy.special import gammaln

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)


def h(alpha, th):
    return th * math.log(alpha / (2 * math.e)) + math.log((1 + (1 - 4 * th) ** -0.5) / 2)


def H(u):
    if u <= 0 or u >= 1:
        return 0.0
    return -(u * math.log(u) + (1 - u) * math.log(1 - u))


def xlx(x):
    return 0.0 if x <= 0 else x * math.log(x)


def Phi(alpha, th, u):
    return (H(u) - math.log(2) + th * math.log(2 * alpha)
            + (u / 2 + th) * math.log(u / 2 + th) - xlx(u / 2) - th)


def g(alpha, th, ngrid=400001):
    u = np.linspace(1e-12, 1 - 1e-12, ngrid)
    Hs = -(u * np.log(u) + (1 - u) * np.log(1 - u))
    v = u / 2 + th
    val = Hs - math.log(2) + th * math.log(2 * alpha) + v * np.log(v) - (u / 2) * np.log(u / 2) - th
    i = int(np.argmax(val))
    return float(val[i]), float(u[i])


def logEXthetan_exact(alpha, th, n):
    """(1/n) log E[(X_n^alpha)^{theta n}] from the exact Gamma sum (eq. 54), in log space."""
    p = th * n
    ms = np.arange(1, n + 1)
    terms = (gammaln(n + 1) - gammaln(ms + 1) - gammaln(n - ms + 1)
             + p * math.log(2 * alpha / n)
             + gammaln(ms / 2.0 + p) - gammaln(ms / 2.0))
    mx = terms.max()
    return (mx + math.log(np.exp(terms - mx).sum()) - n * math.log(2)) / n


res = {"note": "computed exactly (no sampling except where stated)"}

# (a) the two rate functions side by side, and the derivative at zero
res["rates"] = {}
for alpha in (1.2, 1.6, 1.9, 2.0, 2.2, 2.6):
    ths = [float(t) for t in np.unique(np.concatenate([np.geomspace(1e-5, 0.002, 40), np.arange(0.002, 0.2401, 0.002)]))]
    hs = [h(alpha, t) for t in ths]
    gs = [g(alpha, t)[0] for t in ths]
    res["rates"][str(alpha)] = {
        "theta": ths, "h": hs, "g": gs,
        "log_alpha_over_2": math.log(alpha / 2),
        "h_slope_at_0": h(alpha, 1e-7) / 1e-7,
        "g_slope_at_0": g(alpha, 1e-7)[0] / 1e-7,
        "best_theta_h": float(ths[int(np.argmin(hs))]), "min_h": float(min(hs)),
        "best_theta_g": float(ths[int(np.argmin(gs))]), "min_g": float(min(gs))}

# (b) alpha = 2 : the 5/2 coefficient, read off the curve
res["at_alpha_2"] = {"theta": [], "h": [], "g": [], "h_over_theta2": [], "g_over_theta2": []}
for t in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1):
    hh, gg = h(2.0, t), g(2.0, t)[0]
    res["at_alpha_2"]["theta"].append(t)
    res["at_alpha_2"]["h"].append(hh); res["at_alpha_2"]["g"].append(gg)
    res["at_alpha_2"]["h_over_theta2"].append(hh / t ** 2)
    res["at_alpha_2"]["g_over_theta2"].append(gg / t ** 2)

# (c) the maximiser drifts:  u_theta = 1/2 + theta/2 + O(theta^2)
res["maximiser"] = [{"theta": t, "u": g(2.0, t)[1], "predicted": 0.5 + t / 2.0}
                    for t in (0.0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2)]

# (d) does the finite-n exact moment really converge to g ?
res["convergence_to_g"] = []
for alpha in (1.6, 2.0):
    for th in (0.05, 0.1):
        row = {"alpha": alpha, "theta": th, "g": g(alpha, th)[0], "by_n": {}}
        for n in (10, 50, 200, 1000, 5000, 20000):
            row["by_n"][str(n)] = logEXthetan_exact(alpha, th, n)
        row["h_bound"] = h(alpha, th)
        res["convergence_to_g"].append(row)

# (e) the exponent budget:  LH(1/L) + log5 + theta*log4 + L*h_alpha(theta)  -- when does it go negative?
res["budget"] = {}
for alpha in (1.0, 1.4, 1.7, 1.9, 1.95, 1.99):
    # theta is searched on a grid, step 0.0005 above 0.002.  Checked against a golden-section search
    # on the same h: the optimum moves by less than 0.0002 and L0 is IDENTICAL for alpha = 1.0, 1.4
    # and 1.7 (237, 897, 4590).  It moves at the top end, where the budget curve is nearly flat:
    # 52 862 -> 52 841 at alpha = 1.9, 237 405 -> 237 381 at 1.95, 7 353 765 -> 7 353 316 at 1.99 -
    # four parts in ten thousand, and always in the conservative direction.  The grid value is kept:
    # the tables were typeset against it, and a shift of one figure in the fourth place between
    # two tables would be worse than the resolution itself.
    ths = np.unique(np.concatenate([np.geomspace(1e-6, 0.002, 200), np.arange(0.002, 0.2401, 0.0005)]))
    hv = np.array([h(alpha, t) for t in ths])
    i = int(np.argmin(hv)); th_best, gamma = float(ths[i]), float(-hv[i])
    rows, L0 = [], None
    for L in list(range(1, 60)) + [80, 120, 200, 400, 800, 1500, 3000, 6000, 12000, 30000, 10**5, 10**6, 10**7]:
        pos = (L * H(1.0 / L) if L >= 2 else math.log(2)) + math.log(5) + th_best * math.log(4)
        tot = pos - L * gamma
        rows.append({"L": L, "positive": pos, "negative": -L * gamma, "total": tot})
        if L0 is None and tot < 0:
            L0 = L
    # L0 taken off the row grid is the first SAMPLED depth below zero, not the crossing: the grid
    # jumps 800 -> 1500, so alpha=1.4 reported 1500 where the exponent actually turns at 897.  Solve
    # it on the continuous budget instead and keep the grid value only as a labelled coarse figure.
    def _tot(L):
        p = (L * H(1.0 / L) if L >= 2 else math.log(2)) + math.log(5) + th_best * math.log(4)
        return p - L * gamma
    L0_exact = None
    if L0 is not None:
        lo, hi = 1, int(L0)
        if _tot(hi) < 0:
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if _tot(mid) > 0:
                    lo = mid
                else:
                    hi = mid
            L0_exact = hi
    res["budget"][str(alpha)] = {"theta_star": th_best, "gamma": gamma,
                                 "L0": L0_exact if L0_exact else L0,
                                 "L0_grid": L0, "rows": rows}

# (f) the quantitative curve  alpha(L) = 2 exp(-sqrt(10(logL+A)/L))
A = 1 + math.log(5) + 0.5
res["quantitative"] = {"A": A, "rows": [
    {"L": L, "delta": math.sqrt(10 * (math.log(L) + A) / L),
     "alpha_lower": 2 * math.exp(-math.sqrt(10 * (math.log(L) + A) / L)),
     "sqrt10_logL_over_L": math.sqrt(10 * math.log(L) / L)}
    for L in (10, 30, 100, 300, 1000, 10000, 10 ** 5, 10 ** 6, 10 ** 8, 10 ** 10)]}

json.dump(res, io.open(os.path.join(OUT, "s4_rates.json"), "w", encoding="utf-8"), indent=1)
print("SLOPE AT ZERO  (should equal log(alpha/2))")
for a, v in res["rates"].items():
    print("  alpha=%-4s  log(a/2)=%+.5f   h'(0)=%+.5f   g'(0)=%+.5f" %
          (a, v["log_alpha_over_2"], v["h_slope_at_0"], v["g_slope_at_0"]))
print("\nAT alpha = 2   (h and g over theta^2 should tend to 5/2 = 2.5)")
for i, t in enumerate(res["at_alpha_2"]["theta"]):
    print("  theta=%-6g  h/th^2=%.4f   g/th^2=%.4f" %
          (t, res["at_alpha_2"]["h_over_theta2"][i], res["at_alpha_2"]["g_over_theta2"][i]))
print("\nMAXIMISER u_theta   (predicted 1/2 + theta/2)")
for m in res["maximiser"]:
    print("  theta=%-5g  u=%.5f   predicted %.5f" % (m["theta"], m["u"], m["predicted"]))
print("\nEXACT MOMENT -> g")
for r in res["convergence_to_g"]:
    ns = r["by_n"]
    print("  alpha=%.1f theta=%.2f  g=%+.5f  h=%+.5f  n=10 %+.5f  n=1000 %+.5f  n=20000 %+.5f"
          % (r["alpha"], r["theta"], r["g"], r["h_bound"], ns["10"], ns["1000"], ns["20000"]))
print("\nDEPTH AT WHICH THE BUDGET TURNS NEGATIVE")
for a, v in res["budget"].items():
    print("  alpha=%-5s theta*=%.3f  gamma=%.5f  L0=%s" % (a, v["theta_star"], v["gamma"], v["L0"]))
print("\nQUANTITATIVE LOWER BOUND")
for r in res["quantitative"]["rows"]:
    print("  L=%-10d delta=%.4f  alpha_L >= %.4f" % (r["L"], r["delta"], r["alpha_lower"]))
print("wrote s4_rates.json")

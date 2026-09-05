# -*- coding: utf-8 -*-
"""Simulation 6 - an audit of the paper's own identities, by direct Monte Carlo.

Three things are tested against brute force:
  (a) Lemma 10.1   E ||J_D v||^{2p} = prod_l E (alpha chi^2_{m_l} / n)^p     for a FIXED mask
  (b) the mixture identity (37)      sum_m C(n,m) E(a chi2_m/n)^p = 2^n E[X_n^p]
  (c) Proposition 11.1               E K^{2p} <= C 5^n 4^p (E X^p)^L         both sides computed
"""
import io, json, math, os
import numpy as np
from scipy.special import gammaln

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(90210)


def chi2_moment(m, r, alpha, n):
    """E[(alpha/n * chi^2_m)^r] exactly, in log space."""
    if m == 0:
        return 0.0
    return math.exp(r * math.log(alpha / n) + r * math.log(2.0)
                    + gammaln(m / 2.0 + r) - gammaln(m / 2.0))


def cover(n, L):
    return 2 * sum(math.comb(n * L - 1, j) for j in range(n))


res = {"note": "brute-force Monte Carlo against the paper's exact formulae"}

# ---- (a) the fixed-mask product formula ------------------------------------------------------
res["lemma_10_1"] = []
for (n, L, masks, p, trials) in ((4, 3, [[1, 1, 0, 1], [1, 0, 0, 1], [1, 1, 1, 0]], 1.0, 200000),
                                 (4, 3, [[1, 1, 0, 1], [1, 0, 0, 1], [1, 1, 1, 0]], 2.0, 200000),
                                 (6, 2, [[1, 0, 1, 1, 0, 1], [0, 1, 1, 0, 1, 1]], 1.5, 200000)):
    alpha = 1.3
    v = np.zeros(n); v[0] = 1.0
    Ds = [np.diag(np.array(m, dtype=float)) for m in masks]
    acc = 0.0
    s = math.sqrt(alpha / n)
    B = 2000
    vals = []
    for _ in range(trials // B):
        x = np.tile(v, (B, 1))
        for l in range(L):
            W = rng.normal(0.0, s, size=(B, n, n))
            x = np.einsum('bij,bj->bi', W, x)
            x = x * np.array(masks[l], dtype=float)
        vals.append(np.sum(x * x, axis=1) ** p)
    emp = float(np.concatenate(vals).mean())
    exact = 1.0
    for m in masks:
        exact *= chi2_moment(int(sum(m)), p, alpha, n)
    res["lemma_10_1"].append({"n": n, "L": L, "p": p, "alpha": alpha, "masks": masks,
                              "empirical": emp, "exact_product": exact,
                              "ratio": emp / exact if exact else None, "trials": trials})

# ---- (b) the mixture identity ------------------------------------------------------------------
res["mixture"] = []
for (n, p, alpha) in ((6, 1.0, 1.3), (6, 2.5, 1.3), (12, 1.7, 2.0)):
    lhs = sum(math.comb(n, m) * chi2_moment(m, p, alpha, n) for m in range(n + 1))
    Z = rng.normal(size=(400000, n)); Bn = rng.integers(0, 2, size=(400000, n))
    X = (alpha / n) * np.sum(Bn * Z * Z, axis=1)
    rhs_emp = float((X ** p).mean()) * 2 ** n
    res["mixture"].append({"n": n, "p": p, "alpha": alpha, "sum_side": lhs,
                           "2n_times_EXp_empirical": rhs_emp,
                           "ratio": rhs_emp / lhs})

# ---- (c) the master inequality, both sides ------------------------------------------------------
def K_sweep_n2(Ws, ngrid=None):
    """Enumerate planar cells; ngrid is ignored for compatibility."""
    from region_geometry import planar_network
    return planar_network(Ws)["K"]


res["master"] = []
for (n, L, alpha, p, trials) in ((2, 2, 1.0, 1.0, 3000), (2, 3, 1.0, 1.5, 2000),
                                 (2, 2, 0.6, 2.0, 3000)):
    Ks = []
    s = math.sqrt(alpha / n)
    for _ in range(trials):
        Ws = [rng.normal(0.0, s, size=(n, n)) for _ in range(L)]
        Ks.append(K_sweep_n2(Ws, 6000))
    lhs = float(np.mean(np.array(Ks) ** (2 * p)))
    EXp = sum(math.comb(n, m) * chi2_moment(m, p, alpha, n) for m in range(n + 1)) / 2 ** n
    rhs = cover(n, L) * 5 ** n * 4 ** p * EXp ** L
    res["master"].append({"n": n, "L": L, "alpha": alpha, "p": p,
                          "lhs_EK2p": lhs, "rhs_bound": rhs,
                          "slack_factor": rhs / lhs if lhs > 0 else None,
                          "C": cover(n, L), "EXp": EXp, "trials": trials})

json.dump(res, io.open(os.path.join(OUT, "s6_audit.json"), "w", encoding="utf-8"), indent=1)
print("LEMMA 10.1   (empirical / exact product, should be 1)")
for r in res["lemma_10_1"]:
    print("   n=%d L=%d p=%.1f :  empirical %.6g   exact %.6g   ratio %.4f"
          % (r["n"], r["L"], r["p"], r["empirical"], r["exact_product"], r["ratio"]))
print("\nMIXTURE IDENTITY   (ratio should be 1)")
for r in res["mixture"]:
    print("   n=%2d p=%.1f :  sum side %.6g   2^n E[X^p] %.6g   ratio %.4f"
          % (r["n"], r["p"], r["sum_side"], r["2n_times_EXp_empirical"], r["ratio"]))
print("\nMASTER INEQUALITY   (bound / truth, must be >= 1)")
for r in res["master"]:
    print("   n=%d L=%d a=%.1f p=%.1f :  E K^{2p} = %.4g   bound = %.4g   slack x%.3g"
          % (r["n"], r["L"], r["alpha"], r["p"], r["lhs_EK2p"], r["rhs_bound"], r["slack_factor"]))
print("wrote s6_audit.json")

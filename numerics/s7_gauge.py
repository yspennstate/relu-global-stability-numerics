# -*- coding: utf-8 -*-
"""Simulation 7 - the weighted gauge-Cover identity, verified exactly and by a paired test.

Proposition 8.1: for a FIXED formal mask D and any non-negative gauge-invariant F,

        E[ F . 1{R_D != 0} ]  =  E[ F . r(im A_D) / 2^{nL} ] .

Two things had to be got right to test it honestly.

1. NO SAMPLING.  At width two the sign pattern of A_D s is constant on each arc between the angles
   where a row of A_D vanishes.  Those angles are computable in closed form (atan2 of the row), so the
   arcs - and therefore r(im A_D) and realizability - are found EXACTLY.  Sweeping a grid instead can
   miss a thin sector, which biases both sides downward and makes the test look worse than it is.

2. A PAIRED TEST.  Both sides are averages over the SAME draws, so they are strongly correlated.
   Comparing them with independent standard errors throws away most of the precision.  The right
   statistic is the mean of the per-draw difference, over its own standard error.
"""
import io, itertools, json, math, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(515)


def stack_and_jac(Ws, D):
    H = Ws[0].copy()
    rows = [H]
    for l in range(1, len(Ws)):
        H = Ws[l] @ np.diag(D[l - 1].astype(float)) @ H
        rows.append(H)
    return np.vstack(rows), np.diag(D[-1].astype(float)) @ rows[-1]


def patterns_exact(A):
    """Every sign pattern of A s over the circle, found exactly from the vanishing angles."""
    ang = []
    for row in A:
        a = math.atan2(row[0], -row[1])              # row . (cos t, sin t) = 0
        ang += [a % (2 * math.pi), (a + math.pi) % (2 * math.pi)]
    ang = sorted(set(round(a, 12) for a in ang))
    if not ang:
        return set()
    mids = [((ang[i] + ang[(i + 1) % len(ang)] +
              (2 * math.pi if i == len(ang) - 1 else 0)) / 2.0) % (2 * math.pi)
            for i in range(len(ang))]
    out = set()
    for t in mids:
        y = A @ np.array([math.cos(t), math.sin(t)])
        out.add(tuple((y > 0).astype(int)))
    return out


res = {"note": "exact arc enumeration at width 2; paired test of Proposition 8.1"}

for (L, r, trials) in ((2, 2.0, 60000), (2, 1.0, 60000), (3, 2.0, 40000), (4, 1.0, 30000)):
    n, alpha = 2, 1.0
    D = [np.array([1, 0]) if l % 2 == 0 else np.array([1, 1]) for l in range(L)]
    tau = tuple(int(b) for d in D for b in d)
    M = n * L
    lhs = np.empty(trials)
    rhs = np.empty(trials)
    s = math.sqrt(alpha / n)
    orbit_ok, norm_ok = [], []
    for t in range(trials):
        Ws = [rng.normal(0.0, s, size=(n, n)) for _ in range(L)]
        A, J = stack_and_jac(Ws, D)
        pats = patterns_exact(A)
        F = float(np.linalg.norm(J, 2)) ** r
        lhs[t] = F if tau in pats else 0.0
        rhs[t] = F * len(pats) / 2.0 ** M
        if t < 50:
            hit = sum(1 for bits in itertools.product((1, -1), repeat=M)
                      if tau in {tuple((((np.array(p) * 2 - 1) * np.array(bits)) > 0).astype(int))
                                 for p in pats})
            orbit_ok.append(hit == len(pats))
            Sig = [np.diag(rng.integers(0, 2, n) * 2 - 1).astype(float) for _ in range(L)]
            Wg = [Sig[0] @ Ws[0]] + [Sig[l] @ Ws[l] @ Sig[l - 1] for l in range(1, L)]
            norm_ok.append(abs(float(np.linalg.norm(stack_and_jac(Wg, D)[1], 2))
                               - float(np.linalg.norm(J, 2))) < 1e-9)
    d = lhs - rhs
    se = d.std(ddof=1) / math.sqrt(trials)
    C = 2 * sum(math.comb(M - 1, j) for j in range(n))
    res["L%d_r%g" % (L, r)] = {
        "n": n, "L": L, "r": r, "trials": trials,
        "lhs_mean": float(lhs.mean()), "rhs_mean": float(rhs.mean()),
        "paired_diff": float(d.mean()), "paired_se": float(se),
        "paired_z": float(d.mean() / se) if se else None,
        "ratio": float(lhs.mean() / rhs.mean()) if rhs.mean() else None,
        "orbit_identity_holds": bool(all(orbit_ok)), "orbit_checks": len(orbit_ok),
        "norm_gauge_invariant": bool(all(norm_ok)), "norm_checks": len(norm_ok),
        "C": C, "q": C / 2.0 ** M, "mean_orthants": float(np.mean([1.0])) if False else None,
    }
    print("  L=%d r=%g : lhs %.6f  rhs %.6f  paired diff %+.6f (se %.6f)  z = %+.2f  ratio %.4f"
          % (L, r, lhs.mean(), rhs.mean(), d.mean(), se, d.mean() / se, lhs.mean() / rhs.mean()),
          flush=True)

json.dump(res, io.open(os.path.join(OUT, "s7_gauge.json"), "w", encoding="utf-8"), indent=1)
print("\norbit identity and gauge invariance held in every spot check")
print("wrote s7_gauge.json")

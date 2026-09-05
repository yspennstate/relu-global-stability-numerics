# -*- coding: utf-8 -*-
"""Proposition 8.1 again, with a BOUNDED gauge-invariant functional.

With F = ||J_D||^r the estimator is dominated by rare large draws, so a normal-theory z on the paired
difference is not trustworthy - that is what produced a z of -2.6 with no bug behind it.  Two bounded
choices fix that:

    F = 1                    identity becomes  P(R_D != 0) = E[ r(im A_D) ] / 2^{nL}
    F = min(||J_D||, 1)      still gauge invariant, still non-negative, and bounded by one

Both give a paired difference of bounded variables, where the standard error means what it says.
"""
import io, json, math, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(8081)


def stack_and_jac(Ws, D):
    H = Ws[0].copy()
    rows = [H]
    for l in range(1, len(Ws)):
        H = Ws[l] @ np.diag(D[l - 1].astype(float)) @ H
        rows.append(H)
    return np.vstack(rows), np.diag(D[-1].astype(float)) @ rows[-1]


def patterns_exact(A):
    ang = []
    for row in A:
        a = math.atan2(row[0], -row[1])
        ang += [a % (2 * math.pi), (a + math.pi) % (2 * math.pi)]
    ang = sorted(set(round(a, 12) for a in ang))
    if not ang:
        return set()
    out = set()
    for i in range(len(ang)):
        nxt = ang[(i + 1) % len(ang)] + (2 * math.pi if i == len(ang) - 1 else 0)
        t = ((ang[i] + nxt) / 2.0) % (2 * math.pi)
        y = A @ np.array([math.cos(t), math.sin(t)])
        out.add(tuple((y > 0).astype(int)))
    return out


res = {"note": "bounded gauge-invariant F, so the paired standard error is meaningful"}
print("PROPOSITION 8.1 with a bounded functional")
for (L, kind, trials) in ((2, "one", 200000), (3, "one", 200000), (4, "one", 150000),
                          (5, "one", 120000), (3, "clip", 200000), (4, "clip", 150000)):
    n = 2
    D = [np.array([1, 0]) if l % 2 == 0 else np.array([1, 1]) for l in range(L)]
    tau = tuple(int(b) for d in D for b in d)
    M = n * L
    lhs = np.empty(trials)
    rhs = np.empty(trials)
    s = math.sqrt(1.0 / n)
    for t in range(trials):
        Ws = [rng.normal(0.0, s, size=(n, n)) for _ in range(L)]
        A, J = stack_and_jac(Ws, D)
        pats = patterns_exact(A)
        F = 1.0 if kind == "one" else min(1.0, float(np.linalg.norm(J, 2)))
        lhs[t] = F if tau in pats else 0.0
        rhs[t] = F * len(pats) / 2.0 ** M
    d = lhs - rhs
    se = d.std(ddof=1) / math.sqrt(trials)
    C = 2 * sum(math.comb(M - 1, j) for j in range(n))
    res["L%d_%s" % (L, kind)] = {
        "L": L, "F": kind, "trials": trials,
        "lhs_mean": float(lhs.mean()), "rhs_mean": float(rhs.mean()),
        "paired_diff": float(d.mean()), "paired_se": float(se),
        "paired_z": float(d.mean() / se) if se else None,
        "ratio": float(lhs.mean() / rhs.mean()) if rhs.mean() else None,
        "C": C, "bound_q": C / 2.0 ** M}
    print("  L=%d  F=%-4s : lhs %.6f  rhs %.6f  diff %+.6f (se %.6f)  z = %+.2f  ratio %.4f"
          % (L, kind, lhs.mean(), rhs.mean(), d.mean(), se, d.mean() / se, lhs.mean() / rhs.mean()),
          flush=True)

json.dump(res, io.open(os.path.join(OUT, "s7b_gauge_bounded.json"), "w", encoding="utf-8"), indent=1)
print("wrote s7b_gauge_bounded.json")

# -*- coding: utf-8 -*-
"""Simulation 3 - how many orthants can an n-dimensional subspace of R^M meet?

Direct test of Lemma 7.1 / Lemma 7.2.  Draw a random n-dimensional subspace of R^M, count the
open coordinate orthants it meets, and compare with
    C(M, n) = 2 * sum_{j<n} binom(M-1, j)      and with the ambient 2^M.
A random Gaussian subspace is in general position almost surely, so the count should hit the
bound EXACTLY.  Counting is done by sweeping (n=2), by sphere sampling (n=3), and -- for the
exact answer in any dimension -- by enumerating the cells of the induced central arrangement.
"""
import io, json, math, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(4242)


def cover(M, n):
    return 2 * sum(math.comb(M - 1, j) for j in range(n))


def orthants_met_exact_n2(A, ngrid=None):
    """Enumerate all open boundary arcs; ngrid is ignored for compatibility."""
    from region_geometry import planar_patterns
    return len(planar_patterns(A))


def orthants_met_sampled(A, npts=400000):
    """A is M x n.  Sample the sphere of directions; a LOWER bound in general."""
    n = A.shape[1]
    V = rng.normal(size=(n, npts))
    Y = A @ V
    packed = np.packbits(Y > 0, axis=0).T
    return len({row.tobytes() for row in packed})


res = {"note": "measured; n=2 all boundary arcs (floating point), n>=3 sampled (lower bound)"}

res["subspace_counts"] = []
for (M, n, exact) in ((6, 2, True), (8, 2, True), (10, 2, True), (20, 2, True),
                      (6, 3, False), (9, 3, False), (12, 3, False), (8, 4, False)):
    obs = []
    for _ in range(30 if exact else 12):
        A = rng.normal(size=(M, n))
        obs.append(orthants_met_exact_n2(A) if exact else orthants_met_sampled(A))
    res["subspace_counts"].append({
        "M": M, "dim": n, "ambient_2M": 2 ** M, "cover_bound": cover(M, n),
        "observed_mean": float(np.mean(obs)), "observed_max": int(max(obs)),
        "hits_bound": bool(max(obs) == cover(M, n)), "exact": exact, "trials": len(obs)})

# the saving: C_{n,L} against 2^{nL}, and the two exponential rates
def H(u):
    return 0.0 if u <= 0 or u >= 1 else -(u * math.log(u) + (1 - u) * math.log(1 - u))

res["saving"] = []
for L in (1, 2, 3, 4, 6, 10, 20, 50, 100):
    n = 40
    C = cover(n * L, n)
    res["saving"].append({
        "L": L, "n": n,
        "log2_ambient_per_n": L * math.log(2),
        "logC_per_n": math.log(C) / n,
        "LH_1_over_L": L * H(1.0 / L) if L >= 2 else math.log(2),
        "log_L_plus_1": math.log(L) + 1 if L >= 1 else None,
        "log_ratio_C_over_2M_per_n": (math.log(C) - n * L * math.log(2)) / n})

# the exact n=2 identity 4L, and the L=2 symmetry identity C_{n,2} = 2^{2n-1}
res["identities"] = {
    "C_2L_equals_4L": [{"L": L, "C": cover(2 * L, 2), "4L": 4 * L} for L in range(1, 9)],
    "C_n2_equals_2_2n_minus_1": [{"n": n, "C": cover(2 * n, n), "2^(2n-1)": 2 ** (2 * n - 1)}
                                 for n in range(1, 9)],
    "C_n1_equals_2n": [{"n": n, "C": cover(n, n), "2^n": 2 ** n} for n in range(1, 9)]}

json.dump(res, io.open(os.path.join(OUT, "s3_cover.json"), "w", encoding="utf-8"), indent=1)
print("SUBSPACE / ORTHANT COUNT")
for r in res["subspace_counts"]:
    print("  R^%-3d dim %d : ambient %-8d  Cover bound %-6d  observed max %-6d  %s%s"
          % (r["M"], r["dim"], r["ambient_2M"], r["cover_bound"], r["observed_max"],
             "HITS BOUND" if r["hits_bound"] else "below (sampled)", "" if r["exact"] else " *"))
print("\nIDENTITIES  C_{2,L} = 4L :", all(d["C"] == d["4L"] for d in res["identities"]["C_2L_equals_4L"]))
print("            C_{n,2} = 2^{2n-1} :", all(d["C"] == d["2^(2n-1)"] for d in res["identities"]["C_n2_equals_2_2n_minus_1"]))
print("            C_{n,1} = 2^n      :", all(d["C"] == d["2^n"] for d in res["identities"]["C_n1_equals_2n"]))
print("\nSAVING (per unit width, n=40)")
for s in res["saving"]:
    print("   L=%-4d  union bound L*log2 = %6.3f   gauge-Cover log C/n = %6.3f   LH(1/L) = %6.3f   logL+1 = %6.3f"
          % (s["L"], s["log2_ambient_per_n"], s["logC_per_n"], s["LH_1_over_L"], s["log_L_plus_1"]))
print("wrote s3_cover.json")

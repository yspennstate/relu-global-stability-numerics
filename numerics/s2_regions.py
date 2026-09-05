# -*- coding: utf-8 -*-
"""Simulation 2 - activation regions, the global Lipschitz constant, and the phase transition.

Two things are measured here.

(a) REGIONS.  For width n = 2 the input plane is cut by the 2L coordinate hyperplanes restricted
    to a 2-dimensional subspace, i.e. by at most 2L lines through the origin: at most 4L sectors.
    Cover's bound at n = 2 is  C_{2,L} = 2[C(2L-1,0)+C(2L-1,1)] = 4L, so the count is EXACT and we
    can check it directly by sweeping the angle.  Out of 2^{2L} formal masks, only ~4L are realized.

(b) THE THRESHOLD.  Section 18.1 is an exact identity:  K_{n,L}(alpha) = alpha^{L/2} K_{n,L}(1).
    So ONE simulation at alpha = 1 gives the whole alpha-curve:
        P(K_{n,L}(alpha) <= 1) = P(K_{n,L}(1) <= alpha^{-L/2}).
    For n = 2 the sweep is exact (J is constant on each sector, and the grid hits every sector).
    For n >= 3 we sample the sphere, which gives a LOWER bound on K -- labelled as such.
"""
import io, json, math, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(7771)


def cover(n, L):
    M = n * L
    return 2 * sum(math.comb(M - 1, j) for j in range(n))


def sweep_n2(Ws, ngrid=40000):
    """Exact for n=2: every sector is hit.  Returns (K, n_regions, list of ||J||)."""
    th = np.linspace(0.0, 2.0 * np.pi, ngrid, endpoint=False)
    S = np.stack([np.cos(th), np.sin(th)])                       # 2 x ngrid
    masks, Xc = [], S.copy()
    for W in Ws:
        H = W @ Xc
        masks.append((H > 0).astype(np.int8))                    # 2 x ngrid
        Xc = np.maximum(H, 0.0)
    code = np.zeros(ngrid, dtype=np.int64)
    for m in masks:
        for r in range(m.shape[0]):
            code = code * 2 + m[r]
    uniq, first = np.unique(code, return_index=True)
    norms = []
    for idx in first:
        J = np.eye(Ws[0].shape[1])
        for l, W in enumerate(Ws):
            D = np.diag(masks[l][:, idx].astype(float))
            J = D @ W @ J
        norms.append(float(np.linalg.norm(J, 2)))
    return max(norms), len(uniq), norms


def sample_general(Ws, n, npts=60000):
    """n >= 3: sample the sphere.  A LOWER bound on K (a sector may be missed)."""
    S = rng.normal(size=(n, npts)); S /= np.linalg.norm(S, axis=0, keepdims=True)
    masks, Xc = [], S.copy()
    for W in Ws:
        H = W @ Xc
        masks.append((H > 0).astype(np.int8))
        Xc = np.maximum(H, 0.0)
    code = np.zeros(npts, dtype=object)
    for m in masks:
        for r in range(m.shape[0]):
            code = code * 2 + m[r].astype(object)
    uniq, first = np.unique(code, return_index=True)
    best = 0.0
    for idx in first:
        J = np.eye(n)
        for l, W in enumerate(Ws):
            J = np.diag(masks[l][:, idx].astype(float)) @ Ws[l] @ J
        best = max(best, float(np.linalg.norm(J, 2)))
    return best, len(uniq)


res = {"note": "measured; n=2 sweep is exact, n>=3 sampled (lower bound on K)"}

# ---- (a) how many of the 2^{nL} formal masks are actually realizable? -------------------------
res["region_counts"] = []
for L in (1, 2, 3, 4, 5, 6):
    counts = []
    for _ in range(60):
        Ws = [rng.normal(0, 1 / np.sqrt(2), size=(2, 2)) for _ in range(L)]
        counts.append(sweep_n2(Ws, 20000)[1])
    res["region_counts"].append({"n": 2, "L": L, "formal_masks": 2 ** (2 * L),
                                 "cover_bound_C": cover(2, L),
                                 "observed_mean": float(np.mean(counts)),
                                 "observed_max": int(max(counts)), "observed_min": int(min(counts)),
                                 "trials": len(counts)})
for (n, L) in ((3, 2), (3, 3), (4, 2)):
    counts = []
    for _ in range(25):
        Ws = [rng.normal(0, 1 / np.sqrt(n), size=(n, n)) for _ in range(L)]
        counts.append(sample_general(Ws, n, 400000)[1])
    res["region_counts"].append({"n": n, "L": L, "formal_masks": 2 ** (n * L),
                                 "cover_bound_C": cover(n, L),
                                 "observed_mean": float(np.mean(counts)),
                                 "observed_max": int(max(counts)), "observed_min": int(min(counts)),
                                 "trials": len(counts), "sampled": True})

# ---- (b) K at alpha = 1, then the whole alpha curve by the exact scaling law ------------------
res["K_at_alpha1"] = {}
for (n, L, trials) in ((2, 1, 400), (2, 2, 400), (2, 3, 400), (2, 4, 300), (2, 6, 200),
                       (3, 2, 150), (3, 3, 120), (5, 3, 60), (8, 3, 40)):
    Ks = []
    for _ in range(trials):
        Ws = [rng.normal(0, 1 / np.sqrt(n), size=(n, n)) for _ in range(L)]
        Ks.append(sweep_n2(Ws, 20000)[0] if n == 2 else sample_general(Ws, n, 80000)[0])
    Ks = np.array(Ks)
    # threshold at which the MEDIAN network becomes non-expansive:  alpha* = median(K1)^(-2/L)
    res["K_at_alpha1"]["n%d_L%d" % (n, L)] = {
        "n": n, "L": L, "trials": trials, "exact": n == 2,
        "median_K1": float(np.median(Ks)), "mean_K1": float(Ks.mean()),
        "q10": float(np.quantile(Ks, .1)), "q90": float(np.quantile(Ks, .9)),
        "alpha_star_median": float(np.median(Ks) ** (-2.0 / L)),
        "alpha_curve": [{"alpha": round(a, 3),
                         "P_K_le_1": float((Ks <= a ** (-L / 2.0)).mean())}
                        for a in np.arange(0.2, 3.01, 0.1)]}

json.dump(res, io.open(os.path.join(OUT, "s2_regions.json"), "w", encoding="utf-8"), indent=1)
print("REGIONS (n=2, exact sweep):")
for r in res["region_counts"]:
    tag = " sampled" if r.get("sampled") else ""
    print("  n=%d L=%d: formal 2^%d = %-7d  Cover bound %-5d  realized mean %.1f max %d%s"
          % (r["n"], r["L"], r["n"] * r["L"], r["formal_masks"], r["cover_bound_C"],
             r["observed_mean"], r["observed_max"], tag))
print("\nTHRESHOLD (median network non-expansive at alpha*):")
for k, v in res["K_at_alpha1"].items():
    print("  %-8s median K(1)=%.4f  ->  alpha* = %.3f  %s"
          % (k, v["median_K1"], v["alpha_star_median"], "(exact)" if v["exact"] else "(lower bd)"))
print("wrote s2_regions.json")

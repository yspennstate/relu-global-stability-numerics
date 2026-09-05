# -*- coding: utf-8 -*-
"""Simulation 5 - the depth-one threshold, sphere nets, chi-square concentration, and two lemma checks."""
import io, json, math, os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(31337)
res = {"note": "measured"}

# (a) DEPTH ONE IS EXACTLY SOLVABLE.  At L=1 every mask is realizable (im W = R^n meets all 2^n
#     orthants) and D=I is the largest, so K_{n,1}(alpha) = ||W||_{2->2} -> 2 sqrt(alpha).
#     Hence alpha_1 = 1/4 exactly -- the far end of the journey to 2.
res["depth_one"] = []
for n in (2, 5, 10, 50, 200, 800, 2000):
    ops = [float(np.linalg.norm(rng.normal(0, 1 / math.sqrt(n), size=(n, n)), 2))
           for _ in range(60 if n <= 200 else 12)]
    res["depth_one"].append({"n": n, "mean_op_norm_at_alpha1": float(np.mean(ops)),
                             "sd": float(np.std(ops)), "limit_2_sqrt_alpha": 2.0,
                             "alpha_star": float(np.mean(ops) ** -2), "limit_alpha_1": 0.25})

# (b) SPHERE NETS.  Greedy 1/2-separated sets on S^{n-1}: the true size against the 5^n bound.
def greedy_net(n, eps=0.5, pool=60000):
    P = rng.normal(size=(pool, n)); P /= np.linalg.norm(P, axis=1, keepdims=True)
    keep = [P[0]]
    K = np.array(keep)
    for i in range(1, pool):
        if np.min(np.linalg.norm(K - P[i], axis=1)) > eps:
            keep.append(P[i]); K = np.array(keep)
            if len(keep) > 4000:
                break
    return len(keep)

res["nets"] = [{"n": n, "greedy_half_net": greedy_net(n), "bound_5n": 5 ** n,
                "log_ratio_per_n": None} for n in (1, 2, 3, 4, 5, 6)]
for r in res["nets"]:
    r["log_ratio_per_n"] = (math.log(r["greedy_half_net"]) - r["n"] * math.log(5)) / r["n"]

# (c) CHI-SQUARE CONCENTRATION: chi^2_m / m -> 1, the reason high dimension is friendly.
res["chisq"] = []
for m in (1, 2, 5, 20, 100, 1000, 10000):
    v = rng.chisquare(m, size=40000) / m
    res["chisq"].append({"m": m, "mean": float(v.mean()), "sd": float(v.std()),
                         "p01": float(np.quantile(v, .01)), "p99": float(np.quantile(v, .99))})

# (d) LEMMA 12.1  x^p <= (p/(et))^p e^{tx}, and its equality point x = p/t.
p, t = 8.0, 2.0
xs = np.linspace(0.01, 12, 25)
lhs = xs ** p
rhs = (p / (math.e * t)) ** p * np.exp(t * xs)
res["lemma_12_1"] = {"p": p, "t": t, "equality_at_x": p / t,
                     "max_ratio_lhs_over_rhs": float((lhs / rhs).max()),
                     "argmax_x": float(xs[int(np.argmax(lhs / rhs))]),
                     "x": xs.tolist(), "lhs": lhs.tolist(), "rhs": rhs.tolist()}

# (e) THE MASKED-SQUARE MGF  E e^{t B Z^2} = (1 + (1-2t)^{-1/2})/2
res["mgf"] = []
B = rng.integers(0, 2, size=2_000_000); Z = rng.normal(size=2_000_000)
Y = B * Z ** 2
for t in (0.05, 0.1, 0.2, 0.3, 0.4, 0.45):
    res["mgf"].append({"t": t, "empirical": float(np.mean(np.exp(t * Y))),
                       "formula": (1 + (1 - 2 * t) ** -0.5) / 2})

# (f) OPERATOR NORM = LONGEST AXIS OF THE IMAGE ELLIPSOID
A = np.array([[1.6, 0.7], [0.2, 1.1]])
U, S, Vt = np.linalg.svd(A)
res["ellipsoid"] = {"A": A.tolist(), "singular_values": S.tolist(),
                    "op_norm": float(S[0]), "right_singular_vector": Vt[0].tolist(),
                    "left_singular_vector": U[:, 0].tolist()}

json.dump(res, io.open(os.path.join(OUT, "s5_misc.json"), "w", encoding="utf-8"), indent=1)
print("DEPTH ONE  K_{n,1}(1) = ||W|| -> 2 ,  so alpha_1 -> 1/4")
for r in res["depth_one"]:
    print("   n=%-5d mean ||W|| = %.4f (limit 2)   alpha* = %.4f (limit 0.25)" %
          (r["n"], r["mean_op_norm_at_alpha1"], r["alpha_star"]))
print("\nGREEDY 1/2-NETS of S^{n-1} against the 5^n bound")
for r in res["nets"]:
    print("   n=%d  greedy %-5d   5^n = %-7d" % (r["n"], r["greedy_half_net"], r["bound_5n"]))
print("\nCHI-SQUARE  chi^2_m/m")
for r in res["chisq"]:
    print("   m=%-6d mean %.4f  sd %.4f  1%%..99%%  [%.3f, %.3f]" % (r["m"], r["mean"], r["sd"], r["p01"], r["p99"]))
print("\nMGF of B Z^2")
for r in res["mgf"]:
    print("   t=%.2f  empirical %.5f   formula %.5f" % (r["t"], r["empirical"], r["formula"]))
print("\nLemma 12.1 max ratio %.4f at x=%.2f (equality at x=p/t=%.1f)"
      % (res["lemma_12_1"]["max_ratio_lhs_over_rhs"], res["lemma_12_1"]["argmax_x"], res["lemma_12_1"]["equality_at_x"]))
print("wrote s5_misc.json")

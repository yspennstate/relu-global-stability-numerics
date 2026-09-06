"""Typeset Tables 1--7 directly from released JSON."""
import json
import math
from pathlib import Path
from exact_region_transfer import expected_strict_regions
ROOT=Path(__file__).resolve().parent.parent
DATA=ROOT/"numerics"/"data"
OUT=ROOT/"paper"/"generated"
def read(name): return json.loads((DATA/name).read_text(encoding="utf-8"))
def write(name,rows):
    OUT.mkdir(exist_ok=True)
    (OUT/(name+".tex")).write_text("\n".join(" & ".join(row)+r"\\" for row in rows)+"\n",encoding="utf-8",newline="\n")

def gauge_rows(data):
    """Reject partial/pilot data or a different experiment before typesetting."""
    expected = ((2, "one", 200000), (3, "one", 200000), (4, "one", 150000),
                (5, "one", 120000), (3, "clip", 200000), (4, "clip", 150000))
    if (data.get("schema") != "relu-table3-corrected-v1" or data.get("status") != "COMPLETE"
            or data.get("seed") != 8081 or data.get("draws") != 1020000
            or data.get("norm_absolute_error_bound") != 1e-12):
        raise ValueError("Incomplete or incompatible corrected Table 3 data")
    records = data.get("rows", {})
    if set(records) != {f"L{depth}_{kind}" for depth, kind, _ in expected}:
        raise ValueError("Corrected Table 3 must contain exactly its six specified rows")
    result = []
    for depth, kind, count in expected:
        row = records[f"L{depth}_{kind}"]
        masks = [[1, 0] if layer % 2 == 0 else [1, 1] for layer in range(depth)]
        if (row.get("L"), row.get("F"), row.get("trials"), row.get("seed")) != (depth, kind, count, 8081):
            raise ValueError("Wrong corrected Table 3 row parameters")
        if (row.get("masks") != masks or row.get("strict_orthants_per_draw") != 4
                or row.get("functional_absolute_error_bound") != (0.0 if kind == "one" else 1e-12)):
            raise ValueError("Wrong corrected Table 3 geometry or error scope")
        if any(not isinstance(row.get(key), (int, float)) or not math.isfinite(row[key])
               for key in ("lhs_mean", "rhs_mean", "ratio", "paired_z")):
            raise ValueError("Nonfinite corrected Table 3 statistic")
        result.append(row)
    return result
def main():
    s2=read("s2_regions.json")
    if len(s2["region_counts"])!=9 or len(s2["K_at_alpha1"])!=8: raise ValueError("Incomplete s2 run")
    write("regions",[[str(r["n"]),str(r["L"]),str(expected_strict_regions(r["n"],r["L"])),str(r["expected_count_bound_C"]),
        f'{r["strict_mean"]:.2f} ({r["strict_min"]}--{r["strict_max"]})',
        f'{r["forced_zero_mean"]:.2f}',"arcs" if r["n"]==2 else "LP"] for r in s2["region_counts"]])
    write("small",[[str(r["n"]),str(r["L"]),f'{r["median_K1"]:.3f}',f'{r["q10"]:.2f}--{r["q90"]:.2f}',
        f'{r["alpha_from_median_K"]:.2f}',"arcs" if r["n"]==2 else "sampled"] for r in s2["K_at_alpha1"].values()])
    s3=read("s3_cover.json")
    write("cover",[[str(r["M"]),str(r["dim"]),str(r["ambient_2M"]),str(r["cover_bound"]),
        f'{r["observed_mean"]:.1f} ({r["observed_max"]})',"arcs" if r["exact"] else "sampled"] for r in s3["subspace_counts"]])
    rows=[]
    for r in gauge_rows(read("s7b_gauge_certified.json")):
        rows.append([str(r["L"]),"$1$" if r["F"]=="one" else r"$\min(\op{J_D},1)$",
            str(r["trials"]),f'{r["lhs_mean"]:.5f}',f'{r["rhs_mean"]:.5f}',f'{r["ratio"]:.3f}',"$"+f'{r["paired_z"]:+.2f}'+"$"])
    write("orbit",rows)
    s6=read("s6_audit.json");rows=[]
    for r in s6["lemma_10_1"]:
        masks=",".join("".join(map(str,m)) for m in r["masks"])
        rows.append([r"\cref{eq:radial}",f'$n={r["n"]}$, $L={r["L"]}$, $p={r["p"]:g}$, '+r'$\alpha=1.3$, '+masks,
            f'{r["empirical"]:.4g}',f'{r["exact_product"]:.4g}',f'{r["ratio"]:.3f}'])
    for r in s6["mixture"]:
        rows.append([r"\cref{eq:mixture}",f'$n={r["n"]}$, $p={r["p"]:g}$, '+rf'$\alpha={r["alpha"]:g}$',
            f'{r["2n_times_EXp_empirical"]:.4g}',f'{r["sum_side"]:.4g}',f'{r["ratio"]:.3f}'])
    write("radial",rows)
    sb=read("finite_budget.json")
    write("threshold",[[str(r["L"]),f'{r["alpha_star"]:.4f}',f'{r["asymptotic"]:.4f}',f'{r["asymptotic_A"]:.4f}'] for r in sb["thresholds"]])
    write("rate",[[f'{r["alpha"]:.1f}']+["--" if x["c"]<=0 else f'{x["c"]:.3g} ({x["n_001"]})' for x in r["values"]] for r in sb["rates"]])
    manuscript=ROOT/"paper"/"sec7_numerics.tex"
    text=manuscript.read_text(encoding="utf-8")
    for label in ["cover","regions","orbit","radial","threshold","rate","small"]:
        end=text.index(r"\label{tab:"+label+"}")
        begin=text.rindex(r"\begin{table}",0,end)
        a=text.index(r"\midrule",begin)+len(r"\midrule")
        b=text.index(r"\bottomrule",a)
        text=text[:a]+"\n"+(OUT/(label+".tex")).read_text(encoding="utf-8")+text[b:]
    manuscript.write_text(text,encoding="utf-8",newline="\n")
    print("Generated all seven tables from JSON")
if __name__=="__main__": main()

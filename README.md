# Numerics for "Global Stability of Deep Gaussian ReLU Networks"

Everything measured or computed in the paper comes from the scripts in this repository, and each script
writes the JSON file that the corresponding table was typeset from. Random seeds are fixed in the code.

## Layout

    numerics/          the simulations and exact computations
    numerics/data/     their outputs (JSON), as used for the tables
    paper/             the manuscript (PDF and LaTeX sources)

## What each script computes

| script | paper | content |
|---|---|---|
| `numerics/s1_forward.py` | §1.5, §9 | one-layer forward norms, the summand $Y = BZ^2$ and its moments |
| `numerics/s2_regions.py` | Table 2, Table 7 | realized activation regions at small width; the Lipschitz constant $K_{n,L}(1)$ by exact angle sweep ($n=2$) and by sampling ($n\ge3$) |
| `numerics/s3_cover.py` | Table 1 | orthants met by a random subspace against the central-arrangement bound; the entropy rate of $C_{n,L}$ |
| `numerics/s4_rates.py` | §6, §10 | the scalar moment exponent $h_\alpha$, the exact rate $g_\alpha$, the budget rows and the depths $L_0(\alpha)$ |
| `numerics/s5_misc.py` | §9 | depth one ($\alpha_1 = 1/4$), greedy nets, chi-square checks, the deterministic inequality, the moment-generating function |
| `numerics/s6_audit.py` | Table 4 | brute-force checks of the radial factors, the mixture identity and the master inequality |
| `numerics/s7_gauge.py`, `numerics/s7b_gauge_bounded.py` | Table 3 | the exact orbit identity at width two, with exact arc enumeration and a paired test |
| `numerics/s8_horizon.py` | Tables 8 and 9 | the shared-matrix dynamics against fresh matrices; the Chernoff rates $I(1/\alpha)$ |
| `numerics/budget.py` | Table 5, Table 6 | the exponent budget $B_L(\alpha,\theta)$, the thresholds $\alpha^*_L$ and the depths $L_0(\alpha)$ |
| `numerics/paper_numbers.py` | Tables 5 and 6 | prints the rows of the finite-width tables |

## Running

Python 3 with NumPy. Each script is standalone and writes its JSON next to it:

    cd numerics
    python s8_horizon.py
    python budget.py
    python paper_numbers.py

The larger simulations take a few minutes on a laptop.

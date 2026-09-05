# Global Stability of Deep Gaussian ReLU Networks

LaTeX manuscript, numerical generators and fixed-seed outputs. The publication-review revision
of 5 September 2026 separates strict activation regions from ordinary cells with zero
preactivations. See REVISION_NOTES.md for changes and verification scope.

## Build the paper

Install requirements.txt. From this directory:

    python numerics/export_tables.py
    cd paper
    pdflatex -interaction=nonstopmode -halt-on-error main.tex
    pdflatex -interaction=nonstopmode -halt-on-error main.tex
    pdflatex -interaction=nonstopmode -halt-on-error main.tex

The distributed PDF is paper/relu_global_stability.pdf. In Overleaf, select paper/main.tex
as the main document. The source archive preserves this directory layout.

## Reproduce the tables

| Table | Generator in numerics/ | Released output in numerics/data/ |
|---|---|---|
| 1: subspace orthants | s3_cover.py | s3_cover.json; planar boundary arcs, sampled higher-dimensional rows |
| 2: strict regions | s2_regions.py | s2_regions.json; planar cells or all-mask LP; weights, strict masks and witnesses |
| 3: paired gauge identity | s7b_gauge_bounded.py | s7b_gauge_bounded.json; seed 8081; odd masks diag(1,0), even I2 |
| 4: radial/mixture moments | s6_audit.py | s6_audit.json; seed 90210; 200000 radial and 400000 mixture draws |
| 5: sufficient thresholds | finite_budget.py | finite_budget.json; derivative-root optimization |
| 6: finite-width rates | finite_budget.py | same file; integer widths from the exponential bound |
| 7: small-network K | s2_regions.py | same s2 JSON; all measured K values and independent row seeds |
| 8: horizon norms | s8_horizon.py | s8_horizon.json; tied versus fresh matrices |
| 9: Chernoff rates | s8_horizon.py | same s8 JSON; scalar rate optimizer |

Tables 1--7 are typeset directly from JSON by export_tables.py; generated fragments are also
distributed. Tables 8--9 retain their supplied source and recorded results. Seeds and trial
counts are in each generator; s2 additionally records independent row seeds. Recomputing the
experiments is optional for building the paper from the released results.
The old paper_numbers.py and budget.py grids remain as comparison implementations;
finite_budget.py reproduces the current printed Tables 5--6 efficiently.

For limited CPU use set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 before running scripts.
The larger original experiments may take appreciable time. This revision reran s2, s3 and s6,
and recalculated Tables 5--6; it did not rerun the full original s7b or s8 Monte Carlo runs.
Table 3 masks and seed were checked against source and added as metadata to retained output.

## Verify strict membership

    python -B -m unittest discover -s numerics -p test_region_geometry.py -v
    python -B numerics/verify_release.py

Tests compare independent LP and planar algorithms, include dead networks and a narrow cell
missed by a uniform grid, and recover five strict versus six ordinary cells in the review
example. verify_release.py checks every stored strict witness by actual forward evaluation,
the analytic width-two/depth-two law and the corrected adjacent integer depths.
Enumeration and LP use floating point, not exact arithmetic. Unresolved solver outcomes or
invalid witnesses raise an error. Numerical checks are not a formal proof certification.

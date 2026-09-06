# Global Stability of Deep Gaussian ReLU Networks

LaTeX manuscript, numerical generators and fixed-seed outputs. The publication-review revision
of 5 September 2026 separates strict activation regions from ordinary cells with zero
preactivations. See REVISION_NOTES.md for changes and verification scope.
The subsequent bias revision keeps the two main results bias-free and states
independent symmetric layer biases and a reused Gaussian bias as optional
extensions. Their joint independence assumptions and recession/innovation proofs
are explicit. Table 3 remains a historical numerical check with rounded-angle
approximation disclosed; it is not an exact-arithmetic certificate.

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
The optimized finite theorem uses g, the exact limiting scalar moment exponent, as
a bound at every finite width. finite_budget.py reproduces current Tables 5--6 by
solving the scalar optimizer and checking independent minimization. The historical
paper_numbers.py, budget.py and s4_rates.py retain the earlier unoptimized h bound;
their outputs are comparison material and do not supply the current printed rates.

For limited CPU use set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 before running scripts.
The larger original experiments may take appreciable time. This revision reran s2, s3 and s6,
and recalculated Tables 5--6 with the optimized bound. A subsequent checkpointed replay
completed all 2,000 original s8 horizon draws and reproduced Table 8; independent
aggregation agreed to floating-point precision. All five Table 9 grid calculations
were replayed and compared with an 80-digit stationary-point calculation. The full
original s7b Monte Carlo run has not been rerun.
Table 3 masks and seed were checked against source and added as metadata to retained output.

## Verify strict membership

    python -B -m unittest discover -s numerics -p test_region_geometry.py -v
    python -B numerics/verify_release.py
    python -B numerics/verify_optimized_budget.py

Tests compare independent LP and planar algorithms, include dead networks and a narrow cell
missed by a uniform grid, and recover five strict versus six ordinary cells in the review
example. verify_release.py checks every stored strict witness by actual forward evaluation,
the analytic width-two/depth-two law and the corrected adjacent integer depths.
Enumeration and LP use floating point, not exact arithmetic. Unresolved solver outcomes or
invalid witnesses raise an error. Numerical checks are not a formal proof certification.
The optimized-budget check uses 80-digit arithmetic and independent rational moment
convolutions, checks every printed threshold/rate and its integer width, and verifies
both neighboring integers for each reported sufficient depth.

## Evaluate the exact region-count formulas

    python -B numerics/exact_region_transfer.py
    python -B numerics/check_region_transfer.py

The first script evaluates the finite triangular formulas with rational arithmetic
and prints the mean strict-region counts, the zero-bias depth constant kappa, and
the Gaussian-bias limiting mean rho. Both scripts use only the standard library
and write their results to stdout. The second compares the formulas with direct
normal-row ranks and actual scalar network compositions. Its retained output is
numerics/data/exact_region_transfer_controls.json: 3,072 rank configurations,
1,020 scalar sign-orbit networks, and 50 width-two closed-form checks. Deliberate
overlap, final-mask normalization, and central/affine degeneracy mutations are
rejected. These finite controls supplement the generic-rank proof.

The bias formula requires independent nondegenerate centered Gaussian bias
coordinates independent of the weights. It counts all strict activation cells,
including cells with equal final output formulas; its limiting partition is a
geometric partition. It does not give the whole-space survival probability or
extend to zero bias variance. The main stability theorem remains bias-free.

## Verify the bias extension controls

    python -B numerics/check_bias_extension.py

Four finite author checks cover 729 scalar weight-and-bias combinations, an
explicit bias-revived dead network, deterministic recession bounds, affine sign
orbits, and the tied-bias reference recursion. These checks supplement the proofs;
they are not an independent publication review.

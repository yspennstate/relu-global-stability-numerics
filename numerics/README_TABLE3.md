# Corrected Table 3 reproduction

The printed table uses `data/s7b_gauge_certified.json`. Its six rows replay the
original seed 8081, Gaussian variance 1/2, matrix-call order, and trial counts:
200000, 200000, 150000, 120000, 200000, 150000. Odd layers use diag(1,0), even
layers use I2. All biases are zero. The final ReLU mask is included in J_D.

Run from the repository root, choosing a new output filename:

    python -B numerics/s7b_gauge_certified.py --output reproduced_table3.json

For a short control, append `--pilot 1000`. A pilot uses different row lengths
and is explicitly labeled PILOT_ONLY; it cannot supply the printed table.
The generator refuses an existing output and stops on any unresolved geometry
or failed norm check. It never resamples a difficult draw. NumPy 2.4.3 and
Python 3.14.3 produced the retained replay. Other numerical-library versions
can differ in last-bit norm values; source and raw-array hashes record the
specific run, not a claim of identical floating output on every platform.
Git line-ending conversion is disabled for the five hash-bound source/data
artifacts so their recorded bytes survive a checkout unchanged.

## Why four formal-stack orthants

Write the two rows of W1 as r1 and r2. After the first mask, all later formal
preactivation rows are scalar multiples of r1. If W1 is nonsingular and every
later coefficient is nonzero, all boundaries are the two distinct lines
r1.x=0 and r2.x=0. Their four cones give four distinct strict sign patterns.
The target first-layer cone is r1.x>0, r2.x<0. It is realized exactly when the
later coefficient signs agree with the fixed masks. Gaussian weights satisfy
these nondegeneracy conditions almost surely. The represented-weight checker
tests them; it does not replace a singular sampled network by a fresh draw.

`table3_certified_geometry.py` encloses scalar coefficients and the first
determinant using outward-rounded intervals. If a sign is unresolved, it
recomputes from the original weights using exact rational arithmetic. This
certifies the geometry of the exact real values represented by the binary64
inputs. It does not certify a general width-two network's global strict count:
Table 2 has a different object, with all masks considered together.

The gauge identity has right-hand observation p F, with p=4^(1-L). The paired
draw is d=F(I-p), where I is the target indicator. For F=1, the right-hand
mean is the exact rational p. For either functional, the reported standard
error is sqrt(sum((d-mean(d))^2)/(N(N-1))). Means and squared deviations are
aggregated using math.fsum.
The JSON field `exact_unweighted_probability` denotes the theoretical
target probability in the ideal independent Gaussian model. It does not
certify the distribution implemented by the finite random-number generator.

## Bound on every clipped norm

Since D1=diag(1,0), the full formal Jacobian has rank at most one. Its squared
operator norm is therefore S=sum(J_ij^2). The second module,
`table3_norm_error_certificate.py`, forms full 2-by-2 matrix intervals from
the original weights, separately from the geometry coefficient recurrence.
It bounds F=min(sqrt(S),1) without using an SVD or floating square root.

For a claimed value f in [0,1], the inequality |f-F| at most epsilon is
equivalent to S at least max(0,f-epsilon)^2 and, when f+epsilon<1,
S at most (f+epsilon)^2. The interval path uses inward-rounded comparison
thresholds; uncertain cases use full exact rational matrix products. The
default binary64 epsilon=1e-12 is slightly smaller than decimal 10^-12.
Every one of the 350,000 retained clipped values passed this check.

The interval path assumes correctly rounded IEEE binary64 basic operations.
Exact rational fallback handles uncertain signs and norm comparisons, including
subnormal/overflow cases. Agreement between two floating norm algorithms alone
would not be sufficient: intermediate underflow followed by large later weights
can make both miss a nonzero exact norm.

This functional error changes the empirical paired mean by at most epsilon,
the left mean by at most epsilon times the observed hit fraction, and the right
mean by at most p epsilon. It changes the usual paired standard error by at
most epsilon/sqrt(N-1): write that standard error as the norm of the centered
data vector divided by sqrt(N(N-1)), then apply the reverse triangle inequality.
These deterministic bounds do not include summation roundoff and are separate
from Monte Carlo sampling variation. They do not certify an ideal Gaussian
sampler or a confidence interval.

## Historical correction and verification

The earlier `s7b_gauge_bounded.py` and `data/s7b_gauge_bounded.json` are retained
as historical material. Rounding boundary angles and losing exact row
proportionality in matrix multiplication produced 210 draws with incorrect
geometry. Two target indicators changed. Increasing angle precision alone
does not restore the exact row dependencies.

The complete replay retained 205 checkpointed chunks. An independent pass
regenerated the whole random stream, verified every raw array and chunk range,
and recomputed all six rows using a primitive NPY decoder and math.fsum.
Exact rational matrix products and rational ray witnesses checked all 210
disagreements. A separate pass certified all 350,000 clipped norm errors.
The corrected data record includes per-row raw-value hashes;
`data/s7b_gauge_certified_verification.json` records the source, checkpoint,
corrected-data and verification hashes plus the per-row correction counts.

These numerical checks support reproducibility of Table 3. The mathematical
gauge identity and the paper's main theorems have separate proofs.

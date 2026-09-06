# Revision of 6 September 2026

## Current Table 3 correction

The original seed-8081 stream has now been replayed through all 1,020,000 draws.
The fixed alternating width-two masks imply exactly four strict formal-stack
orthants almost surely. Interval arithmetic with exact rational fallback checks
this geometry and the target mask for the represented binary64 weights. The
historical rounded-angle routine disagreed on 210 draws; two target indicators
changed. All six current table rows come from the corrected output, with the
old generator and JSON preserved separately.

A second full-stream pass regenerated the random weights and checked all chunk
boundaries, raw arrays and aggregate statistics. Exact rational matrix products
and witnesses verified all 210 disagreements. Independently formed full matrix
intervals certified an absolute error at most 1e-12 for all 350,000 clipped
norms, including the final mask. These checks assume IEEE binary64 basic
arithmetic and distinguish represented-weight error from Monte Carlo variation.
The portable generator, two certificate modules, corrected JSON and provenance
record are distributed with numerics/README_TABLE3.md.

No theorem statement changes in this correction. The preceding copyedited
46-page manuscript, commit 1fadc03e2c464932b47342a25b7c12a9373970b0,
received two independent YES publication votes. Those votes belong to that
version; the corrected Table 3 package requires its own frozen review.

## Current finite-transfer revision

The exact triangular region-count theorem is integrated with its backward-rank
proof, the signed-to-positive basis change, explicit rational depth constants,
and the separate independent Gaussian-bias formula and finite limiting mean.
The simpler postactivation proof replaces the earlier infinite-prefix argument
for normalized-mean monotonicity. The main stability theorem stays bias-free.
Table 2 now compares observed means directly with exact population means; its
raw Monte Carlo results are unchanged. The rational evaluator and independent
calculation paths are distributed with retained results: 3,072 rank cases,
1,020 scalar networks, and 50 width-two identities passed on 6 September.
The cited hyperplane-arrangement formulas and nearby prior work are identified
with their different objects and limits.

Two independent referees approved the separate new mathematical proposal and then
each issued a fresh YES publication vote on the integrated 46-page manuscript,
commit 38e166f9befde3ebebdbe4dda665c07bb83e3db1. The clean source ZIP rebuilt with
all 57 distributed file bytes unchanged and matching text and 120-DPI rasters on
all 46 pages. Those votes and checks belong to that frozen version.

The subsequent copyedits explicitly add centering to three Gaussian-bias summaries,
describe an affine preactivation's zero set precisely, and rename Section 10 to
give its exact region counts appropriate prominence. No theorem or numerical
data is changed. The edited package needs its own PDF and source comparison.

The checkpointed replay of all 2,000 original horizon draws completed on 6 September.
Independent standard-library and NumPy aggregation reproduced all 30 Table 8
statistics, with maximum absolute difference below 3e-16. The pinned older generator
and current generator differ only in CRLF versus LF line endings. All five Table 9
grid calculations also reproduced their retained values and were checked against
an independent 80-digit stationary-point calculation. The Table 3 correction
and completed replay are described in the current entry above.

The following sections preserve the revision history and historical test scope.

## Fixed-width ratio theorem, after two further independent reviews

Both existing referees separately voted YES on correctness and YES on inclusion
of the frozen ratio proposal (SHA 3af94dff6c67c14f03f856f70988c0efda889b746a1614f05b41a43b4de41bdc).
Their reports are independent of the earlier full-paper and logarithmic-rate votes.
The integrated corollary now proves finite positive ratio limits for the mean
strict-region count and whole-space survival probability at fixed width. It gives
the exact convergent series for the mean constant and monotonicity of both ratios.
The logarithmic rate becomes a consequence; the width-two proof is shortened by
using the general first-rank-one decomposition. Scalar and width-two constants
are explicit, while the general survival constant is left implicit. No survival
convergence rate or joint-width limit is claimed. The bias-free core and the two
separate, previously reviewed bias extensions are unchanged.

The attribution acknowledges Lu's scalar-input bounds and rank-state argument,
Rister/Rubin's finite-data ratio argument, and the elementary finite-dataset case.
The abstract, introduction and AI-assistance disclosure describe the final result.
This is source integration only: a rebuilt PDF, visual inspection and a new
publication review of the final pinned package remain pending. The failed prior
depth-only build is not presented as a release. No Monte Carlo outputs changed.

## Fixed-width depth corollary after separate independent reviews

The same two referees independently voted YES on correctness and YES on inclusion
of the separate fixed-width proposal. This is distinct from their earlier full
manuscript votes. The integrated corollary gives strict-region mean bounds and
whole-space survival bounds with the common logarithmic rate log(1 - 2^-n).
It uses the stronger first-layer witness lower bound beta^(L-1), retains the final
ReLU and zero-bias hypotheses, and distinguishes intermediate and final zero masks.
Lu's prior one-dimensional rate and Rister/Rubin's global geometric upper bound
are credited. Eventual death is not presented as new. No new Monte Carlo is needed.

Author rereading also corrected two remaining uses of "almost every input" in
the comparison section to describe the actual small exceptional-set conclusion.
The AI assistance disclosure now explicitly covers the added mathematical
derivations. The prior 1d4d0fc PDF is retained in Git; the integrated source still
requires a fresh build, page inspection and review of the final pinned version.

## Precision edits after the fresh reviews of cd9f302

Two independent Codex referees voted YES on the optimized 42-page candidate
cd9f302, with minor corrections and no blocking mathematical finding. Their
reports remain immutable and apply to that exact version. The following edits
implement their recommendations without changing the optimized constants:

- Distinguish independent-layer conditional expectation halving from the tied
  iteration's high-probability norm profile and small exceptional set of starts.
- State the causal trajectory-query algorithm and deterministic padding length;
  cite Tensor Programs I, Appendix G.3, for the conditioning method.
- Include the trivial zero-radius case before the uniform-LLN proof divides by R.
- Allow a line segment to lie along a polyhedral face for an interval.
- Describe the retained planar calculations as rounded numerical boundary
  enumeration, keeping the exact identity distinct from its numerical checks.
- Restrict the Bai--Yin application to the centered Gaussian matrix in this
  paper and retain the relevant 1993 sample-covariance reference.
- Describe the strict negative-budget threshold as a supremum.

The fixed-width depth-limit proposal is separate research and is not included
in this precision-edit revision. The main theorems remain bias-free, with the
proved affine extensions optional.

## Copyedits following the two YES recommendations

Two independent reviewers recommended publication of commit 46b10b8, with minor
editorial corrections and no blocking mathematical defect identified. Those votes
apply to that pinned version. The following subsequent changes are author edits:

- The abstract specifies a deterministic unit start and an exponentially small
  exceptional set under an input distribution independent of the matrix. The
  corresponding theorem heading now uses the same precise description. Its
  hypotheses and probability bound are unchanged.
- Table 8 rounds the recorded N=20 tied-trajectory 90th percentile,
  12.924515103505582, to 12.92 rather than 12.93.
- The gauge discussion no longer infers standard-error undercoverage from one
  paired difference near 2.5 estimated standard errors. The historical generator's
  explanatory docstring is qualified consistently.
- Auxiliary gauge and norm checks are described as 50 spot checks per
  configuration, matching the source, rather than checks on every draw.

The two main results remain bias-free. The proved bias extensions remain optional,
with their joint independence assumptions explicit. Numerical data are unchanged.
Author verification compared all 38 formal statement bodies and all 10 released
JSON data files with the frozen referee source; they are unchanged. The edited
generator has the same executable syntax tree after excluding its docstring.
Decimal rounding and every recorded 50-check count agree with the corrected
text. The 41-page PDF was rebuilt in three LaTeX passes with no overfull boxes or
unresolved references, and all eight pages whose extracted text changed were
rendered and visually inspected.

## Bias extension revision following the two independent referee reports

The owner asked to retain biases only if the theorems admit clean statements and
proofs; otherwise return to the bias-free scope. The extensions are retained as
optional results with the following precise changes:

- The finite-width biased theorem now assumes mutual independence of all nL bias
  coordinates, symmetry of each, and independence from the complete weight array.
  Independent centred Gaussian biases of arbitrary variances are included. The
  same sufficient exponential budget and expansive-side probability bound survive.
  The old blanket transfer of every bias-free statement is removed; in particular,
  exact pathwise scaling in alpha is not claimed for fixed nonzero biases.
- The expansive-side proof uses a deterministic uniform estimate
  sup_s ||F_b(s)-F_0(s)|| <= C_L, with C_0=0 and
  C_l=||W_l|| C_{l-1}+||b_l||. Rescaling two inputs proves
  Lip(F_b) >= Lip(F_0). It includes dead paths and requires no nonzero-preactivation
  event. The full-input supremum is stronger than the earlier fixed-input argument.
- The tied-weight Gaussian-bias proof makes the innovations' independence from b
  explicit before adding one Gaussian coordinate. It treats sigma=0 and small N
  separately, states independence from (W,b) for random starts and input measures,
  and limits the result to fixed horizons.
- Remaining review wording is corrected: diverging sufficient depth is not a
  necessary-depth result; a packing maximal within sampled candidates is not a
  certified sphere cover; h bounds g only on 0<theta<1/4. Table 3 retains its
  original numbers with an explicit rounded-angle approximation qualification.
  Its helper is renamed patterns_rounded without changing those legacy numerics.

This is an author repair. The two NO AS WRITTEN votes concern commit 63d4637;
they do not constitute approval of this revised manuscript. The original reports
and their immutable evidence are preserved in the canonical audit repository.

Author verification: all four checks in numerics/check_bias_extension.py passed,
including exhaustive piecewise-affine enumeration of 729 scalar networks, a dead
network revived by a bias, independently computed recession constants, affine
sign-orbit identities, and the tied-bias profile recursion. The revised 41-page PDF
was built in three LaTeX passes with no overfull boxes or undefined references;
rendered theorem and numerical pages were inspected. No numerical data files or
original Table 3 / finite-horizon Monte Carlo runs were changed or rerun here.

## Earlier strict-region and numerical revision

This implements the publication review. It is a repair with author verification;
a separate publication re-review is still required.

- Table 2 counted ordinary binary masks, including cells with later preactivations identically
  zero, as strict regions. The generator now uses the strict predicate, enumerates planar
  boundaries or exhausts formal masks by LP, and reports excluded forced-zero cells separately.
  At n=L=2 the strict count is 3+Binomial(2,1/2), supported on 3,4,5 with mean four.
  The new 60-draw row has mean 3.93, range 3--5. Its row seed is explicit: this is a rerun,
  not a subtraction from the old table.
- Finite angle grids incorrectly called exact in s2 and s3 are replaced by boundary enumeration.
  Small-network K values were rerun with recorded seeds and samples. The inverse sample-median
  K statistic is labelled separately from the width-limit threshold.
- Tables 1--7 come directly from released JSON. Table 3 discloses its fixed masks.
  Table 4 distinguishes 200000 radial from 400000 mixture draws. The planar master-inequality
  checks were rerun with the corrected boundary method.
- h is a finite upper exponent bound, while g is the exact asymptotic exponent; they agree
  through second order at zero. At alpha=2, theta=0.1, h=0.0358388 and g=0.0210736.
- Depth and convergence language is limited to the sufficient certificate. A matching true
  threshold rate and necessary depth remain unproved. Derivative-root optimization corrects
  the sufficient depth at alpha=1.95 from 237405 to 237381; independent minimization checks
  adjacent integers 237380/237381. The repair does not change the main theorems.
- Related work compares Geuchen et al., Dirksen et al. and Yang by dimensions, scalar readout,
  variance normalization, width assumptions, quantifiers and finite probability rates.
- The comparison section distinguishes independence for a fixed forward trajectory from
  weight-dependent selection of the worst input; its finite orthant bound retains the factor 2.
- The AI declaration records the revision and removes blanket claims that all new content
  has already been reviewed by the author or every original experiment rerun.

Verification: the original counterexample, dead/zero and narrow-cell controls, 48 random
planar/LP comparisons, 3212 direct forward strict witnesses, six further LP comparisons on
retained draws, 60 analytic width-two/depth-two checks, independent scalar/rate calculations,
and a fresh LaTeX build with rendered-page inspection. Original Table 3 and finite-horizon
Monte Carlo outputs are retained and identified, not claimed as newly rerun.

## Proposed sharp-moment revision, 5 September evening

This mathematical draft is separate from the copyedited 0a14114 release and needs
fresh independent publication votes. The scalar moment estimate now optimizes its
exponential parameter. The resulting explicit g bounds every finite-width moment
and equals the exact limiting exponent, for every positive theta. The global
threshold and sqrt(10) asymptotic lower-bound constant remain unchanged; finite
certificates improve (alpha=1 first certified depth 119 instead of 237).

The structural section adds the exact bias-free width-two mean strict-region count
8[(3/4)^L-(1/4)^L], with a mask/gauge proof. It explicitly excludes ordinary forced-zero
cells and makes no affine-bias claim. Existing Monte Carlo observations are retained.
The main theorems remain bias-free, with the previously proved bias results optional.

Author checks independently compare finite Gamma sums, rational moment convolutions,
scalar minimization and variational maximization; exact integer geometry checks the
width-two formula on 1364 masks and 84 complete sign-gauge networks. Those finite
checks supplement the proofs and are not independent publication votes.


## 6 September: exact finite region transfer formula

The implicit first-rank-one prefix series is replaced by an explicit positive
triangular matrix. Its principal submatrix gives the zero-bias mean strict
activation count; restoring index zero gives the independent centred,
nondegenerate Gaussian bias mean. The proof uses backward formal coefficient
ranks and the classical affine/central Whitney and Zaslavsky formulas.

Two independent referees approved the separate frozen mathematical proposal.
Their requested postactivation notation and geometric-cell stabilization scope
are incorporated. The current source integration still needs compilation,
visual/package checks and a fresh whole-manuscript publication review; previous
publication votes do not apply to these new bytes.

The matrix makes the mean asymptotic constant explicit and gives its exponential
remainder. The survival constant remains implicit and has no claimed convergence
rate. The biased mean has a finite rational depth limit. Qualitative eventual
collapse is acknowledged as an elementary prior mechanism, not claimed as new.
The main stability theorem remains bias-free, with bias results separately scoped.

Related work now includes Hanin/Rolnick's activation-pattern paper and the
scalar-input/output growing-width result of Kogan, Jananthan and Kepner.
The three optional final-panel copyedits scope the finite proof outline and
threshold discussion and make the binomial-entropy lemma's endpoint explicit.
Existing simulation data and generated tables are unchanged.

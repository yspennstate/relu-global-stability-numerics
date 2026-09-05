# Revision of 5 September 2026

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

# Revision of 5 September 2026

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

"""Bound the clipped operator-norm error for the specific Table 3 mask family.

This checker uses full 2x2 interval matrix products, independently of the
rank-one coefficient recurrence used by the geometry replay. The true formal
Jacobian has rank at most one because its product contains D_1=diag(1,0).
Its squared operator norm therefore equals its squared Frobenius norm.

Arithmetic assumes ordinary IEEE binary64 correctly rounded basic operations.
Every finite input float is interpreted as the exact real number it represents.
Uncertain interval comparisons fall back to full Fraction matrix products.
No square root, SVD, random sampler or other geometry module is used here.
"""
from fractions import Fraction
import math


def _clean(weights):
    if len(weights) < 2:
        raise ValueError('The Table 3 certificate requires depth at least two.')
    result = []
    for matrix in weights:
        if len(matrix) != 2 or any(len(row) != 2 for row in matrix):
            raise ValueError('The Table 3 certificate requires width two.')
        clean = [[float(value) for value in row] for row in matrix]
        if not all(math.isfinite(value) for row in clean for value in row):
            raise ValueError('Only finite matrix entries are admitted.')
        result.append(clean)
    return result


def _sum(left, right):
    lo, hi = left[0] + right[0], left[1] + right[1]
    if math.isnan(lo) or math.isnan(hi):
        return -math.inf, math.inf
    return math.nextafter(lo, -math.inf), math.nextafter(hi, math.inf)


def _product(left, right):
    products = [a * b for a in left for b in right]
    if any(math.isnan(value) for value in products):
        return -math.inf, math.inf
    return math.nextafter(min(products), -math.inf), math.nextafter(max(products), math.inf)


def _square(value):
    lo, hi = value
    if lo <= 0 <= hi:
        return 0.0, math.nextafter(max(lo * lo, hi * hi), math.inf)
    ends = (lo * lo, hi * hi)
    return max(0.0, math.nextafter(min(ends), -math.inf)), math.nextafter(max(ends), math.inf)


def interval_norm_squared(weights):
    """Full matrix product followed by a sum of four enclosed squares."""
    h = [[(value, value) for value in row] for row in weights[0]]
    for layer, weight in enumerate(weights[1:], 1):
        active = range(1) if (layer - 1) % 2 == 0 else range(2)
        product = []
        for i in range(2):
            row = []
            for j in range(2):
                terms = [_product((weight[i][k], weight[i][k]), h[k][j]) for k in active]
                row.append(terms[0] if len(terms) == 1 else _sum(*terms))
            product.append(row)
        h = product
    final_rows = range(1) if (len(weights) - 1) % 2 == 0 else range(2)
    total = (0.0, 0.0)
    for i in final_rows:
        for entry in h[i]:
            total = _sum(total, _square(entry))
    return max(0.0, total[0]), total[1]


def exact_norm_squared(weights):
    """Exact represented-real full matrix product; does not take a square root."""
    matrices = [[[Fraction(value) for value in row] for row in weight] for weight in weights]
    h = matrices[0]
    for layer, weight in enumerate(matrices[1:], 1):
        active = range(1) if (layer - 1) % 2 == 0 else range(2)
        h = [[sum(weight[i][k] * h[k][j] for k in active) for j in range(2)] for i in range(2)]
    final_rows = range(1) if (len(weights) - 1) % 2 == 0 else range(2)
    return sum(value * value for i in final_rows for value in h[i])


def certify_clipped_value(weights, value, absolute_error=1e-12):
    """Certify |value-min(||J_D||_2,1)| <= absolute_error; never silently repair."""
    weights = _clean(weights)
    value, absolute_error = float(value), float(absolute_error)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('The claimed clipped value must lie in [0,1].')
    if not math.isfinite(absolute_error) or not 0 < absolute_error < 1:
        raise ValueError('A positive absolute-error budget below one is required.')
    lo, hi = interval_norm_squared(weights)
    # Round the lower threshold UP and the upper threshold DOWN: this makes
    # interval acceptance stronger than the exact requested inequalities.
    lower = max(0.0, math.nextafter(value - absolute_error, math.inf))
    upper = math.nextafter(value + absolute_error, -math.inf)
    lower_squared = math.nextafter(lower * lower, math.inf) if lower else 0.0
    upper_squared = math.nextafter(upper * upper, -math.inf)
    lower_ok = lo >= lower_squared
    upper_ok = upper >= 1 or hi <= upper_squared
    if lower_ok and upper_ok:
        return {'certified': True, 'method': 'FULL_MATRIX_INTERVAL', 'absolute_error': absolute_error}
    squared = exact_norm_squared(weights)
    wanted_lo = max(Fraction(0), Fraction(value) - Fraction(absolute_error))
    wanted_hi = Fraction(value) + Fraction(absolute_error)
    passed = squared >= wanted_lo ** 2 and (wanted_hi >= 1 or squared <= wanted_hi ** 2)
    return {'certified': passed, 'method': 'FULL_MATRIX_RATIONAL', 'absolute_error': absolute_error,
            'exact_squared_numerator': str(squared.numerator),
            'exact_squared_denominator': str(squared.denominator)}

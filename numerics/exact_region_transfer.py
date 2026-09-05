"""Exact rational means for the manuscript's strict activation-region counts.

Square width n, depth L, full R^n input, independent centered Gaussian weights
with positive variance, and a final ReLU. The affine case also requires
independent centered Gaussian bias coordinates with positive variances,
independent of the weights. These are activation cells, not maximal domains
on which the final output has one affine formula.

Run this file to print the first four widths' exact constants and finite means.
It uses only the Python standard library and writes no files. Independent
rank and actual scalar-composition controls are in check_region_transfer.py.
"""
from fractions import Fraction
from math import comb
import json


def binomial_tail(d, t):
    if not isinstance(d, int) or not isinstance(t, int) or not 0 <= t <= d:
        raise ValueError("Require integers 0 <= t <= d")
    return sum(comb(d, j) for j in range(t, d + 1))


def transfer_matrix(n):
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("Width must be a positive integer")
    return [
        [
            2 ** (r - t) * comb(n, r - t) * binomial_tail(n - r + t, t)
            if r >= t else 0
            for t in range(n + 1)
        ]
        for r in range(n + 1)
    ]


def expected_strict_regions(n, depth, *, gaussian_bias=False):
    """Return the exact finite mean as a Fraction, with all activation cells kept."""
    if not isinstance(depth, int) or isinstance(depth, bool) or depth < 1:
        raise ValueError("Depth must be a positive integer")
    matrix = transfer_matrix(n)
    if gaussian_bias:
        vector = [comb(n, t) for t in range(n + 1)]
        factor = 1
    else:
        matrix = [row[1:] for row in matrix[1:]]
        vector = [comb(n - 1, t) for t in range(n)]
        factor = 2
    for _ in range(depth - 1):
        vector = [sum(a * b for a, b in zip(row, vector)) for row in matrix]
    return Fraction(factor * sum(vector), 2 ** (n * (depth - 1)))


def depth_constants(n):
    """Return kappa_n for the normalized zero-bias mean and the affine limit rho_n."""
    matrix = transfer_matrix(n)
    q, m = 2 ** n, 2 ** n - 1
    u = [Fraction(1)]
    for r in range(2, n + 1):
        u.append(
            sum(matrix[r][t] * u[t - 1] for t in range(1, r))
            / (m - matrix[r][r])
        )
    z = [Fraction(1)]
    for r in range(1, n + 1):
        z.append(sum(matrix[r][t] * z[t] for t in range(r)) / (q - matrix[r][r]))
    return Fraction(2 * q, m) * sum(u), sum(z)


if __name__ == "__main__":
    rows = []
    for width in range(1, 5):
        kappa, rho = depth_constants(width)
        rows.append({
            "width": width, "kappa": str(kappa), "affine_limit": str(rho),
            "depths": [
                {"depth": depth,
                 "zero_bias_mean": str(expected_strict_regions(width, depth)),
                 "gaussian_bias_mean": str(expected_strict_regions(
                     width, depth, gaussian_bias=True))}
                for depth in range(1, 7)
            ],
        })
    print(json.dumps(rows, indent=2))

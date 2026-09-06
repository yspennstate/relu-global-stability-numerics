"""Independent strict geometry for the rank-one masks in historical Table 3.

Certifies signs for the exact real values represented by finite binary floats.
No angle calculation, sampled inputs, Gaussian RNG, or author geometry import.
The accompanying derivation states the limited mask family and hypotheses.
"""
from fractions import Fraction
import math


def point(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('Only finite matrix entries are admitted.')
    return value, value


def add(left, right):
    lo, hi = left[0] + right[0], left[1] + right[1]
    if math.isnan(lo) or math.isnan(hi):
        return -math.inf, math.inf
    return math.nextafter(lo, -math.inf), math.nextafter(hi, math.inf)


def neg(value):
    return -value[1], -value[0]


def mul(left, right):
    products = [a * b for a in left for b in right]
    if any(math.isnan(value) for value in products):
        return -math.inf, math.inf
    return (math.nextafter(min(products), -math.inf),
            math.nextafter(max(products), math.inf))


def interval_sign(value):
    if value[0] > 0:
        return 1
    if value[1] < 0:
        return -1
    return None


def exact_sign(value):
    return (value > 0) - (value < 0)


def _exact_coefficients(matrices):
    vector = [Fraction(1), Fraction(0)]
    signs = []
    for layer in range(1, len(matrices)):
        # Previous one-based layer is odd exactly when this zero-based index is odd.
        previous = [vector[0], Fraction(0)] if layer % 2 else vector
        matrix = matrices[layer]
        vector = [sum(Fraction.from_float(float(matrix[i][j])) * previous[j]
                      for j in range(2)) for i in range(2)]
        signs.append([exact_sign(value) for value in vector])
    return signs


def exact_geometry(matrices):
    """Slow independent rational oracle, primarily for adversarial controls."""
    matrices = _validate(matrices)
    first = [[Fraction.from_float(x) for x in row] for row in matrices[0]]
    determinant = first[0][0] * first[1][1] - first[0][1] * first[1][0]
    signs = _exact_coefficients(matrices)
    return _conclusion(exact_sign(determinant), signs, used_rational=True)


def _validate(matrices):
    if len(matrices) < 2:
        raise ValueError('The certified Table 3 family has depth at least two.')
    clean = []
    for matrix in matrices:
        if len(matrix) != 2 or any(len(row) != 2 for row in matrix):
            raise ValueError('The certified Table 3 family has width exactly two.')
        clean.append([[point(value)[0] for value in row] for row in matrix])
    return clean


def _conclusion(determinant_sign, signs, *, used_rational):
    if determinant_sign == 0:
        return {'status': 'DEGENERATE_FIRST_MATRIX', 'region_count': None,
                'target_realized': None, 'used_rational': used_rational,
                'coefficient_signs': signs}
    if any(sign == 0 for pair in signs for sign in pair):
        return {'status': 'DEGENERATE_FUTURE_ROW', 'region_count': None,
                'target_realized': None, 'used_rational': used_rational,
                'coefficient_signs': signs}
    # signs[0] belongs to one-based layer two, whose target is (+,+).
    realized = all(pair == ([1, 1] if index % 2 == 0 else [1, -1])
                   for index, pair in enumerate(signs))
    return {'status': 'CERTIFIED', 'region_count': 4,
            'target_realized': realized, 'used_rational': used_rational,
            'coefficient_signs': signs}


def certified_geometry(matrices):
    """Outward interval certificate, falling back to exact represented arithmetic."""
    matrices = _validate(matrices)
    first = matrices[0]
    determinant = add(mul(point(first[0][0]), point(first[1][1])),
                      neg(mul(point(first[0][1]), point(first[1][0]))))
    determinant_sign = interval_sign(determinant)
    vector = [point(1), point(0)]
    signs = []
    uncertain = determinant_sign is None
    for layer in range(1, len(matrices)):
        previous = [vector[0], point(0)] if layer % 2 else vector
        matrix = matrices[layer]
        vector = [add(mul(point(matrix[i][0]), previous[0]),
                      mul(point(matrix[i][1]), previous[1])) for i in range(2)]
        pair = [interval_sign(value) for value in vector]
        uncertain = uncertain or None in pair
        signs.append(pair)
    if uncertain:
        return exact_geometry(matrices)
    return _conclusion(determinant_sign, signs, used_rational=False)


def rank_one_norm(matrices):
    """Independent floating-point norm, explicitly not an interval certificate."""
    matrices = _validate(matrices)
    vector = [1.0, 0.0]
    for layer in range(1, len(matrices)):
        previous = [vector[0], 0.0] if layer % 2 else vector
        matrix = matrices[layer]
        vector = [math.fsum(matrix[i][j] * previous[j] for j in range(2))
                  for i in range(2)]
    if len(matrices) % 2:
        vector[1] = 0.0
    return math.hypot(*vector) * math.hypot(*matrices[0][0])

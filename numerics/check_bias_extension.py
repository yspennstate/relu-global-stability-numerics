"""Finite author-verification controls; these do not certify the full theorems."""
import itertools
import math
import unittest
import numpy as np


def middle(lo, hi):
    if math.isinf(lo) and math.isinf(hi):
        return 0.
    if math.isinf(lo):
        return hi - max(1., abs(hi))
    if math.isinf(hi):
        return lo + max(1., abs(lo))
    return (lo + hi) / 2


def scalar_lipschitz(weights, biases):
    """Enumerate all affine intervals directly, including zero-slope pieces."""
    pieces = [(-math.inf, math.inf, 1., 0.)]
    for w, b in zip(weights, biases):
        following = []
        for lo, hi, slope, offset in pieces:
            a, c = w*slope, w*offset+b
            cuts = [lo, hi]
            if a and lo < -c/a < hi:
                cuts.insert(1, -c/a)
            for left, right in zip(cuts, cuts[1:]):
                active = a*middle(left, right)+c > 0
                following.append((left, right, a if active else 0., c if active else 0.))
        pieces = following
    return max(abs(p[2]) for p in pieces)


def forward(weights, biases, s):
    x = np.asarray(s, dtype=float)
    for w, b in zip(weights, biases):
        x = np.maximum(w @ x+b, 0.)
    return x


def affine_stack(weights, biases, mask):
    h, c = weights[0], biases[0]
    hs, cs = [h], [c]
    for i in range(1, len(weights)):
        h = weights[i]*mask[i-1]*h
        c = weights[i]*mask[i-1]*c+biases[i]
        hs.append(h)
        cs.append(c)
    return np.array(hs), np.array(cs), mask[-1]*h


def affine_patterns(h, c):
    if np.any((h == 0) & (c == 0)):
        return set()
    roots = sorted(set(-offset/slope for slope, offset in zip(h, c) if slope))
    cuts = [-math.inf, *roots, math.inf]
    return {tuple((h*middle(lo, hi)+c > 0).astype(int))
            for lo, hi in zip(cuts, cuts[1:])}


class BiasControls(unittest.TestCase):
    def test_scalar_lipschitz_domination(self):
        checked = 0
        for weights in itertools.product((-1., 0., 1.), repeat=3):
            k0 = scalar_lipschitz(weights, (0.,)*3)
            for biases in itertools.product((-2., 0., 3.), repeat=3):
                self.assertGreaterEqual(scalar_lipschitz(weights, biases), k0)
                checked += 1
        self.assertEqual(checked, 729)
        # A dead zero-bias network can acquire a nonconstant bounded piece.
        self.assertEqual(scalar_lipschitz((1., -1.), (0., 0.)), 0.)
        self.assertEqual(scalar_lipschitz((1., -1.), (0., 1.)), 1.)

    def test_uniform_recession_bound(self):
        rng = np.random.default_rng(20260905)
        for dead_first in (False, True):
            for _ in range(12):
                ws = [rng.normal(size=(3, 3)) for _ in range(4)]
                if dead_first:
                    ws[0] *= 0
                bs = [rng.normal(size=3) for _ in ws]
                norms = [np.linalg.norm(w, 2) for w in ws]
                bound = 0.
                for op, b in zip(norms, bs):
                    bound = op*bound+np.linalg.norm(b)
                expanded = sum(np.linalg.norm(bs[i])*math.prod(norms[i+1:]) for i in range(4))
                self.assertAlmostEqual(bound, expanded, places=10)
                for s in rng.normal(size=(8, 3)):
                    f0 = forward(ws, [np.zeros(3)]*4, s)
                    for scale in (1., 10., 1e4):
                        fb = forward(ws, bs, scale*s)
                        self.assertLessEqual(np.linalg.norm(fb/scale-f0), bound/scale+1e-9)

    def test_full_affine_sign_orbit(self):
        weights = np.array([.7, -1.3, 2.1])
        for biases in (np.zeros(3), np.array([.5, -.8, .2])):
            for mask in itertools.product((0, 1), repeat=3):
                h, c, j = affine_stack(weights, biases, mask)
                patterns = affine_patterns(h, c)
                hits = 0
                for signs in itertools.product((-1., 1.), repeat=3):
                    sig = np.array(signs)
                    wg = sig*weights*np.r_[1., sig[:-1]]
                    hg, cg, jg = affine_stack(wg, sig*biases, mask)
                    np.testing.assert_allclose(hg, sig*h, atol=1e-14)
                    np.testing.assert_allclose(cg, sig*c, atol=1e-14)
                    self.assertAlmostEqual(jg, sig[-1]*j)
                    pg = affine_patterns(hg, cg)
                    self.assertEqual(len(pg), len(patterns))
                    hits += mask in pg
                self.assertEqual(hits, len(patterns))
        # The old correlated-bias law is excluded by the new hypothesis.
        b = np.array([0., 1., 1.])
        self.assertNotEqual(b[1]*b[2], (-b[1])*b[2])

    def test_reference_profile_including_zero_bias(self):
        for sigma in (0., .5, 1., 2.):
            q = 1.
            for t in range(9):
                self.assertAlmostEqual(q, sigma*sigma+(1-sigma*sigma)*2.**(-t))
                q = (q+sigma*sigma)/2


if __name__ == '__main__':
    unittest.main(verbosity=2)

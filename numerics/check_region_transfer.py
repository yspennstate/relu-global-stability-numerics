"""Exact geometry controls for the finite region-transfer theorem.

This ports the retained author checks to a portable standard-library script.
Direct normal ranks and actual scalar compositions are independent calculation
paths from matrix multiplication; finite checks do not prove genericity.
Run from any directory. The JSON record goes to stdout; no files are written.
"""
from fractions import Fraction as F
import datetime as dt
import hashlib
import itertools
import json
from pathlib import Path
import random
import time

from exact_region_transfer import expected_strict_regions, depth_constants


def rank(rows, n):
    a = [[F(x) for x in row] for row in rows]
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        for i in range(r + 1, len(a)):
            if a[i][c]:
                ratio = a[i][c] / a[r][c]
                a[i] = [x - ratio * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == n:
            break
    return r


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def predicted_rank(masks, selections, *, remove_overlap=True):
    result = selections[-1].bit_count()
    for j in range(len(masks) - 1, -1, -1):
        retained = masks[j] & ~selections[j] if remove_overlap else masks[j]
        result = selections[j].bit_count() + min(result, retained.bit_count())
    return result


def scalar_regions(weights, biases):
    """Compose exact affine pieces, retaining earlier activation distinctions."""
    pieces = [(None, None, F(1), F(0))]
    for w, b in zip(weights, biases):
        new = []
        for lo, hi, a, c in pieces:
            p, q = w * a, w * c + b
            if not p:
                if not q:
                    continue  # A forced zero violates strictness at this layer.
                new.append((lo, hi, F(0), max(F(0), q)))
                continue
            root = -q / p
            inside = (lo is None or lo < root) and (hi is None or root < hi)
            intervals = [(lo, root), (root, hi)] if inside else [(lo, hi)]
            for left, right in intervals:
                sample = ((left + right) / 2 if left is not None and right is not None
                          else left + 1 if left is not None
                          else right - 1 if right is not None else F(0))
                value = p * sample + q
                assert value != 0
                new.append((left, right, p if value > 0 else F(0), q if value > 0 else F(0)))
        pieces = new
    return len(pieces)


def check():
    start = time.perf_counter()
    rank_checks, overlap_mutation_failures = 0, 0
    geometry = []
    for seed in (2026090601, 2026090602):
        rng = random.Random(seed)
        for n, depth in ((2, 3), (3, 2)):
            weights = [
                [[rng.randrange(-100000, 100001) for _ in range(n)] for _ in range(n)]
                for _ in range(depth)
            ]
            assert all(rank(w, n) == n for w in weights)
            central = affine = 0
            zero_mask_affine_counts = []
            for masks in itertools.product(range(2 ** n), repeat=depth - 1):
                block = weights[0]
                blocks = [block]
                for j, mask in enumerate(masks):
                    block = multiply(weights[j + 1], [
                        row if mask >> i & 1 else [0] * n for i, row in enumerate(block)
                    ])
                    blocks.append(block)
                affine_for_mask = central_for_mask = 0
                for selections in itertools.product(range(2 ** n), repeat=depth):
                    rows = [
                        row for block, selection in zip(blocks, selections)
                        for i, row in enumerate(block) if selection >> i & 1
                    ]
                    actual = rank(rows, n)
                    assert actual == predicted_rank(masks, selections)
                    rank_checks += 1
                    overlap_mutation_failures += (
                        actual != predicted_rank(masks, selections, remove_overlap=False)
                    )
                    central_for_mask += (-1) ** (len(rows) - actual)
                    affine_for_mask += int(len(rows) == actual)
                if masks[0] == 0:
                    assert central_for_mask == 0 and affine_for_mask == 2 ** n
                    zero_mask_affine_counts.append(affine_for_mask)
                central += central_for_mask
                affine += affine_for_mask
            zero_mean = F(central, 2 ** (n * (depth - 1)))
            bias_mean = F(affine, 2 ** (n * (depth - 1)))
            assert zero_mean == expected_strict_regions(n, depth)
            assert bias_mean == expected_strict_regions(n, depth, gaussian_bias=True)
            assert zero_mean != F(central, 2 ** (n * depth))  # lost final-mask factor
            assert zero_mask_affine_counts  # The empty-affine versus zero-row branch was reached.
            geometry.append({
                "seed": seed, "n": n, "depth": depth, "weights": weights,
                "central_mean": str(zero_mean), "affine_mean": str(bias_mean),
                "zero_mask_affine_counts": zero_mask_affine_counts,
            })
    assert overlap_mutation_failures > 0

    weights = (2, -3, 5, -7, 11, -13, 17, -19)
    biases = (23, 29, 31, 37, 41, 43, 47, 53)
    scalar, networks = [], 0
    for depth in range(1, 9):
        for affine_case in (False, True):
            total = 0
            for signs in itertools.product((-1, 1), repeat=depth):
                ws = [
                    signs[j] * weights[j] * (signs[j - 1] if j else 1)
                    for j in range(depth)
                ]
                bs = [signs[j] * biases[j] if affine_case else 0 for j in range(depth)]
                total += scalar_regions(ws, bs)
                networks += 1
            observed = F(total, 2 ** depth)
            closed_form = 3 - F(2) ** (1 - depth) if affine_case else F(2) ** (2 - depth)
            assert observed == closed_form == expected_strict_regions(
                1, depth, gaussian_bias=affine_case
            )
            scalar.append({"depth": depth, "gaussian_bias": affine_case,
                           "networks": 2 ** depth, "regions_sum": total, "mean": str(observed)})

    for depth in range(1, 26):
        a, b = F(3, 4) ** depth, F(1, 4) ** depth
        assert expected_strict_regions(2, depth) == 8 * (a - b)
        assert expected_strict_regions(2, depth, gaussian_bias=True) == 21 - 24 * a + 4 * b
    constants = [
        {"n": n, "kappa": str(depth_constants(n)[0]),
         "affine_limit": str(depth_constants(n)[1])}
        for n in range(1, 9)
    ]
    assert [(row["kappa"], row["affine_limit"]) for row in constants[:3]] == [
        ("4", "3"), ("8", "21"), ("240/7", "279")
    ]
    here = Path(__file__).resolve().parent
    return {
        "status": "AUTHOR_EXACT_GEOMETRY_CONTROLS_PASS",
        "observed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter() - start,
        "rank_checks": rank_checks, "actual_scalar_networks": networks,
        "width_two_closed_form_checks": 50,
        "mutations_rejected": {
            "forget_selected_retained_overlap": overlap_mutation_failures,
            "lose_final_mask_factor": len(geometry),
            "treat_empty_affine_constraint_as_central_zero_row": len(geometry),
        },
        "geometry": geometry, "scalar_networks": scalar, "constants": constants,
        "source_hashes": {
            name: hashlib.sha256((here / name).read_bytes()).hexdigest()
            for name in ("exact_region_transfer.py", "check_region_transfer.py")
        },
        "limits": [
            "Finite integer tuples and sign orbits do not prove Gaussian genericity.",
            "Exact normal ranks and scalar compositions are distinct from transfer multiplication.",
            "These are author controls, separate from independent referee derivations.",
            "Only strict activation cells are counted; maximal output-affine regions may merge.",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(check(), indent=2))

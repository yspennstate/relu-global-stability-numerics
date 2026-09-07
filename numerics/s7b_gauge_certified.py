"""Reproduce the corrected Table 3 gauge-identity experiment.

The seed, six row sizes, alternating masks and Gaussian calls match the retained
historical s7b generator. Geometry is certified for the exact real values of the
sampled binary64 weights. Each clipped functional is evaluated numerically and
then checked against an absolute error budget of 1e-12. Neither statement is a
certificate for the ideal Gaussian sampler or for Monte Carlo sampling error.

Usage: python -B numerics/s7b_gauge_certified.py --output reproduced_table3.json
For a short controls run, add --pilot 1000; this changes the row sizes and is
explicitly labeled PILOT_ONLY. The historical JSON is never overwritten.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys

os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
import numpy as np
from table3_certified_geometry import certified_geometry
from table3_norm_error_certificate import certify_clipped_value

ROWS = ((2, 'one', 200000), (3, 'one', 200000), (4, 'one', 150000),
        (5, 'one', 120000), (3, 'clip', 200000), (4, 'clip', 150000))
SEED = 8081
NORM_ABSOLUTE_ERROR = 1e-12


def jacobian(weights, masks):
    # Preserve the historical floating matrix-product order for the functional.
    h = weights[0].copy()
    for layer in range(1, len(weights)):
        h = weights[layer] @ np.diag(masks[layer - 1].astype(float)) @ h
    return np.diag(masks[-1].astype(float)) @ h


def statistics(functional, target, depth):
    count = len(functional)
    if count < 2 or len(target) != count:
        raise ValueError('Paired arrays must have equal length at least two.')
    probability = 4.0 ** (1 - depth)
    left = [float(value) * bool(hit) for value, hit in zip(functional, target)]
    right = [float(value) * probability for value in functional]
    difference = [a - b for a, b in zip(left, right)]
    lhs, rhs = math.fsum(left) / count, math.fsum(right) / count
    mean = math.fsum(difference) / count
    se = math.sqrt(math.fsum((value - mean) ** 2 for value in difference) / (count * (count - 1)))
    return {'lhs_mean': lhs, 'rhs_mean': rhs, 'paired_diff': mean, 'paired_se': se,
            'paired_z': mean / se if se else None, 'ratio': lhs / rhs if rhs else None,
            'functional_mean': math.fsum(float(value) for value in functional) / count,
            'functional_second_moment': math.fsum(float(value) ** 2 for value in functional) / count,
            'target_successes': sum(bool(hit) for hit in target)}


def experiment(pilot=None, progress=None):
    if pilot is not None and (type(pilot) is not int or not 2 <= pilot <= 10000):
        raise ValueError('A pilot uses an integer from 2 through 10000 draws per row.')
    rng = np.random.default_rng(SEED)
    results = {}
    for depth, kind, original_count in ROWS:
        count = original_count if pilot is None else pilot
        masks = [np.array([1, 0]) if layer % 2 == 0 else np.array([1, 1]) for layer in range(depth)]
        functional, target = np.empty(count), np.empty(count, dtype=np.bool_)
        geometry_fallbacks = norm_fallbacks = 0
        for draw in range(count):
            weights = [rng.normal(0.0, math.sqrt(0.5), size=(2, 2)) for _ in range(depth)]
            geometry = certified_geometry(weights)
            if geometry['status'] != 'CERTIFIED' or geometry['region_count'] != 4:
                raise ArithmeticError(f'Unresolved represented geometry at depth={depth}, F={kind}, draw={draw}; no resampling is allowed.')
            geometry_fallbacks += geometry['used_rational']
            target[draw] = geometry['target_realized']
            value = 1.0 if kind == 'one' else min(1.0, float(np.linalg.norm(jacobian(weights, masks), 2)))
            if kind == 'clip':
                certificate = certify_clipped_value(weights, value, NORM_ABSOLUTE_ERROR)
                if not certificate['certified']:
                    raise ArithmeticError(f'Clipped norm exceeds error budget at depth={depth}, draw={draw}; no resampling is allowed.')
                norm_fallbacks += certificate['method'] == 'FULL_MATRIX_RATIONAL'
            functional[draw] = value
            if progress is not None and (draw + 1) % 5000 == 0:
                progress({'L': depth, 'F': kind, 'completed': draw + 1, 'trials': count})
        key = f'L{depth}_{kind}'
        results[key] = {'L': depth, 'F': kind, 'trials': count, 'seed': SEED,
            'masks': [mask.tolist() for mask in masks], 'strict_orthants_per_draw': 4,
            'exact_unweighted_probability': f'1/{4 ** (depth - 1)}',
            'geometry_rational_fallbacks': geometry_fallbacks,
            'norm_rational_fallbacks': norm_fallbacks,
            'functional_absolute_error_bound': 0.0 if kind == 'one' else NORM_ABSOLUTE_ERROR,
            'functional_binary64_sha256': hashlib.sha256(functional.astype('<f8').tobytes()).hexdigest(),
            'target_uint8_sha256': hashlib.sha256(target.astype('u1').tobytes()).hexdigest(),
            **statistics(functional, target, depth)}
        print(json.dumps({'row': key, **results[key]}), flush=True)
    return {'schema': 'relu-table3-corrected-v1', 'status': 'COMPLETE' if pilot is None else 'PILOT_ONLY',
            'seed': SEED, 'numpy_version': np.__version__, 'python_version': sys.version,
            'norm_absolute_error_bound': NORM_ABSOLUTE_ERROR,
            'geometry_scope': 'Strict formal-stack orthants for the exact real values represented by binary64 weights, with the fixed alternating width-two masks only.',
            'numerical_scope': 'Clipped functional error at most 1e-12, certified for represented weights under IEEE binary64 basic arithmetic; paired standard errors measure Monte Carlo variation.',
            'draws': sum(row['trials'] for row in results.values()), 'rows': results,
            'rng_state_after': rng.bit_generator.state}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--pilot', type=int)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error('The output already exists; choose a new filename to preserve prior results.')
    if not output.parent.is_dir():
        parser.error('The output directory must already exist.')
    result = experiment(args.pilot, progress=lambda row: print(json.dumps({'progress': row}), flush=True))
    # Exclusive creation protects an output that appeared while the experiment ran.
    with output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()

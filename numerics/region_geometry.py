"""Strict cell enumeration, using floating-point boundary arcs or LP.

Strict means every preactivation is nonzero. Ordinary full-dimensional cells
allow identically zero preactivations with inactive mask bits. Boundary-only
patterns belong to neither count. No finite angle grid is called exact.
"""
import itertools
import math
import warnings
import numpy as np
from scipy.optimize import linprog, OptimizeWarning

TAU = 2*math.pi


def cuts_for(A, a=0., b=TAU):
    cuts = {a,b}
    for row in A:
        if not np.any(row):
            continue
        t = math.atan2(row[0],-row[1]) % TAU
        for q in (t,(t+math.pi)%TAU):
            if a < q < b:
                cuts.add(q)
    return sorted(cuts)


def midpoint_sign(A,a,b):
    t = a+(b-a)/2
    if not a < t < b:
        raise ArithmeticError('Unresolved narrow angular cell')
    x = np.array([math.cos(t),math.sin(t)])
    y = A@x
    zero = np.all(A==0,axis=1)
    if np.any(y[~zero]==0):
        raise ArithmeticError('Boundary midpoint requires higher precision')
    return tuple((y>0).astype(int)),x


def planar_patterns(A):
    A = np.asarray(A,dtype=float)
    if np.any(np.all(A==0,axis=1)):
        return set()
    cuts = cuts_for(A)
    return {midpoint_sign(A,a,b)[0] for a,b in zip(cuts,cuts[1:])}


def planar_network(Ws):
    cells = [(0.,TAU,np.eye(2),(),True,2)]
    for W in Ws:
        new = []
        for a,b,J,mask,strict,rank in cells:
            H = W@J
            # A rank-one ReLU image is a fixed positive ray on this cell;
            # later layers cannot split it. Recomputing inherited angles
            # would manufacture rounding-width slivers.
            cuts = cuts_for(H,a,b) if rank==2 else [a,b]
            for lo,hi in zip(cuts,cuts[1:]):
                d,_ = midpoint_sign(H,lo,hi)
                next_rank = min(rank,sum(d),2 if np.linalg.det(W)!=0 else 1)
                new.append((lo,hi,np.asarray(d)[:,None]*H,mask+d,
                            strict and not bool(np.any(np.all(H==0,axis=1))),next_rank))
        cells = new
    ordinary,strict_masks,witnesses = {},{},{}
    for a,b,J,mask,strict,rank in cells:
        ordinary[mask] = J
        t = a+(b-a)/2
        witnesses[mask] = [math.cos(t),math.sin(t)]
        if strict:
            strict_masks[mask] = J
    return dict(strict=strict_masks,ordinary=ordinary,witnesses=witnesses,
                K=max((float(np.linalg.norm(J,2)) for J in ordinary.values()),default=0.))


def formal_rows(Ws,mask):
    n = Ws[0].shape[0]
    J,rows = np.eye(n),[]
    for ell,W in enumerate(Ws):
        H = W@J
        d = np.asarray(mask[ell*n:(ell+1)*n])
        rows.extend((2*d-1)[:,None]*H)
        J = d[:,None]*H
    return np.asarray(rows),J


def lp_network(Ws):
    """Exhaust masks; a strict homogeneous system is scalable to A x >= 1.

    Normalize nonzero rows, check returned witnesses, fail on unresolved LPs.
    Zero rows exclude strictness; ordinary cells permit them only at bit zero.
    """
    n,L = Ws[0].shape[0],len(Ws)
    strict_masks,ordinary,witnesses = {},{},{}
    for mask in itertools.product((0,1),repeat=n*L):
        A,J = formal_rows(Ws,mask)
        norms = np.linalg.norm(A,axis=1)
        zero = norms==0
        if np.any(np.asarray(mask)[zero]!=0):
            continue
        B = A[~zero]/norms[~zero,None]
        if len(B):
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',category=OptimizeWarning,message='Unrecognized options detected')
                sol = linprog(np.zeros(n),A_ub=-B,b_ub=-np.ones(len(B)),
                              bounds=[(None,None)]*n,method='highs',
                              options={'threads':1,'primal_feasibility_tolerance':1e-9})
            if sol.status==2:
                continue
            if not sol.success or np.min(B@sol.x)<1-1e-7:
                raise ArithmeticError('LP unresolved or witness violates strict margin')
            x = sol.x.tolist()
        else:
            x = [1.]+[0.]*(n-1)
        ordinary[mask],witnesses[mask] = J,x
        if not zero.any():
            strict_masks[mask] = J
    return dict(strict=strict_masks,ordinary=ordinary,witnesses=witnesses,
                K=max((float(np.linalg.norm(J,2)) for J in ordinary.values()),default=0.))

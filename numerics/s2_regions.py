"""Tables 2 and 7: strict cells and Lipschitz constants (revision 2026-09-05).
Boundary arcs for width two; exhaustive strict LP for higher-dimensional
region counts. Threshold samples use 40000 directions for n>=3.
Every row has a separate seed and per-network counts or K values.
"""
import json, math
from pathlib import Path
import numpy as np
from region_geometry import planar_network, lp_network
OUT = Path(__file__).resolve().parent / 'data'
SEED = 20260905

def cover(n,L):
    return 2*sum(math.comb(n*L-1,j) for j in range(n))

def sample_K(Ws,rng,npts=40000):
    n=Ws[0].shape[0]
    X=rng.normal(size=(n,npts)); masks=[]
    for W in Ws:
        H=W@X; masks.append(H>0); X=np.maximum(H,0)
    code=np.zeros(npts,dtype=np.uint64)
    for bit in np.concatenate(masks,axis=0): code=(code<<1)|bit.astype(np.uint64)
    _,first=np.unique(code,return_index=True)
    best=0.
    for idx in first:
        J=np.eye(n)
        for W,D in zip(Ws,masks): J=D[:,idx,None]*(W@J)
        best=max(best,float(np.linalg.norm(J,2)))
    return best

def main():
    OUT.mkdir(exist_ok=True)
    res=dict(schema='strict-regions-v2',seed=SEED,
        note='Enumeration uses floating point, not exact arithmetic. Strict requires all preactivations nonzero. Forced-zero full-dimensional cells are separate. C bounds expected strict count.',
        region_counts=[],K_at_alpha1={})
    for n,L,trials in [(2,L,60) for L in range(1,7)]+[(3,2,25),(3,3,25),(4,2,25)]:
        seed=SEED+100*n+L; rng=np.random.default_rng(seed); draws=[]
        for _ in range(trials):
            Ws=[rng.normal(0,1/math.sqrt(n),(n,n)) for _ in range(L)]
            ans=planar_network(Ws) if n==2 else lp_network(Ws)
            keys=sorted(ans['strict'])
            draws.append(dict(weights=[W.tolist() for W in Ws],
                strict_masks=[[int(b) for b in m] for m in keys],
                strict_witnesses=[ans['witnesses'][m] for m in keys],
                strict_count=len(keys),ordinary_count=len(ans['ordinary']),
                forced_zero_count=len(ans['ordinary'])-len(keys)))
        counts=[d['strict_count'] for d in draws]
        row=dict(n=n,L=L,trials=trials,seed=seed,
            method='boundary arcs' if n==2 else 'all-mask LP',formal_masks=2**(n*L),
            expected_count_bound_C=cover(n,L),strict_mean=float(np.mean(counts)),
            strict_min=min(counts),strict_max=max(counts),
            forced_zero_mean=float(np.mean([d['forced_zero_count'] for d in draws])),draws=draws)
        res['region_counts'].append(row)
        (OUT/'s2_regions.json').write_text(json.dumps(res,indent=1,allow_nan=False)+'\n',encoding='utf-8')
        print('regions',n,L,row['strict_mean'],row['strict_min'],row['strict_max'],flush=True)
    for n,L,trials in [(2,1,400),(2,2,400),(2,3,400),(2,4,300),(3,2,150),(3,3,120),(5,3,60),(8,3,40)]:
        seed=SEED+10000+100*n+L; rng=np.random.default_rng(seed); Ks=[]
        for _ in range(trials):
            Ws=[rng.normal(0,1/math.sqrt(n),(n,n)) for _ in range(L)]
            Ks.append(planar_network(Ws)['K'] if n==2 else sample_K(Ws,rng))
        median=float(np.median(Ks))
        res['K_at_alpha1'][f'n{n}_L{L}']=dict(n=n,L=L,trials=trials,seed=seed,
            method='boundary arcs' if n==2 else '40000 sampled directions',
            median_K1=median,q10=float(np.quantile(Ks,.1)),q90=float(np.quantile(Ks,.9)),
            alpha_from_median_K=median**(-2/L) if median else None,K_values=Ks)
        print('K',n,L,median,flush=True)
        (OUT/'s2_regions.json').write_text(json.dumps(res,indent=1,allow_nan=False)+'\n',encoding='utf-8')
    (OUT/'s2_regions.json').write_text(json.dumps(res,indent=1,allow_nan=False)+'\n',encoding='utf-8')
if __name__=='__main__': main()

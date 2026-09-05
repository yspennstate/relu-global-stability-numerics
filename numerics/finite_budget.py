"""Tables 5/6 from the optimized finite-width moment exponent g.

The old grid scripts retain the earlier unoptimized h bound for comparison.
This implementation solves for the exact scalar optimizer in z and checks
independent minimization in theta. No Monte Carlo is involved.
"""
import json
import math
from pathlib import Path
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gammaln, logsumexp
import numpy as np

def H(u):
    return -u*math.log(u)-(1-u)*math.log1p(-u) if 0<u<1 else 0.

def h(alpha,t):
    """Historical, unoptimized bound; not the current theorem's exponent."""
    return t*math.log(alpha/(2*math.e))+math.log((1+(1-4*t)**-.5)/2)

def g(alpha,t):
    if alpha<=0 or t<0: raise ValueError('alpha must be positive and theta nonnegative')
    if t<1e-4:
        # Avoid cancellation of the opposite linear terms at tiny theta.
        # Coefficients follow by expanding the explicit radical formula.
        return t*math.log(alpha/2)+t*t*(2.5+t*(-31/6+t*(209/12+t*(-1471/20+t*2125/6))))
    z=2/(1+math.sqrt(1+8*t))
    if t<0.01:
        # 1-z=2 theta z^2 avoids subtracting almost equal floating-point values.
        w=2*t*z*z
        lz=math.log1p(-w)
        lhalf=math.log1p(-w/2)
        return t*math.log(alpha/2)+t*(-1-2*lz-lhalf)+lhalf-lz
    # log1p preserves accuracy when theta is small and z is close to one.
    return t*(math.log(alpha)-1-2*math.log(z)-math.log1p(z))+math.log1p((1-z)/(2*z))

def dg(alpha,t):
    z=2/(1+math.sqrt(1+8*t))
    return math.log(alpha)-2*math.log(z)-math.log1p(z)

def B(L,alpha,t):
    return L*H(1/L)+math.log(5)+t*math.log(4)+L*g(alpha,t)

def optimizer(L,alpha):
    if L<2 or alpha<=0: raise ValueError('L>=2 and alpha>0 required')
    if math.log(4)+L*math.log(alpha/2)>=0:return 0.
    target=alpha*math.exp(math.log(4)/L)
    z=brentq(lambda z:z*z*(1+z)-target,0.,1.,xtol=5e-15)
    return (1-z)/(2*z*z)

def rate(L,alpha):
    return -B(L,alpha,optimizer(L,alpha))

def alpha_star(L):
    entropy=L*H(1/L)+math.log(5)
    z=brentq(lambda z:(1-z)/(2*z*z)-math.log1p((1-z)/(2*z))-entropy/L,1e-12,1.,xtol=5e-15)
    return z*z*(1+z)*math.exp(-math.log(4)/L)

def depth(alpha):
    if not 0<alpha<2:raise ValueError('A finite sufficient depth requires 0<alpha<2')
    if rate(2,alpha)>0:return 2
    lo,hi=2,4
    while rate(hi,alpha)<=0: lo,hi=hi,2*hi
    while hi-lo>1:
        mid=(lo+hi)//2
        if rate(mid,alpha)>0: hi=mid
        else: lo=mid
    return hi

def exact_g(alpha,t):
    def phi(u):
        return H(u)-math.log(2)+t*math.log(2*alpha)+(u/2+t)*math.log(u/2+t)-(u/2)*math.log(u/2)-t
    opt=minimize_scalar(lambda u:-phi(u),bounds=(1e-12,1-1e-12),method='bounded',options={'xatol':1e-14})
    return -opt.fun

def main():
    depths=[10,30,100,300,1000,3000,10000,100000,1000000]
    res={'schema':'optimized-finite-g-v1','thresholds':[], 'rates':[], 'depths':{str(a):depth(a) for a in [1,1.4,1.7,1.9,1.95,1.99]}}
    for L in depths:
        res['thresholds'].append(dict(L=L,alpha_star=alpha_star(L),
            asymptotic=2*math.exp(-math.sqrt(10*math.log(L)/L)),
            asymptotic_A=2*math.exp(-math.sqrt(10*(math.log(L)+2+math.log(5))/L))))
    for a in [.5,1.,1.5,1.8,1.9]:
        vals=[]
        for L in [100,300,1000,3000,10000,100000]:
            c=rate(L,a)
            vals.append(dict(L=L,c=c,n_001=math.ceil(math.log(200)/c) if c>0 else None))
        res['rates'].append(dict(alpha=a,values=vals))
    for L,a in [(100,1),(1000,1),(100000,1.9)]:
        opt=minimize_scalar(lambda t:B(L,a,t),bounds=(1e-12,10),method='bounded',options={'xatol':1e-14})
        assert abs(rate(L,a)+opt.fun)<1e-7
    n=10000;t=.1;m=np.arange(1,n+1,dtype=float)
    finite=logsumexp(gammaln(n+1)-gammaln(m+1)-gammaln(n-m+1)-n*math.log(2)+t*n*math.log(4/n)+gammaln(m/2+t*n)-gammaln(m/2))/n
    res['scalar_check']=dict(alpha=2,theta=t,legacy_h=h(2,t),g=g(2,t),variational_g=exact_g(2,t),finite_n=n,finite_exponent=float(finite))
    assert abs(g(2,t)-exact_g(2,t))<1e-12
    assert res['scalar_check']['legacy_h']>res['scalar_check']['g']+.01
    assert abs(finite-exact_g(2,t))<3e-5
    out=Path(__file__).resolve().parent/'data'/'finite_budget.json'
    out.write_text(json.dumps(res,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'depths':res['depths'],'scalar':res['scalar_check'],'rate_1000_1':rate(1000,1)},indent=2))

if __name__=='__main__': main()

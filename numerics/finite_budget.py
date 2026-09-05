"""Fast derivative-root reproduction of Tables 5/6 and scalar h-versus-g.

The old grid scripts remain available for comparison. This implementation
finds stationary points with bracketing roots and checks a separate scalar
minimization for representative rates. No Monte Carlo is involved.
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
    return t*math.log(alpha/(2*math.e))+math.log((1+(1-4*t)**-.5)/2)

def dh(alpha,t):
    z=math.sqrt(1-4*t)
    return math.log(alpha/(2*math.e))+2/(z*z*(1+z))

def B(L,alpha,t):
    return L*H(1/L)+math.log(5)+t*math.log(4)+L*h(alpha,t)

def rate(L,alpha):
    if math.log(4)+L*dh(alpha,0)>=0:
        return -L*H(1/L)-math.log(5)
    t=brentq(lambda t:math.log(4)+L*dh(alpha,t),0,.249999999,xtol=1e-14)
    return -B(L,alpha,t)

def alpha_star(L):
    entropy=L*H(1/L)+math.log(5)
    t=brentq(lambda t:L*(t*dh(2,t)-h(2,t))-entropy,1e-10,.249999999,xtol=1e-14)
    return 2*math.exp(-B(L,2,t)/(L*t))

def depth(alpha):
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
    res={'thresholds':[], 'rates':[], 'depths':{str(a):depth(a) for a in [1,1.4,1.7,1.9,1.95,1.99]}}
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
        opt=minimize_scalar(lambda t:B(L,a,t),bounds=(1e-10,.249999),method='bounded',options={'xatol':1e-14})
        assert abs(rate(L,a)+opt.fun)<1e-7
    n=10000;t=.1;m=np.arange(1,n+1,dtype=float)
    finite=logsumexp(gammaln(n+1)-gammaln(m+1)-gammaln(n-m+1)-n*math.log(2)+t*n*math.log(4/n)+gammaln(m/2+t*n)-gammaln(m/2))/n
    res['scalar_check']=dict(alpha=2,theta=t,h=h(2,t),g=exact_g(2,t),finite_n=n,finite_exponent=float(finite))
    assert res['scalar_check']['h']>res['scalar_check']['g']+.01
    assert abs(finite-exact_g(2,t))<3e-5
    out=Path(__file__).resolve().parent/'data'/'finite_budget.json'
    out.write_text(json.dumps(res,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'depths':res['depths'],'scalar':res['scalar_check'],'rate_1000_1':rate(1000,1)},indent=2))

if __name__=='__main__': main()

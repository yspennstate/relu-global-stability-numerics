"""Independent high-precision numerical checks; these supplement the proofs."""
import json
import math
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from finite_budget import g

mp.mp.dps=80
DATA=Path(__file__).resolve().parent/'data'

def root(f,a,b):
    a,b=mp.mpf(a),mp.mpf(b)
    fa=f(a)
    assert fa*f(b)<=0
    for _ in range(285):
        c=(a+b)/2;fc=f(c)
        if (fa>0)==(fc>0):a,fa=c,fc
        else:b=c
    return (a+b)/2

def phi(a,t,u):
    q=u/2
    return (-u*mp.log(u)-(1-u)*mp.log(1-u)-mp.log(2)+t*mp.log(2*a)
            +(q+t)*mp.log(q+t)-q*mp.log(q)-t)

def variational(a,t):
    # Solve the active-fraction equation, independently of the radical formula.
    u=root(lambda u:mp.log((1-u)/u)+mp.log((u/2+t)/(u/2))/2,
           mp.mpf('1e-70'),1-mp.mpf('1e-70'))
    return phi(a,t,u)

def exact_log(n,a,p):
    vals=[mp.loggamma(n+1)-mp.loggamma(m+1)-mp.loggamma(n-m+1)
          +mp.loggamma(mp.mpf(m)/2+p)-mp.loggamma(mp.mpf(m)/2) for m in range(1,n+1)]
    top=max(vals)
    return -n*mp.log(2)+p*mp.log(2*a/n)+top+mp.log(mp.fsum(mp.exp(v-top) for v in vals))

def budget_at_z(L,a,z):
    t=(1-z)/(2*z*z)
    H=-mp.log(mp.mpf(1)/L)/L-(1-mp.mpf(1)/L)*mp.log(1-mp.mpf(1)/L)
    value=t*mp.log(a/(mp.e*z*z*(1+z)))+mp.log((1+z)/(2*z))
    return L*H+mp.log(5)+t*mp.log(4)+L*value

def mp_rate(L,a):
    a=mp.mpf(str(a))
    if mp.log(4)+L*mp.log(a/2)>=0:return -budget_at_z(L,a,mp.mpf(1))
    z=root(lambda z:z*z*(1+z)-a*mp.power(4,mp.mpf(1)/L),0,1)
    return -budget_at_z(L,a,z)

def main():
    finite=0; rational=0; tiny=0; rates=0; thresholds=0; adjacent=[]
    for txt in ['1e-10','1e-6','.001','.1','.2499','.25','.5','1','10','100']:
        t=mp.mpf(txt);value=variational(mp.mpf(2),t)
        assert abs(mp.mpf(g(2,float(t)))-value)<mp.mpf('2e-13')*(1+abs(value))
        for n in [1,2,3,5,10,20,100,500]:
            assert exact_log(n,mp.mpf(2),t*n)<=n*value+mp.mpf('1e-65')
            finite+=1
    for txt in ['1e-20','1e-16','1e-12','1e-8','0.00009999','0.0001','0.00010001']:
        t=mp.mpf(txt);value=variational(mp.mpf(2),t)
        assert abs(mp.mpf(g(2,float(t)))-value)<mp.mpf('1e-14')*t+mp.mpf('1e-65')
        tiny+=1
    for n in [1,2,3,5,10]:
        mu=[Fraction(1)]+[Fraction(math.prod(range(1,2*r,2)),2) for r in range(1,13)]
        moments=[Fraction(1)]+[Fraction(0)]*12
        for _ in range(n):
            moments=[sum(Fraction(math.comb(r,j))*moments[j]*mu[r-j] for j in range(r+1)) for r in range(13)]
        for r in range(1,13):
            exact=moments[r]*Fraction(2,n)**r
            assert abs(mp.log(exact.numerator)-mp.log(exact.denominator)-exact_log(n,mp.mpf(2),mp.mpf(r)))<mp.mpf('1e-70')
            rational+=1
    data=json.loads((DATA/'finite_budget.json').read_text())
    assert data['schema']=='optimized-finite-g-v1'
    for row in data['thresholds']:
        L=row['L']
        H=-mp.log(mp.mpf(1)/L)/L-(1-mp.mpf(1)/L)*mp.log(1-mp.mpf(1)/L)
        z=root(lambda z:(1-z)/(2*z*z)-mp.log((1+z)/(2*z))-H-mp.log(5)/L,'1e-30',1)
        a=z*z*(1+z)*mp.power(4,-mp.mpf(1)/L)
        assert abs(a-row['alpha_star'])<mp.mpf('2e-12')
        assert abs(budget_at_z(L,a,z))<mp.mpf('1e-65')
        thresholds+=1
    for row in data['rates']:
        for cell in row['values']:
            c=mp_rate(cell['L'],row['alpha'])
            assert abs(c-cell['c'])<mp.mpf('1e-8')
            assert cell['n_001']==(int(mp.ceil(mp.log(200)/c)) if c>0 else None)
            rates+=1
    for alpha,L in data['depths'].items():
        c0,c1=mp_rate(L-1,alpha),mp_rate(L,alpha)
        assert c0<=0<c1
        adjacent.append({'alpha':alpha,'depth':L,'previous_rate':mp.nstr(c0,30),'first_positive_rate':mp.nstr(c1,30)})
    # This deliberately incorrect scale must be rejected by a finite exact sum.
    assert exact_log(100,mp.mpf(2),mp.mpf(10))>100*variational(mp.mpf(1),mp.mpf('.1'))
    print(json.dumps({'status':'PASS','precision_digits':80,'finite_moment_checks':finite,
        'rational_convolution_checks':rational,'tiny_theta_checks':tiny,'threshold_rows':thresholds,
        'rate_cells_and_integer_widths':rates,'adjacent_depth_checks':adjacent,
        'scope':'Numerical verification of finite formulas; not proof certification or an independent publication vote.'},indent=2))

if __name__=='__main__':main()

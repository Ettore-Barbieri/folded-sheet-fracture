"""Williams eigenvalues for a Kirchhoff plate tearing crack whose tip touches a line.
Crack on theta=pi, line on theta=pi/2, antisymmetry (w odd) on theta=0.
Modes: continuous (flat plate, expect 1/2), support (w=0 on the line, slope and
moment continuous: the long-wave limit of the ridge), hinge (w=0, M=0 both sides).
w = r^(lambda+1) F(theta); moments ~ r^(lambda-1). nu set in the file."""
import sympy as sp, numpy as np
from scipy.optimize import brentq
r,th,l=sp.symbols('r theta lambda')
nu=sp.Rational(3,10)
def ops(w):
    lap=sp.diff(w,r,2)+sp.diff(w,r)/r+sp.diff(w,th,2)/r**2
    M=-(sp.diff(w,r)/r+sp.diff(w,th,2)/r**2+nu*sp.diff(w,r,2))
    Mrt=-(1-nu)*sp.diff(sp.diff(w,th)/r,r)
    V=-sp.diff(lap,th)/r+sp.diff(Mrt,r)
    return M,V
b=[sp.cos((l+1)*th),sp.sin((l+1)*th),sp.cos((l-1)*th),sp.sin((l-1)*th)]
odd=[b[1],b[3]]; h=sp.pi/2
def build(mode):
    R=[]
    for f in b:
        w=r**(l+1)*f; M,V=ops(w)
        col=[M.subs(th,sp.pi),V.subs(th,sp.pi)]
        if mode=='support': col+=[w.subs(th,h),0,sp.diff(w,th).subs(th,h),sp.diff(w,th,2).subs(th,h)]
        elif mode=='hinge': col+=[w.subs(th,h),0,sp.diff(w,th,2).subs(th,h),0]
        else: col+=[sp.diff(w,th,k).subs(th,h) for k in range(4)]
        R.append(col)
    for f in odd:
        w=r**(l+1)*f
        if mode=='support': col=[0,0,0,w.subs(th,h),-sp.diff(w,th).subs(th,h),-sp.diff(w,th,2).subs(th,h)]
        elif mode=='hinge': col=[0,0,0,w.subs(th,h),0,sp.diff(w,th,2).subs(th,h)]
        else: col=[0,0]+[-sp.diff(w,th,k).subs(th,h) for k in range(4)]
        R.append(col)
    return sp.lambdify(l,sp.Matrix(R).T.subs(r,1),'numpy')
for mode in ('continuous','support','hinge'):
    F=build(mode); d=lambda x: np.linalg.det(np.array(F(x),dtype=float))
    xs=np.linspace(0.005,2.2,2000); v=[d(x) for x in xs]
    print(mode,[round(brentq(d,xs[i],xs[i+1]),4) for i in range(len(xs)-1) if np.sign(v[i])!=np.sign(v[i+1])],flush=True)

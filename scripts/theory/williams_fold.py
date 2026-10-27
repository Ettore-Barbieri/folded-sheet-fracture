"""Williams eigenvalues (run: python3 williams_fold.py)

Results: Mode I with the hinge: lambda^2 + cos(pi lambda) = 0, lambda = 0.6298.
Mode II with the hinge: (1 + cos pi lambda)(lambda^2 + cos pi lambda - 1) = 0,
lambda = 1 (zero-stress mechanism) and 1.2798 (first mode with stress).

Williams eigenvalues for a crack whose tip sits on a hinge fold
(long-wave limit, lc << r << a). Membrane (Airy) problem in the half plane
y >= 0 by symmetry about the crack line; ridge at theta = pi/2.
Symmetric (Mode I) and antisymmetric (Mode II) families."""
import numpy as np
import sympy as sp
from scipy.optimize import brentq

lam, th = sp.symbols('lambda theta', real=True)
A, B, C, D = sp.symbols('A B C D')
p, m = lam + 1, lam - 1

# Region 1 (cracked quadrant): traction-free crack face at theta = pi
f1 = A*(sp.cos(p*(th-sp.pi)) - sp.cos(m*(th-sp.pi))) \
   + B*(m*sp.sin(p*(th-sp.pi)) - p*sp.sin(m*(th-sp.pi)))

def system(family):
    if family == 'I':      # symmetric: sigma_rtheta = 0, u_theta = 0 on theta = 0 -> even f
        f2 = C*sp.cos(p*th) + D*sp.cos(m*th)
    elif family == 'II':   # antisymmetric: sigma_thetatheta = 0, u_r = 0 on theta = 0 -> odd f
        f2 = C*sp.sin(p*th) + D*sp.sin(m*th)
    elif family == 'flat_I':
        return None
    h = sp.pi/2
    eqs = [f1.subs(th, h),                                   # sigma_xx = 0 on the ridge, side 1
           f2.subs(th, h),                                   # sigma_xx = 0 on the ridge, side 2
           (sp.diff(f1, th) - sp.diff(f2, th)).subs(th, h),  # shear continuous
           (sp.diff(f1, th, 2) - sp.diff(f2, th, 2)).subs(th, h)]  # u_r (along ridge) continuous
    M = sp.Matrix([[sp.diff(e, v) for v in (A, B, C, D)] for e in eqs])
    return sp.lambdify(lam, M.det(), 'numpy')

for fam in ('I', 'II'):
    F = system(fam)
    xs = np.linspace(0.01, 1.99, 4000)
    vals = np.array([F(x) for x in xs])
    roots = [brentq(F, xs[i], xs[i+1]) for i in range(len(xs)-1)
             if np.sign(vals[i]) != np.sign(vals[i+1])]
    print(fam, 'roots in (0,2):', [round(r, 5) for r in roots])
print('check lambda^2 + cos(pi lambda) = 0 root:', brentq(lambda l: l*l+np.cos(np.pi*l), 0.3, 0.9))

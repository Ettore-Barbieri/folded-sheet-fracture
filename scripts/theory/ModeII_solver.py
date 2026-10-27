"""Mode II (antisymmetric family) for a crack perpendicular to a rooftop fold.
Membrane glide dislocations b = d[u_x]/dt and plate dislocations g = d^2[w]/dt^2
(w odd in y: tearing). Lengths in units of lc. Kirchhoff plate.

Equations on the crack (-delta-2a < x < -delta), X=|x|, T=|t|, S=X+T, c=(1-nu)/2:
 membrane: 1 + (1/pi) int f_b [1/(t-x) - X(M1 + T M2)]
              + (4/pi) int f_g X[(1-c) M2 - c T M3] = 0
 plate:   -(c3/pi) int f_g/(t-x)
          + (1/(2 pi)) int f_b [(1-nu) X (M2 + T M3) - 2 (M1 + T M2)]
          + (2/pi) int f_g [(1-nu)(1-c) X M1 - (1-nu) c X T M2 - 2(1-c) M0 + 2 c T M1]
          = C0 (unknown constant; the plate condition V_y = 0 is imposed in integrated form)
 closure: int f_b = 0, int f_g = 0, int t f_g = 0.
 c3 = (1-nu)(3+nu).  G/G_flat = [phi_b^2 + 4 c3 phi_g^2] / a^2.
"""
import numpy as np
from StageB_solver import moments, cheb_end

def solve2(a, delta, nu=0.3, n=200, ridge=True, plate=True):
    i=np.arange(1,n+1); s=np.cos((2*i-1)*np.pi/(2*n))
    j=np.arange(1,n);   r=np.cos(j*np.pi/n)
    m=-delta-a; t=m+a*s; x=m+a*r
    X,Tt=np.meshgrid(x,t,indexing='ij')
    Xa=np.abs(X); Ta=np.abs(Tt); S=Xa+Ta
    M0,M1,M2,M3,M4=moments(S)
    c=(1-nu)/2; c3=(1-nu)*(3+nu)
    C=1.0/(Tt-X)
    w=np.pi/n
    z=0.0 if ridge else 0.0*X
    Kbb = -Xa*(M1+Ta*M2)
    Kbg = Xa*((1-c)*M2 - c*Ta*M3)
    Kgb = (1-nu)*Xa*(M2+Ta*M3) - 2*(M1+Ta*M2)
    Kgg = (1-nu)*(1-c)*Xa*M1 - (1-nu)*c*Xa*Ta*M2 - 2*(1-c)*M0 + 2*c*Ta*M1
    if not ridge:
        Kbb=Kbg=Kgb=Kgg=0*X
    A11=(w/np.pi)*(C+Kbb)
    if not plate:
        A=np.vstack([A11,np.ones((1,n))]); rhs=np.concatenate([-np.ones(n-1),[0]])
        pf=np.linalg.solve(A,rhs); pg=np.zeros(n)
    else:
        A12=(4*w/np.pi)*Kbg
        A21=(w/(2*np.pi))*Kgb
        A22=-(c3*w/np.pi)*C + (2*w/np.pi)*Kgg
        # unknowns: pf (n), pg (n), C0 (1)
        top=np.hstack([A11,A12,np.zeros((n-1,1))])
        mid=np.hstack([A21,A22,-np.ones((n-1,1))])
        c1=np.concatenate([np.ones(n),np.zeros(n),[0]])
        c2=np.concatenate([np.zeros(n),np.ones(n),[0]])
        c4=np.concatenate([np.zeros(n),t,[0]])
        A=np.vstack([top,mid,c1,c2,c4])
        rhs=np.concatenate([-np.ones(n-1),np.zeros(n-1),[0,0,0]])
        sol=np.linalg.solve(A,rhs); pf=sol[:n]; pg=sol[n:2*n]
    Pf=cheb_end(pf,s,1.0); Pg=cheb_end(pg,s,1.0)
    Gm=Pf**2/a**2; Gb=4*c3*Pg**2/a**2
    return Gm+Gb, Gm, Gb

def conv(a, d, tol=1e-4, **kw):
    n=200; g=solve2(a,d,n=n,**kw)
    while n<3200:
        n*=2; g2=solve2(a,d,n=n,**kw)
        if abs(g2[0]-g[0])<tol*abs(g2[0]): return g2
        g=g2
    return g

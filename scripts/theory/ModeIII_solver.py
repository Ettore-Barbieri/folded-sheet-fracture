"""Mode III (tearing) for a crack perpendicular to a rooftop fold.
Antisymmetric family: same kernels and unknowns as mode II (ModeII_solver).
Load: the two crack faces are torn apart by a uniform transverse shear force
V0 per unit length (V_y = -V0 on the faces); no remote membrane load.
In the integrated plate condition this gives a right-hand side linear in x.
Scaling: as mode II with S = V0/alpha, so G_ref = pi V0^2 a/(alpha^2 E h^2);
results are reported as G/G0, with G0 the same crack in a planar sheet
(computed with the ridge switched off)."""
import numpy as np
from StageB_solver import moments, cheb_end

def solve3(a, delta, nu=0.3, n=200, ridge=True, want_I=False):
    i=np.arange(1,n+1); s=np.cos((2*i-1)*np.pi/(2*n))
    j=np.arange(1,n);   r=np.cos(j*np.pi/n)
    m=-delta-a; t=m+a*s; x=m+a*r
    X,Tt=np.meshgrid(x,t,indexing='ij')
    Xa=np.abs(X); Ta=np.abs(Tt); S=Xa+Ta
    M0,M1,M2,M3,M4=moments(S)
    c=(1-nu)/2; c3=(1-nu)*(3+nu); C=1.0/(Tt-X); w=np.pi/n
    if ridge:
        Kbb=-Xa*(M1+Ta*M2); Kbg=Xa*((1-c)*M2-c*Ta*M3)
        Kgb=(1-nu)*Xa*(M2+Ta*M3)-2*(M1+Ta*M2)
        Kgg=(1-nu)*(1-c)*Xa*M1-(1-nu)*c*Xa*Ta*M2-2*(1-c)*M0+2*c*Ta*M1
    else:
        Kbb=Kbg=Kgb=Kgg=0*X
    A11=(w/np.pi)*(C+Kbb); A12=(4*w/np.pi)*Kbg
    A21=(w/(2*np.pi))*Kgb; A22=-(c3*w/np.pi)*C+(2*w/np.pi)*Kgg
    top=np.hstack([A11,A12,np.zeros((n-1,1))])
    mid=np.hstack([A21,A22,-np.ones((n-1,1))])
    c1=np.concatenate([np.ones(n),np.zeros(n),[0]])
    c2=np.concatenate([np.zeros(n),np.ones(n),[0]])
    c4=np.concatenate([np.zeros(n),t,[0]])
    A=np.vstack([top,mid,c1,c2,c4])
    rhs=np.concatenate([np.zeros(n-1),-x,[0,0,0]])   # membrane: no load; plate: Q = C0 - x
    sol=np.linalg.solve(A,rhs); pf=sol[:n]; pg=sol[n:2*n]
    Pf=cheb_end(pf,s,1.0); Pg=cheb_end(pg,s,1.0)
    Gm=Pf**2/a**2; Gt=4*c3*Pg**2/a**2
    out=(Gm+Gt,Gm,Gt)
    if want_I:   # work of the face load ~ int [w] dx = 1/2 int t^2 g dt
        out=out+(0.5*(np.pi/n)*np.sum(pg*t**2),)
    return out

def ratio(a,d,n=400,**kw):
    G=solve3(a,d,n=n,**kw); G0=solve3(a,d,n=n,ridge=False)[0]
    return G[0]/G0, G[1]/G0, G[2]/G0

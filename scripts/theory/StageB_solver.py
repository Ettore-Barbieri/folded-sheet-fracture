"""Stage B: crack perpendicular to a rooftop ridge, coupled membrane+bending
singular integral equations, Gauss-Chebyshev. Lengths in units of l_c."""
import numpy as np
from scipy import special

def fg(p):
    si,ci=special.sici(p); s2=si-np.pi/2
    return ci*np.sin(p)-s2*np.cos(p), -ci*np.cos(p)-s2*np.sin(p)

def moments(p):
    f,g=fg(p)
    return f, g, 1/p-f, 1/p**2-g, 2/p**3-1/p+f

def kernels(x,t,nu):
    X=np.abs(x); T=np.abs(t); S=X+T; c=(1-nu)/2
    M0,M1,M2,M3,M4=moments(S)
    Jbb = X*M2 - M1
    JNb = X*M2 - M1 + c*T*(X*M3 - M2)
    JMb = (1+nu)*M2 + (1-nu)*X*M3
    JMB = (1+nu)*M0 + (1-nu)*X*M1 + c*T*((1+nu)*M1 + (1-nu)*X*M2)
    return Jbb,JNb,JMb,JMB

def cheb_end(phi, s, end=1.0):
    # Chebyshev interpolation through nodes s_i, evaluate at end
    n=len(s); th=np.arccos(s)
    k=np.arange(n)
    T=np.cos(np.outer(th,k))
    coef=np.linalg.solve(T,phi)
    return np.sum(coef*np.cos(k*np.arccos(end)))

def solve(a, delta, nu=0.3, n=200, ridge=True, bending=True):
    """crack on (-delta-2a, -delta); returns G/G_flat at near tip (and parts)"""
    i=np.arange(1,n+1); s=np.cos((2*i-1)*np.pi/(2*n))
    j=np.arange(1,n);   r=np.cos(j*np.pi/n)
    m=-delta-a
    t=m+a*s; x=m+a*r
    X,Tt=np.meshgrid(x,t,indexing='ij')
    C=1.0/(Tt-X)                       # 1/(t-x)
    w=np.pi/n
    if ridge:
        Jbb,JNb,JMb,JMB=kernels(X,Tt,nu)
    else:
        Jbb=JNb=JMb=JMB=np.zeros_like(X)
    sg=np.sign(Tt)
    # unknowns [phi_f (n), phi_g (n)]
    A11=(w/np.pi)*(C + Tt*Jbb)
    A12=(4*w/np.pi)*sg*JNb
    A21=-(w/(2*np.pi))*Tt*JMb
    A22=-((1-nu)*(3+nu)*w/np.pi)*C + (2*w/np.pi)*sg*JMB
    if not bending:
        A=np.vstack([A11, np.ones((1,n))]); rhs=np.concatenate([-np.ones(n-1),[0]])
        pf=np.linalg.solve(A,rhs); pg=np.zeros(n)
    else:
        top=np.hstack([A11,A12]); bot=np.hstack([A21,A22])
        c1=np.concatenate([np.ones(n),np.zeros(n)])
        c2=np.concatenate([np.zeros(n),np.ones(n)])
        A=np.vstack([top,bot,c1,c2])
        rhs=np.concatenate([-np.ones(n-1),np.zeros(n-1),[0,0]])
        sol=np.linalg.solve(A,rhs); pf=sol[:n]; pg=sol[n:]
    # phi(s)=a*psi? f = phi/(a sqrt(1-s^2)) relative to t; phi as defined acts
    # in f(t)=phi/sqrt((t-c1)(c2-t)) -> our quadrature used phi directly
    Pf=cheb_end(pf,s,1.0); Pg=cheb_end(pg,s,1.0)
    Gm=Pf**2/a**2; Gb=4*(1-nu)*(3+nu)*Pg**2/a**2
    return Gm+Gb, Gm, Gb


def solve_converged(a, delta, nu=0.3, tol=1e-4, n0=200, nmax=3200, **kw):
    """solve() with the number of collocation points doubled until
    G changes by less than tol (relative). Large cracks (a/lc >= 1000)
    near the fold need n of order 1000 or more."""
    n = n0
    prev = solve(a, delta, nu=nu, n=n, **kw)
    while n < nmax:
        n *= 2
        cur = solve(a, delta, nu=nu, n=n, **kw)
        if abs(cur[0] - prev[0]) <= tol * abs(cur[0]):
            return cur
        prev = cur
    return prev

"""Crack-face profiles (jumps across the crack) from the dislocation densities.
Mode I : [u_y] (opening) and [w_,y] (face rotation)
Mode II: [u_x] (sliding) and [w] (tearing)
Mode III: [w] (tearing) and [u_x] (sliding)
Units: U0 = 4 F lc/(E h) (F = load scale). In-plane jumps in U0; out-of-plane:
lc*[w_,y] and [w] are (4/alpha)*U0 * J, reported here as alpha*[.]/(4 U0) = J."""
import numpy as np, sys
sys.path.insert(0,'..')
import StageB_solver as S, ModeII_solver as M2, ModeIII_solver as M3
def patch(mod, fn, ret):
    src=open(mod.__file__).read().replace(ret, ret.replace('return','__R=(') + ',pf,pg,s);return __R' if False else ret)
def nodes_solve(modname, a, d, n=800, ridge=True):
    mod={'I':S,'II':M2,'III':M3}[modname]
    src=open(mod.__file__).read()
    if modname=='I':
        src=src.replace('return Gm+Gb, Gm, Gb','return Gm+Gb, Gm, Gb, pf, pg, s',1); fn='solve'
    elif modname=='II':
        src=src.replace('return Gm+Gb, Gm, Gb','return Gm+Gb, Gm, Gb, pf, pg, s',1); fn='solve2'
    else:
        src=src.replace('out=(Gm+Gt,Gm,Gt)','out=(Gm+Gt,Gm,Gt,pf,pg,s)',1); fn='solve3'
    ns={'__name__':'x'}; import os; cwd=os.getcwd(); os.chdir(os.path.dirname(mod.__file__))
    exec(compile(src,mod.__file__,'exec'),ns); os.chdir(cwd)
    r=ns[fn](a,d,n=n,ridge=ridge)
    return r[3],r[4],r[5]
def cheb_coeffs(phi,s):
    n=len(s); th=np.arccos(s); k=np.arange(n)
    T=np.cos(np.outer(th,k)); return np.linalg.lstsq(T,phi,rcond=None)[0]
def jump(phi,s,sq):   # J(s)=int_{-1}^{s} phi/sqrt(1-s^2) ds' = -(sum_k c_k sin(k th)/k) + c0 (pi-th)
    c=cheb_coeffs(phi,s); th=np.arccos(sq); k=np.arange(1,len(c))
    return c[0]*(np.pi-th) - (np.sin(np.outer(th,k))/k)@c[1:]
def profiles(mode,a,d,n=800,ridge=True):
    pf,pg,s=nodes_solve(mode,a,d,n,ridge)
    sq=np.linspace(-1,1,801); x=-d-a+a*sq
    Jm=jump(pf,s,sq); Jp=jump(pg,s,sq)          # [u]/U0 and (alpha/4)*lc*[w_,*]/U0 (ell_c units)
    if mode in ('II','III'):                    # integrate [w_,x] once more for [w]
        W=np.concatenate([[0],np.cumsum(0.5*(Jp[1:]+Jp[:-1])*np.diff(x))]); Jp=W
    return x,Jm,Jp

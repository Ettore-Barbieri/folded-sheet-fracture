"""Mode I stress intensity factors with the tip on the fold.
K_I/K0 = membrane intensity (sigma_yy ~ K_I/sqrt(2 pi r)), K0 = sigma_inf sqrt(pi a).
K_b/K0 = bending intensity of the surface stress 6 M_yy/h^2 ~ K_b/sqrt(2 pi r).
Kirchhoff: G_b = (1+nu) K_b^2 / (3(3+nu) E)  ->  K_b/K0 = sqrt(3(3+nu)/(1+nu) * Gb/G0).
Signs: the overall sign of the solver is flipped so that K_I > 0."""
import numpy as np, StageB_solver as S
src=open(S.__file__).read().replace('return Gm+Gb, Gm, Gb','return Gm+Gb, Gm, Gb, Pf/a, Pg/a')
ns={}; exec(compile(src,'StageB_K','exec'),ns)
def K_at(a,delta=1e-3,nu=0.3,tol=1e-4):
    n=200; prev=None
    while True:
        r=ns['solve'](a,delta,nu=nu,n=n)
        if (prev and abs(r[0]-prev[0])<tol*r[0]) or n>=3200: break
        prev=r; n*=2
    G,Gm,Gb,pf,pg=r; s=np.sign(pf)
    KI=s*pf; Kb=s*np.sign(pg)*np.sqrt(3*(3+nu)/(1+nu)*Gb)
    return G,Gm,Gb,KI,Kb
if __name__=='__main__':
    rows=[]
    for a in (1,3,10,30,100,300,1000,3000):
        G,Gm,Gb,KI,Kb=K_at(a); rows.append((a,G,Gm,Gb,KI,Kb,KI+Kb,KI-Kb))
        print(a,[round(v,4) for v in rows[-1][1:]],flush=True)
    np.savetxt('../../results/theory/K_at_fold_vs_a.csv',np.array(rows),delimiter=',',fmt='%.6g',comments='',
      header='a_over_lc,G_over_G0,Gm,Gb,KI_over_K0,Kb_over_K0,Kplus,Kminus')

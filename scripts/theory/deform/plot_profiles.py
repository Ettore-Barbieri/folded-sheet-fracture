import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from profiles import profiles
a=100; ds=[(1e-3,'tip on the fold','C3','-'),(3,r'$\delta = 3\ell_c$','C1','--'),(100,r'$\delta = a$','C0','-.')]
names={'I':(r'opening $[u_y]$',r'face rotation $\frac{\alpha}{4}\,\ell_c[w_{,y}]$'),
       'II':(r'sliding $[u_x]$',r'tearing $\frac{\alpha}{4}[w]$'),
       'III':(r'sliding $[u_x]$',r'tearing $\frac{\alpha}{4}[w]$')}
fig,ax=plt.subplots(3,2,figsize=(10,10.5))
for i,m in enumerate(('I','II','III')):
    x0,m0,p0=profiles(m,a,1e-3,ridge=False)
    ref = np.abs(m0).max() if m!='III' else np.abs(p0).max()
    for d,lab,c,ls in ds:
        x,jm,jp=profiles(m,a,d)
        xi=(x-x[0])/(2*a)
        ax[i,0].plot(xi,jm/ref,c=c,ls=ls,lw=2,label=lab); ax[i,1].plot(xi,jp/ref,c=c,ls=ls,lw=2,label=lab)
    xi0=(x0-x0[0])/(2*a)
    ax[i,0].plot(xi0,m0/ref,'k:',lw=1.5,label='planar sheet'); ax[i,1].plot(xi0,p0/ref,'k:',lw=1.5,label='planar sheet')
    ax[i,0].set_title(f'mode {m}: '+names[m][0]); ax[i,1].set_title(f'mode {m}: '+names[m][1])
    for j in (0,1): ax[i,j].grid(alpha=.3); ax[i,j].axhline(0,c='k',lw=.5)
    ax[i,0].set_ylabel('jump / planar maximum')
ax[0,0].legend(fontsize=9); 
for j in (0,1): ax[2,j].set_xlabel(r'position along the crack, $(x-c_1)/2a$  (near tip at 1)')
fig.suptitle(rf'Crack-face jumps, $a/\ell_c = {a}$ (theory). Out-of-plane jumps scaled by $\alpha/4$', fontsize=12)
fig.tight_layout(); fig.savefig('crack_profiles_a100.png',dpi=130)

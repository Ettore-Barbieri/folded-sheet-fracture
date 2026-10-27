import numpy as np, matplotlib, sys; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.interpolate import griddata
import os
# COMSOL exports, relative to this script: <repo>/results/comsol/
U=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","results","comsol")+os.sep
SYM={'I':'sym','II':'anti','III':'anti'}
def load(m):
    P={l.split()[0]:float(l.split()[1]) for l in open(U+f'fig_Mode{m}_params.txt')}
    F={k:np.loadtxt(U+f'fig_Mode{m}_{k}.txt',comments='%') for k in ('u','v','w','Ws','mises')}
    return P,F
def remove_remote(P,xyz,vals):
    a=P['a_crack']; xc=-P['delta_tip']-a; r=np.hypot(xyz[:,0]-xc,xyz[:,1]); out=vals.copy()
    for side in (xyz[:,0]<0, xyz[:,0]>=0):
        far=side&(r>8*a)
        A=np.c_[np.ones(far.sum()),xyz[far,0],xyz[far,1]]; Aall=np.c_[np.ones(side.sum()),xyz[side,0],xyz[side,1]]
        out[side]-=Aall@np.linalg.lstsq(A,vals[far],rcond=None)[0]
    return out
def gridf(P,arr,X,Y,win):
    a=P['a_crack']; m=(arr[:,0]>win[0]*a*1.4)&(arr[:,0]<win[1]*a*1.4+1e-9)&(arr[:,1]<win[2]*a*1.4)
    return griddata(arr[m][:,:2],arr[m][:,3],(X,Y),method='linear')
def render(m,ax,field,win=(-2.8,1.4,1.4),amp_frac=0.22,cmap='YlOrRd',vmin=-2.5):
    P,F=load(m); a=P['a_crack']; al=P['alpha_g']
    X,Y=np.meshgrid(np.linspace(win[0]*a,win[1]*a,300),np.linspace(0,win[2]*a,150))
    Z=-al*np.abs(X)
    d={}
    for k in ('u','v','w'):
        arr=F[k].copy(); arr[:,3]=remove_remote(P,arr[:,:3],arr[:,3]); d[k]=gridf(P,arr,X,Y,win)
    c=gridf(P,F['Ws' if field=='plain' else field],X,Y,win)
    dmax=np.nanmax(np.sqrt(d['u']**2+d['v']**2+d['w']**2)); amp=amp_frac*a/dmax
    s=(1,-1,1) if SYM[m]=='sym' else (-1,1,-1)
    cn=np.log10(np.clip(c/np.nanmax(c),10**vmin,1)); cmapf=colors.LinearSegmentedColormap.from_list('trunc',plt.get_cmap('inferno')(np.linspace(0.28,1,256))); norm=colors.Normalize(vmin,0)
    allx=[];ally=[];allz=[]
    for sg,(su,sv,sw) in ((1,(1,1,1)),(-1,s)):
        xs=(X+amp*su*d['u'])/a; ys=(sg*Y+amp*sv*d['v'])/a; zs=(Z+amp*sw*d['w'])/a
        if field=='plain':
            from matplotlib.colors import LightSource
            ax.plot_surface(xs,ys,zs,color='#c9d3e0',rstride=2,cstride=2,linewidth=0,antialiased=False,shade=True,lightsource=LightSource(azdeg=315,altdeg=35))
        else:
            ax.plot_surface(xs,ys,zs,facecolors=cmapf(norm(np.nan_to_num(cn,nan=vmin))),rstride=2,cstride=2,linewidth=0,antialiased=False,shade=False)
        # deformed ridge (column nearest x = 0), drawn slightly above the surface
        j=np.argmin(np.abs(X[0]))
        ax.plot(xs[:,j],ys[:,j],zs[:,j]+0.05,color='#2b8cff',lw=2.2,ls='--',zorder=10)
        allx.append(xs);ally.append(ys);allz.append(zs)
        c1=(-P['delta_tip']-2*a)/a; c2=-P['delta_tip']/a; xr=X[0]/a
        kk=(xr>=c1)&(xr<=c2)
        ax.plot(xs[0,kk],ys[0,kk],zs[0,kk]+0.02,color='#d62728',lw=2.0,zorder=11)
        xo=np.r_[X[0],X[-1][::-1]]/a; yo=np.r_[sg*Y[0],sg*Y[-1][::-1]]/a
        pass
    X3=np.concatenate([q.ravel() for q in allx]);Y3=np.concatenate([q.ravel() for q in ally]);Z3=np.concatenate([q.ravel() for q in allz])
    lims=[(np.nanmin(q),np.nanmax(q)) for q in (X3,Y3,Z3)]
    ax.set_xlim(*lims[0]); ax.set_ylim(*lims[1]); ax.set_zlim(*lims[2])
    ax.set_box_aspect(tuple(l[1]-l[0] for l in lims),zoom=0.95); ax.view_init(elev=24,azim=-58)
    ax.set_axis_off(); pass
    print(m,'magnification',amp)
    return norm,cmapf
for field,lab in (('plain',''),('Ws','strain energy density / maximum')):
    fig=plt.figure(figsize=(15,4.6))
    for i,m in enumerate(('I','II','III')):
        ax=fig.add_subplot(1,3,i+1,projection='3d'); norm,cm=render(m,ax,field)
    if field=='plain':
        fig.subplots_adjust(left=0,right=1,top=0.92,bottom=0.05,wspace=-0.05); fig.savefig('rooftop3d_plain.png',dpi=300); continue
    cax=fig.add_axes([0.35,0.13,0.3,0.035]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cm),cax=cax,orientation='horizontal')
    cb.set_ticks([-2,-1,0]); cb.set_ticklabels(['$10^{-2}$','$10^{-1}$','1']); cb.set_label(lab)
    pass
    fig.subplots_adjust(left=0,right=1,top=0.92,bottom=0.2,wspace=-0.05)
    fig.savefig(f'rooftop3d_{field}.png',dpi=130)

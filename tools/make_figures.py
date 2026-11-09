import numpy as np, matplotlib, csv
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.rcParams.update({'font.size':8,'axes.titlesize':8,'axes.labelsize':8,
 'legend.fontsize':6.2,'xtick.labelsize':7,'ytick.labelsize':7,
 'lines.linewidth':1.2,'figure.dpi':150,'savefig.bbox':'tight'})
OUT='/home/user/GPS_spoofing_distribuito/paper/figs/'
C={'h0':'#777777','atk':'#c62828','blue':'#1f4e9c','grn':'#2e7d32','org':'#e8871a','pur':'#6a1b9a'}
R_E=6371000.0

# ======================= Fig. 1: scenario, ENU km, equal aspect ============
traj=np.loadtxt('scn4_traj.csv',delimiter=',')    # lon,lat
sp=np.loadtxt('scn4_spoofers.csv',delimiter=',')
dec=np.loadtxt('scn4_spoof.csv',delimiter=',')[0]
gt=np.loadtxt('scn4_grid_true.csv',delimiter=','); ge=np.loadtxt('scn4_grid_est.csv',delimiter=',')
lat0=float(traj[:,1].mean()); lon0=float(traj[:,0].mean())
def enu(a):
    a=np.atleast_2d(a)
    e=(a[:,0]-lon0)*np.pi/180*R_E*np.cos(lat0*np.pi/180)/1000
    n=(a[:,1]-lat0)*np.pi/180*R_E/1000
    return np.column_stack([e,n])
T=enu(traj); S=enu(sp); Dc=enu(dec)[0]; GT=enu(gt); GE=enu(ge)
fig,axs=plt.subplots(1,2,figsize=(3.5,2.1),gridspec_kw={'width_ratios':[1.45,1]})
ax=axs[0]
ax.plot(T[:,0],T[:,1],'k-',lw=1.2,label='true trajectory')
ax.plot(S[:,0],S[:,1],'^',color=C['atk'],ms=4,mfc='none',label='transmitters (24)')
ax.plot(Dc[0],Dc[1],'*',color=C['atk'],ms=9,label='decoy')
ctr=GT.mean(0); est=GE.mean(0)
ax.annotate('',xy=(est[0],est[1]),xytext=(ctr[0],ctr[1]),arrowprops=dict(arrowstyle='->',color=C['atk'],lw=1.0,ls='--'))
ax.plot(ctr[0],ctr[1],'o',color='k',ms=3,label='centre (true)')
ax.plot(est[0],est[1],'s',mfc='none',mec=C['atk'],ms=4,label='centre (computed)')
ax.set_xlabel('East [km]'); ax.set_ylabel('North [km]'); ax.set_aspect('equal'); ax.grid(alpha=.3)
ax.legend(loc='upper left',bbox_to_anchor=(-0.02,1.02),framealpha=.92,fontsize=5.4,handlelength=1.4,borderpad=.3,labelspacing=.25)
ax=axs[1]
gtt=GT-ctr+est
ax.scatter((gtt[:,0]-est[0])*1000,(gtt[:,1]-est[1])*1000,s=9,marker='o',facecolors='none',edgecolors='k',linewidths=.7,label='true, translated')
ax.scatter((GE[:,0]-est[0])*1000,(GE[:,1]-est[1])*1000,s=10,marker='s',facecolors='none',edgecolors=C['atk'],linewidths=.8,label='computed')
for p,col,ls in [(gtt,'k','--'),(GE,C['atk'],'-')]:
    x=(p[:,0]-est[0])*1000; y=(p[:,1]-est[1])*1000
    ax.plot([x.min(),x.max(),x.max(),x.min(),x.min()],[y.min(),y.min(),y.max(),y.max(),y.min()],ls,color=col,lw=.8)
ax.set_aspect('equal'); ax.grid(alpha=.3); ax.set_xlabel('East [m]'); ax.set_ylabel('North [m]')
ax.legend(loc='upper center',bbox_to_anchor=(0.5,1.28),fontsize=5.6,ncol=1,framealpha=.92,handlelength=1.2,borderpad=.25,labelspacing=.2)
fig.tight_layout(pad=.35); fig.savefig(OUT+'fig_scenario.pdf')

# ======================= Fig. 2: time series ==============================
ts={k:np.loadtxt(f'ts3_{k}.csv',delimiter=',') for k in ('h0','power')}
fig,axs=plt.subplots(2,1,figsize=(3.4,3.0),sharex=True)
axs[0].plot(ts['power'][:,0]/60,ts['power'][:,1]/1000,color=C['atk'])
axs[0].set_ylabel('centre displacement [km]'); axs[0].grid(alpha=.3); axs[0].set_ylim(bottom=0)
axs[1].plot(ts['h0'][:,0]/60,100*ts['h0'][:,2],'--',color=C['h0'],label='no attack')
axs[1].plot(ts['power'][:,0]/60,100*ts['power'][:,2],color=C['atk'],label='under attack')
axs[1].set_ylabel('distortion $E$ [%]'); axs[1].set_xlabel('time [min]'); axs[1].grid(alpha=.3); axs[1].legend(ncol=2)
fig.savefig(OUT+'fig_timeseries.pdf')

# ======================= Fig. 3: distortion vs spacing ====================
er=np.loadtxt('err4_runs.csv',delimiter=','); gs=sorted(set(er[:,0]))
fig,ax=plt.subplots(figsize=(3.2,2.3))
for col,ci,lab,ls,mk in [(C['h0'],2,'no attack (noise floor)','--','s'),(C['atk'],3,'under attack','-','o')]:
    med=[np.median(er[er[:,0]==g][:,ci]) for g in gs]
    lo=[np.percentile(er[er[:,0]==g][:,ci],10) for g in gs]; hi=[np.percentile(er[er[:,0]==g][:,ci],90) for g in gs]
    ax.fill_between(gs,100*np.array(lo),100*np.array(hi),color=col,alpha=.18,lw=0)
    ax.plot(gs,100*np.array(med),ls+mk,color=col,ms=4,label=lab)
ax.set_xlabel('formation spacing $g$ [m]'); ax.set_ylabel('distortion $E$ [%]'); ax.set_xticks(gs)
ax.grid(alpha=.3); ax.legend(); fig.savefig(OUT+'fig_distortion.pdf')

# ======================= Fig. 4: Pd vs spacing, 6 panels ==================
det=['Jansen min','Jansen max','Jansen majority','Worst-pair','Two-sided majority','Two-sided majority, 1-min']
cols=[1,2,3,4,5,7]
pd=np.loadtxt('pd4_ci.csv',delimiter=',')
fig,axs=plt.subplots(2,3,figsize=(3.5,2.5),sharex=True,sharey=True)
for a,(k,lab) in enumerate(zip(cols,det)):
    ax=axs.flat[a]; sub=pd[pd[:,0]==k]
    ax.fill_between(sub[:,1],sub[:,3],sub[:,4],color=C['atk'],alpha=.2,lw=0)
    ax.plot(sub[:,1],sub[:,2],'-o',color=C['atk'],ms=2.8)
    ax.axhline(0.01,ls='--',color='k',lw=.7)
    ax.set_title(lab,fontsize=6.3,pad=2); ax.grid(alpha=.3); ax.set_ylim(-.04,1.04); ax.set_xticks([5,10,15,20])
for ax in axs[1]: ax.set_xlabel('$g$ [m]')
for ax in axs[:,0]: ax.set_ylabel('$P_D$')
fig.tight_layout(pad=.3); fig.savefig(OUT+'fig_pd.pdf')

# ======================= Fig. 5: ROC, log x-axis ==========================
alarm_small=[True,True,True,False,False,False]
fig,axs=plt.subplots(2,3,figsize=(3.5,2.7),sharex=True,sharey=True)
for a,(k,lab,small) in enumerate(zip(cols,det,alarm_small)):
    ax=axs.flat[a]
    t0=np.loadtxt(f'roc4_g10_k{k}_s1.csv'); tk=np.loadtxt(f'roc4_g10_k{k}_s2.csv')
    t0=t0[~np.isnan(t0)]; tk=tk[~np.isnan(tk)]
    thr=np.unique(np.concatenate(([-1e9],np.sort(t0),[1e9])))
    if small: pfa=np.array([(t0<x).mean() for x in thr]); pdv=np.array([(tk<x).mean() for x in thr])
    else:     pfa=np.array([(t0>x).mean() for x in thr]); pdv=np.array([(tk>x).mean() for x in thr])
    o=np.argsort(pfa); pfa=np.clip(pfa[o],1e-4,1); pdv=pdv[o]
    ax.plot(pfa,pdv,color=C['atk']); ax.plot([1e-4,1],[1e-4,1],'k--',lw=.8)
    ax.set_xscale('log'); ax.set_xlim(1e-4,1); ax.set_ylim(-.04,1.04)
    ax.set_xticks([1e-4,1e-2,1]); ax.set_title(lab,fontsize=6.3,pad=2); ax.grid(alpha=.3,which='both')
for ax in axs[1]: ax.set_xlabel('$P_{FA}$')
for ax in axs[:,0]: ax.set_ylabel('$P_D$')
fig.tight_layout(pad=.3); fig.savefig(OUT+'fig_roc.pdf')

# ======================= Fig. 6: model validation =========================
v=np.loadtxt('valid_g20.csv',delimiter=',')   # t, E_meas, E_pred3D, E_simple, dil, E_pred_horiz
fig,axs=plt.subplots(1,2,figsize=(3.5,1.95))
ax=axs[0]
ax.scatter(100*v[:,2],100*v[:,1],s=6,color=C['atk'],alpha=.6,label='first-order model')
ax.scatter(100*v[:,3],100*v[:,1],s=6,color=C['blue'],alpha=.5,marker='x',label=r'crude scale $\kappa_0$')
lim=[0,max(100*v[:,1].max(),100*v[:,3].max())*1.05]
ax.plot(lim,lim,'k--',lw=.8); ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal')
ax.set_xlabel('predicted $E$ [%]'); ax.set_ylabel('measured $E$ [%]'); ax.grid(alpha=.3); ax.legend(loc='lower right')
ax=axs[1]
ax.plot(v[:,0]/60,100*v[:,1],color=C['atk'],label='measured')
ax.plot(v[:,0]/60,100*v[:,2],'k--',lw=.9,label='first-order model')
ax.set_xlabel('time [min]'); ax.set_ylabel('$E$ [%]'); ax.grid(alpha=.3); ax.legend(loc='upper right')
fig.tight_layout(pad=.3); fig.savefig(OUT+'fig_validation.pdf')

# ======================= Fig. 7: altitude sweep, 3 panels =================
cfgs=['ground','alt500','alt1000','alt2000','alt3500','alt5000','alt7000','alt10000']
cells=['g10','alt500','alt1000','alt2000','alt3500','alt5000','alt7000','alt10000']
hk=np.array([0.1,0.5,1,2,3.5,5,7,10])
th=[];dmed=[];dlo=[];dhi=[];Emed=[];Epred=[]
for c in cfgs:
    rows=list(csv.DictReader(open(f'geom_{c}.csv')))
    d=np.array([float(r['dil_signed']) for r in rows])
    th.append(np.median([float(r['theta_k_med']) for r in rows])); dmed.append(np.median(d))
    dlo.append(np.percentile(d,10)); dhi.append(np.percentile(d,90))
    Emed.append(np.median([float(r['E_meas']) for r in rows])); Epred.append(np.median([float(r['E_pred']) for r in rows]))
th=np.array(th);dmed=np.array(dmed);dlo=np.array(dlo);dhi=np.array(dhi);Emed=np.array(Emed);Epred=np.array(Epred)
sw={}
for r in csv.DictReader(open('sweep5.csv')): sw.setdefault(r['cell'],{})[int(r['k'])]=r
fig,axs=plt.subplots(3,1,figsize=(3.4,6.2))
ax=axs[0]
ax.fill_between(th,100*dlo,100*dhi,color=C['blue'],alpha=.12,lw=0,label='$D$, 10th–90th pct. over missions')
ax.plot(th,100*dmed,'-s',color=C['blue'],ms=3.5,label='$D$, median over missions')
ax.plot(th,100*Emed,'o',color=C['atk'],ms=3.5,label='$E$, measured')
ax.plot(th,100*Epred,'k--',lw=.9,label='$E$, first-order model')
ax.axhline(0,color='k',lw=.5); ax.set_ylim(-32,95); ax.set_xlim(-5,52)
ax.set_ylabel('[% of baseline], noiseless'); ax.grid(alpha=.3); ax.legend(loc='upper right',framealpha=.95)
ax2=ax.twiny(); ax2.set_xlim(ax.get_xlim()); ax2.set_xticks(th[[0,2,3,4,5,6,7]])
ax2.set_xticklabels([f'{x:g}' for x in hk[[0,2,3,4,5,6,7]]],fontsize=6); ax2.set_xlabel('transmitter altitude [km]',fontsize=7)
ax=axs[1]
for k,lab,col,mk in [(3,'Jansen majority',C['h0'],'^'),(4,'worst-pair',C['org'],'v'),(6,'worst-pair, 1-min',C['pur'],'d'),(5,'two-sided majority',C['atk'],'o'),(7,'two-sided maj., 1-min',C['grn'],'s')]:
    pdv=np.array([float(sw[c][k]['pd']) for c in cells]); lo=np.array([float(sw[c][k]['lo']) for c in cells]); hi=np.array([float(sw[c][k]['hi']) for c in cells])
    ax.fill_between(th,lo,hi,color=col,alpha=.12,lw=0); ax.plot(th,pdv,'-'+mk,color=col,ms=3,label=lab)
ax.axhline(0.01,ls='--',color='k',lw=.7); ax.set_ylim(-.04,1.12); ax.set_xlim(-5,52); ax.grid(alpha=.3)
ax.set_xlabel(r'transmitter elevation $\theta_k$ [deg]'); ax.set_ylabel('$P_D$ at $P_{FA}=10^{-2}$')
ax.legend(loc='upper center',ncol=2,framealpha=.95)
geo={}
allc=cfgs+['loit_circle','loit_los2000','loit_los5000','loit_los10000']
for c in allc:
    for r in csv.DictReader(open(f'geom_{c}.csv')): geo[('g10' if c=='ground' else c,int(r['mission']))]=float(r['dil_signed'])
P=[]
for r in csv.DictReader(open('permission_pd.csv')):
    k=(r['cell'],int(r['mission']))
    if k in geo: P.append((geo[k],float(r['pd5']),float(r['pd3'])))
P=np.array(P)
ax=axs[2]
ax.axvspan(-8,5,color='k',alpha=.07,lw=0)
ax.scatter(100*P[:,0],P[:,1],s=6,color=C['atk'],alpha=.3)
ax.scatter(100*P[:,0],P[:,2],s=6,color=C['h0'],alpha=.3,marker='^')
edges=np.array([-35,-20,-12,-8,-5,-3,-1.5,0,1.5,3,5,8,12,20,35,60,100])/100
xc=[];y5=[];y3=[]
for lo,hi in zip(edges[:-1],edges[1:]):
    m=(P[:,0]>=lo)&(P[:,0]<hi)
    if m.sum()>=3: xc.append(100*P[m,0].mean()); y5.append(P[m,1].mean()); y3.append(P[m,2].mean())
ax.plot(xc,y5,'-',color=C['atk'],lw=1.6,label='two-sided majority')
ax.plot(xc,y3,'-',color=C['h0'],lw=1.6,label='Jansen majority')
ax.annotate('blind window',xy=(-1.5,0.45),xytext=(-30,0.62),fontsize=6,
            arrowprops=dict(arrowstyle='->',lw=.7,color='k'))
ax.axhline(0.01,ls='--',color='k',lw=.7)
ax.set_xscale('symlog',linthresh=10); ax.set_xticks([-30,-10,-3,0,3,10,30,100])
ax.set_xticklabels(['-30','-10','-3','0','3','10','30','100'])
ax.set_xlim(-45,115); ax.set_ylim(-.04,1.04); ax.grid(alpha=.3)
ax.set_xlabel('signed dilation $D$ of the mission [%], noiseless'); ax.set_ylabel('$P_D$ of the mission')
ax.legend(loc='upper center',ncol=2,framealpha=.95)
fig.tight_layout(pad=.4); fig.savefig(OUT+'fig_sweep.pdf')
print('figures 1-7 written')

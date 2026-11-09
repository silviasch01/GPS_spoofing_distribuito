import numpy as np, matplotlib, csv
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':8,'axes.titlesize':8,'axes.labelsize':8,
 'legend.fontsize':6.4,'xtick.labelsize':7,'ytick.labelsize':7,
 'lines.linewidth':1.2,'figure.dpi':150,'savefig.bbox':'tight'})
OUT='/home/user/GPS_spoofing_distribuito/paper/figs/'
C={'h0':'#777777','atk':'#c62828','blue':'#1f4e9c','grn':'#2e7d32','org':'#e8871a','pur':'#6a1b9a'}
R_E=6371000.0

# ---------------- Fig. 1: scenario, lat/lon + formation zoom in metres ----
traj=np.loadtxt('scn4_traj.csv',delimiter=',')   # lon,lat
sp=np.loadtxt('scn4_spoofers.csv',delimiter=',')
dec=np.loadtxt('scn4_spoof.csv',delimiter=',')[0]
gt=np.loadtxt('scn4_grid_true.csv',delimiter=','); ge=np.loadtxt('scn4_grid_est.csv',delimiter=',')
lat0=float(gt[:,1].mean()); lon0=float(gt[:,0].mean())
fig,axs=plt.subplots(1,2,figsize=(3.5,2.05),gridspec_kw={'width_ratios':[1.5,1]})
ax=axs[0]
ax.plot(traj[:,0],traj[:,1],'k-',lw=1.2,label='true trajectory')
ax.plot(sp[:,0],sp[:,1],'^',color=C['atk'],ms=4,mfc='none',label='transmitters (24)')
ax.plot(dec[0],dec[1],'*',color=C['atk'],ms=9,label='decoy')
ctr=gt.mean(0); est=ge.mean(0)
ax.annotate('',xy=est,xytext=ctr,arrowprops=dict(arrowstyle='->',color=C['atk'],lw=1.0,ls='--'))
ax.plot(*ctr,'o',color='k',ms=3,label='centre (true)')
ax.plot(*est,'s',mfc='none',mec=C['atk'],ms=4,label='centre (computed)')
ax.set_xlabel('longitude [deg]'); ax.set_ylabel('latitude [deg]')
ax.set_aspect(1/np.cos(np.deg2rad(lat0))); ax.grid(alpha=.3)
ax.legend(loc='upper left',bbox_to_anchor=(-0.02,1.02),framealpha=.92,fontsize=5.4,
          handlelength=1.3,borderpad=.3,labelspacing=.25)
ax=axs[1]
def to_m(a,c):
    a=np.atleast_2d(a)
    return np.column_stack([(a[:,0]-c[0])*np.pi/180*R_E*np.cos(np.deg2rad(lat0)),
                            (a[:,1]-c[1])*np.pi/180*R_E])
GTm=to_m(gt,ctr); GEm=to_m(ge,est)
ax.scatter(GTm[:,0],GTm[:,1],s=9,marker='o',facecolors='none',edgecolors='k',linewidths=.7,label='true')
ax.scatter(GEm[:,0],GEm[:,1],s=10,marker='s',facecolors='none',edgecolors=C['atk'],linewidths=.8,label='computed')
for p,col,ls in [(GTm,'k','--'),(GEm,C['atk'],'-')]:
    ax.plot([p[:,0].min(),p[:,0].max(),p[:,0].max(),p[:,0].min(),p[:,0].min()],
            [p[:,1].min(),p[:,1].min(),p[:,1].max(),p[:,1].max(),p[:,1].min()],ls,color=col,lw=.8)
ax.set_aspect('equal'); ax.grid(alpha=.3)
ax.set_xlabel('East [m] from centre'); ax.set_ylabel('North [m]')
ax.legend(loc='upper center',bbox_to_anchor=(0.5,1.26),fontsize=5.8,framealpha=.92,
          handlelength=1.1,borderpad=.25,labelspacing=.2,ncol=2,columnspacing=.8)
fig.tight_layout(pad=.35); fig.savefig(OUT+'fig_scenario.pdf')

# ---------------- Fig. 2: distortion over one mission (single panel) ------
ts={k:np.loadtxt(f'ts3_{k}.csv',delimiter=',') for k in ('h0','power')}
fig,ax=plt.subplots(figsize=(3.3,1.9))
ax.plot(ts['h0'][:,0]/60,100*ts['h0'][:,2],'--',color=C['h0'],label='no attack')
ax.plot(ts['power'][:,0]/60,100*ts['power'][:,2],color=C['atk'],label='under attack')
ax.set_xlabel('time [min]'); ax.set_ylabel('distortion $E$ [%]')
ax.grid(alpha=.3); ax.legend(ncol=2); ax.set_ylim(bottom=0)
fig.tight_layout(pad=.3); fig.savefig(OUT+'fig_timeseries.pdf')

# ---------------- Fig. 3: model validation --------------------------------
v=np.loadtxt('valid_g20.csv',delimiter=',')
fig,axs=plt.subplots(1,2,figsize=(3.5,1.9))
ax=axs[0]
ax.scatter(100*v[:,2],100*v[:,1],s=6,color=C['atk'],alpha=.6,label='first-order model')
ax.scatter(100*v[:,3],100*v[:,1],s=6,color=C['blue'],alpha=.5,marker='x',label=r'crude scale $\kappa_0$')
lim=[0,max(100*v[:,1].max(),100*v[:,3].max())*1.05]
ax.plot(lim,lim,'k--',lw=.8); ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal')
ax.set_xlabel('predicted $E$ [%]'); ax.set_ylabel('measured $E$ [%]'); ax.grid(alpha=.3)
ax.legend(loc='lower right')
ax=axs[1]
ax.plot(v[:,0]/60,100*v[:,1],color=C['atk'],label='measured')
ax.plot(v[:,0]/60,100*v[:,2],'k--',lw=.9,label='first-order model')
ax.set_xlabel('time [min]'); ax.set_ylabel('$E$ [%]'); ax.grid(alpha=.3); ax.legend(loc='upper right')
fig.tight_layout(pad=.3); fig.savefig(OUT+'fig_validation.pdf')

# ---------------- Fig. 4: ROC at 10 m, the two majority tests -------------
fig,ax=plt.subplots(figsize=(3.0,2.2))
for k,lab,col,small in [(3,'Jansen majority (collapse)',C['h0'],True),
                        (5,'two-sided majority',C['atk'],False)]:
    t0=np.loadtxt(f'roc4_g10_k{k}_s1.csv'); tk=np.loadtxt(f'roc4_g10_k{k}_s2.csv')
    t0=t0[~np.isnan(t0)]; tk=tk[~np.isnan(tk)]
    thr=np.unique(np.concatenate(([-1e9],np.sort(t0),[1e9])))
    if small: pfa=np.array([(t0<x).mean() for x in thr]); pdv=np.array([(tk<x).mean() for x in thr])
    else:     pfa=np.array([(t0>x).mean() for x in thr]); pdv=np.array([(tk>x).mean() for x in thr])
    o=np.argsort(pfa); ax.plot(np.clip(pfa[o],1e-4,1),pdv[o],color=col,label=lab)
ax.plot([1e-4,1],[1e-4,1],'k--',lw=.8,label='chance')
ax.set_xscale('log'); ax.set_xlim(1e-4,1); ax.set_ylim(-.03,1.03)
ax.set_xlabel('$P_{FA}$ per epoch'); ax.set_ylabel('$P_D$')
ax.grid(alpha=.3,which='both'); ax.legend(loc='upper left')
fig.tight_layout(pad=.3); fig.savefig(OUT+'fig_roc.pdf')

# ---------------- Fig. 5: altitude sweep, two panels ----------------------
cfgs=['ground','alt500','alt1000','alt2000','alt3500','alt5000','alt7000','alt10000']
cells=['g10','alt500','alt1000','alt2000','alt3500','alt5000','alt7000','alt10000']
hk=np.array([0.1,0.5,1,2,3.5,5,7,10])
th=[];dmed=[];dlo=[];dhi=[]
for c in cfgs:
    rows=list(csv.DictReader(open(f'geom_{c}.csv')))
    d=np.array([float(r['dil_signed']) for r in rows])
    th.append(np.median([float(r['theta_k_med']) for r in rows]))
    dmed.append(np.median(d)); dlo.append(np.percentile(d,10)); dhi.append(np.percentile(d,90))
th=np.array(th);dmed=np.array(dmed);dlo=np.array(dlo);dhi=np.array(dhi)
sw={}
for r in csv.DictReader(open('sweep5.csv')): sw.setdefault(r['cell'],{})[int(r['k'])]=r
fig,axs=plt.subplots(2,1,figsize=(3.4,3.9),sharex=True)
ax=axs[0]
ax.axhspan(-8,5,color='k',alpha=.10,lw=0)
ax.text(49,-1.5,'blind window\nof the two-sided test',fontsize=5.8,ha='right',va='center',color='#444444')
ax.fill_between(th,100*dlo,100*dhi,color=C['blue'],alpha=.15,lw=0,label='10th–90th pct. over missions')
ax.plot(th,100*dmed,'-s',color=C['blue'],ms=3.5,label='median over missions')
ax.axhline(0,color='k',lw=.6)
ax.set_ylim(-32,95); ax.set_xlim(-4,51)
ax.set_ylabel('signed dilation $D$ [%]'); ax.grid(alpha=.3); ax.legend(loc='upper right',framealpha=.95)
ax2=ax.twiny(); ax2.set_xlim(ax.get_xlim()); ax2.set_xticks(th[[0,2,3,4,5,6,7]])
ax2.set_xticklabels([f'{x:g}' for x in hk[[0,2,3,4,5,6,7]]],fontsize=6)
ax2.set_xlabel('transmitter altitude [km]',fontsize=7)
ax=axs[1]
for k,lab,col,mk in [(3,'Jansen majority (collapse)',C['h0'],'^'),
                     (4,'worst-pair',C['org'],'v'),
                     (5,'two-sided majority',C['atk'],'o')]:
    pdv=np.array([float(sw[c][k]['pd']) for c in cells])
    lo=np.array([float(sw[c][k]['lo']) for c in cells]); hi=np.array([float(sw[c][k]['hi']) for c in cells])
    ax.fill_between(th,lo,hi,color=col,alpha=.13,lw=0); ax.plot(th,pdv,'-'+mk,color=col,ms=3.2,label=lab)
ax.axhline(0.01,ls='--',color='k',lw=.7)
ax.set_ylim(-.04,1.08); ax.set_xlim(-4,51); ax.grid(alpha=.3)
ax.set_xlabel(r'transmitter elevation $\theta_k$ seen from the swarm [deg]')
ax.set_ylabel('$P_D$ at $P_{FA}=10^{-2}$')
ax.legend(loc='upper left',framealpha=.95)
fig.tight_layout(pad=.35); fig.savefig(OUT+'fig_sweep.pdf')
print('slim figure set written (5 figures)')

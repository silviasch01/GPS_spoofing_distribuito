import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':8,'figure.dpi':150,'savefig.bbox':'tight'})
C={'atk':'#c62828','blue':'#1f4e9c','gry':'#555555'}
fig,ax=plt.subplots(figsize=(3.4,2.05))
ax.plot([-0.6,5.0],[0,0],color='#aaaaaa',lw=.9)
ax.fill_between([-0.6,5.0],[-0.3,-0.3],[0,0],color='#f0f0f0',lw=0)
pc=np.array([1.15,0.95]); du=np.array([1.05,0.0]); pu=pc+du
el=np.deg2rad(45); es=np.array([np.cos(el),np.sin(el)])
q=np.array([4.6,0.10]); vk=(q-pc)/np.linalg.norm(q-pc)
# unit vectors from the centre
ax.annotate('',xy=pc+1.45*es,xytext=pc,arrowprops=dict(arrowstyle='-|>',color=C['blue'],lw=1.5))
ax.annotate('',xy=pc+1.45*vk,xytext=pc,arrowprops=dict(arrowstyle='-|>',color=C['atk'],lw=1.5))
ax.text(*(pc+1.52*es+np.array([-0.30,0.10])),r'$\hat{\mathbf{e}}_s$',color=C['blue'],fontsize=10)
ax.text(*(pc+1.52*es+np.array([0.02,0.30])),'to satellite $s$',color=C['blue'],fontsize=6.5,ha='center')
ax.text(*(pc+1.50*vk+np.array([-0.10,-0.28])),r'$\hat{\mathbf{v}}_k$',color=C['atk'],fontsize=10)
ax.plot(*q,'^',color=C['atk'],ms=8,mfc='none',mew=1.3)
ax.text(q[0],q[1]-0.34,'transmitter $k$',color=C['atk'],fontsize=6.5,ha='center')
# swarm centre, displaced UAV
ax.annotate('',xy=pu,xytext=pc,arrowprops=dict(arrowstyle='-|>',color='k',lw=1.2))
ax.plot(*pc,'o',color='k',ms=5.5); ax.plot(*pu,'o',color='k',ms=4.5,mfc='white',mew=1.1)
ax.text(pc[0]-0.12,pc[1]+0.04,r'$\mathbf{p}^{\mathrm{ctr}}$',fontsize=8,ha='right',va='center')
ax.text(pc[0]+du[0]/2,pc[1]+0.14,r'$\boldsymbol{\delta}_u$',fontsize=10,ha='center')
ax.text(pu[0]+0.12,pu[1]+0.16,r'UAV $u$',fontsize=7,va='center')
# elevation angles
def arc(r,a0,a1,col,lab,off):
    t=np.linspace(a0,a1,60); ax.plot(pc[0]+r*np.cos(t),pc[1]+r*np.sin(t),color=col,lw=.9)
    am=(a0+a1)/2; p=pc+(r+off[0])*np.array([np.cos(am),np.sin(am)])+np.array([0,off[1]])
    ax.text(p[0],p[1],lab,color=col,fontsize=8.5,ha='center',va='center')
arc(0.72,0.30,el,C['blue'],r'$\mathrm{el}_s$',(0.26,0.0))
arc(1.18,np.arctan2(vk[1],vk[0]),0,C['atk'],r'$\theta_k$',(0.16,-0.04))
ax.text(2.2,-0.78,r'$\epsilon_{u,s}=(\hat{\mathbf{e}}_s-\hat{\mathbf{v}}_k)\cdot\boldsymbol{\delta}_u$'
        '\n'r'$\approx-(\hat{\mathbf{a}}_s\cdot\boldsymbol{\delta}_u)(\cos\theta_k-\cos\mathrm{el}_s)<0$',
        fontsize=7.2,ha='center',va='top',color=C['gry'])
ax.set_xlim(-0.6,5.2); ax.set_ylim(-1.55,2.45); ax.set_aspect('equal'); ax.axis('off')
fig.tight_layout(pad=.15); fig.savefig('/home/user/GPS_spoofing_distribuito/paper/figs/fig_geometry.pdf')
print('fig_geometry written')

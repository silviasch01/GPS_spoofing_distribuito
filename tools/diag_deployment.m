1;
% run from the repository root (needs physconst, e.g. Phased Array / Antenna Toolbox)
function D = pd(C, Nu)
  R=6371000; P=zeros(Nu,3);
  for u=1:Nu
    lat=C(3*u-2)*pi/180; lon=C(3*u-1)*pi/180; a=C(3*u);
    P(u,:)=(R+a)*[cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)];
  end
  D=[]; for i=1:Nu, for j=i+1:Nu, D(end+1)=norm(P(i,:)-P(j,:)); end, end
end
R=6371000;
for Rs=[2500 20000 50000]
  rng(11); o=GPS_spoofing_sim_20250617(struct('PLOT',false,'grid_step',20,'sigma_pr',3,'max_spoofed_sats',inf,'R_spoofers',Rs));
  Nt=size(o.coords_est,1); err=zeros(Nt,1); dist=zeros(Nt,1);
  cc=o.center_coords;
  for t=1:Nt
    De=pd(o.coords_est(t,:),25); Dt=pd(o.coords_grid(t,2:end),25); err(t)=median(abs(De-Dt)./Dt);
    c=o.coords(t,2:3); dist(t)=R*norm([(c(1)-cc(1))*pi/180,(c(2)-cc(2))*pi/180*cos(cc(1)*pi/180)]);
  end
  inside=dist<Rs;
  printf('R=%5d m | swarm inside circle %3.0f%% of epochs | dist to circle centre: min %5.0f max %5.0f m | distortion: overall %5.1f%%, inside %5.1f%%, outside %5.1f%% | corr(dist,err)=%.2f\n', ...
    Rs, 100*mean(inside), min(dist), max(dist), 100*median(err), 100*median(err(inside)), 100*median(err(~inside)), corr(dist,err));
end

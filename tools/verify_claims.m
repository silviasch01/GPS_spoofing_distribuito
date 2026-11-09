1;
% Verifications requested by the reviewer:
%  (1) D as in eq. (11): median over epochs of the per-epoch median over pairs
%  (2) M as in eq. (10): [(G'G)^-1 G' B]_{1:3}
%  (3) sign of the clock term in eq. (4)
%  (4) what happens when a satellite is re-routed to another transmitter
%  (5) depression-angle claim in III-C (|cos theta_k| = cos el_s)
%  (6) 2.5 km array: why the deformation inside is smaller
% run from the repository root
function p = frame(t,lat,lon,alt)
  R_E=6371000; T_E=86164; th0=0.23*2*pi; la=lat*pi/180; lo=lon*pi/180;
  p=(R_E+alt)*[cos(la)*cos(2*pi*t/T_E+th0+lo); cos(la)*sin(2*pi*t/T_E+th0+lo); sin(la)];
end
function G = enu(p)
  u1=p/norm(p); uN=[0;0;1]; e=uN-(u1'*uN)*u1; uN1=e/norm(e); uE=cross(uN1,u1); G=[uE,uN1,u1];
end

disp('--- (3) sign of the clock term: estimated clock vs delta_t ---')
rng(1001);
o = GPS_spoofing_sim_20250617(struct('PLOT',false,'grid_step',10,'sigma_pr',0,'Pt_spoof',0));
% with no attack and no noise the solver must return the true positions;
% recompute the fix for one drone and compare the clock state with delta_t
printf('H0 noiseless: max |lat err| %.2e deg, delta_t(1) = %+0.4f s\n', ...
  max(max(abs(o.coords_est(:,1:3:end)-o.coords_grid(:,2:3:end)))), o.delta_t(1,1));

disp('--- (1)+(2) D and M: recompute from scratch and compare ---')
rng(1001);
o = GPS_spoofing_sim_20250617(struct('PLOT',false,'grid_step',10,'sigma_pr',0,'max_spoofed_sats',inf,'R_spoofers',10000));
Nt=size(o.coords_est,1); Nu=25; uc=o.u_center; sc=reshape(o.spoofer_coords,3,[])';
[I,J]=find(triu(true(Nu),1));
med_ratio=zeros(Nt,1); med_diff=zeros(Nt,1); Epred=zeros(Nt,1); thk=zeros(Nt,1); elmin=zeros(Nt,1);
for t=1:Nt
  tt=o.coords_grid(t,1); pc=frame(tt,o.coords(t,2),o.coords(t,3),o.coords(t,4)); G=enu(pc);
  vis=find(o.elevation(t,:,uc)>0); el=o.elevation(t,vis,uc); az=o.azimuth(t,vis,uc);
  E=[cos(el).*cos(az); cos(el).*sin(az); sin(el)]; V=zeros(3,numel(vis)); th=zeros(1,numel(vis));
  for n=1:numel(vis)
    k=o.spoof_for_sat(t,vis(n)); q=frame(tt,sc(k,1),sc(k,2),sc(k,3)); v=G'*(q-pc); V(:,n)=v/norm(v); th(n)=asin(V(3,n));
  end
  Gm=[-E', ones(numel(vis),1)]; B=(E-V)';           % eq. (10): G rows [-e', 1], B rows (e-v)'
  X=(Gm'*Gm)\(Gm'*B); M=X(1:3,1:3);                  % M = [(G'G)^-1 G' B]_{1:3}
  P=zeros(Nu,3); Pe=zeros(Nu,3);
  for u=1:Nu
    P(u,:)=(G'*(frame(tt,o.coords_grid(t,2+3*(u-1)),o.coords_grid(t,3+3*(u-1)),o.coords_grid(t,4+3*(u-1)))-pc))';
    Pe(u,:)=(G'*(frame(tt,o.coords_est(t,1+3*(u-1)),o.coords_est(t,2+3*(u-1)),o.coords_est(t,3+3*(u-1)))-pc))';
  end
  D3=P(J,:)-P(I,:); L=sqrt(sum(D3.^2,2));
  Dn=D3+(M*D3')'; Epred(t)=median(abs(sqrt(sum(Dn.^2,2))-L)./L);
  de=sqrt(sum((Pe(I,:)-Pe(J,:)).^2,2)); dt=sqrt(sum((P(I,:)-P(J,:)).^2,2));
  med_ratio(t)=median(de./dt)-1; med_diff(t)=median((de-dt)./dt);
  thk(t)=median(th)*180/pi; elmin(t)=min(el)*180/pi;
end
printf('D (median over epochs of median over pairs of (dhat-d)/d) = %+0.4f\n', median(med_diff));
printf('same from the ratio form (median r - 1)                   = %+0.4f  [identical: %d]\n', median(med_ratio), isequal(round(1e9*med_ratio),round(1e9*med_diff)));
printf('M as in eq. (10) with full 3D pair vectors: E_pred %.4f vs measured E %.4f\n', median(Epred), median(abs(med_diff)));

disp('--- (4) satellite re-routing between transmitters ---')
sw=0; tot=0;
for t=2:Nt
  v=find(o.elevation(t,:,uc)>0 & o.elevation(t-1,:,uc)>0);
  sw=sw+sum(o.spoof_for_sat(t,v)~=o.spoof_for_sat(t-1,v)); tot=tot+numel(v);
end
printf('routing switches: %d over %d satellite-epochs (%.2f%%), i.e. %.2f per epoch\n', sw, tot, 100*sw/tot, sw/(Nt-1));
d=abs(diff(med_diff)); [mx,it]=max(d);
printf('largest epoch-to-epoch jump of the per-epoch dilation: %.4f (epoch %d); median |jump| %.4f\n', mx, it+1, median(d));

disp('--- (5) depression angle: does the residual vanish for theta_k = -el_s? ---')
% place the transmitters 1 km BELOW the swarm at short horizontal range, so
% that |theta_k| approaches the satellite elevations from below
for Rs=[1000 2000 5000]
  rng(1001);
  ob = GPS_spoofing_sim_20250617(struct('PLOT',false,'grid_step',10,'sigma_pr',0,'max_spoofed_sats',inf, ...
       'R_spoofers',Rs,'spoofer_alt',-900,'channel','freespace','Pt_spoof',1e-2,'traj_scale',0.1));
  Nt2=size(ob.coords_est,1); scb=reshape(ob.spoofer_coords,3,[])'; md=zeros(Nt2,1); tk=zeros(Nt2,1);
  for t=1:Nt2
    tt=ob.coords_grid(t,1); pc=frame(tt,ob.coords(t,2),ob.coords(t,3),ob.coords(t,4)); G=enu(pc);
    vis=find(ob.elevation(t,:,uc)>0); th=zeros(1,numel(vis));
    for n=1:numel(vis)
      k=ob.spoof_for_sat(t,vis(n)); q=frame(tt,scb(k,1),scb(k,2),scb(k,3)); v=G'*(q-pc); th(n)=asin(v(3)/norm(v));
    end
    P=zeros(Nu,3); Pe=zeros(Nu,3);
    for u=1:Nu
      P(u,:)=(G'*(frame(tt,ob.coords_grid(t,2+3*(u-1)),ob.coords_grid(t,3+3*(u-1)),ob.coords_grid(t,4+3*(u-1)))-pc))';
      Pe(u,:)=(G'*(frame(tt,ob.coords_est(t,1+3*(u-1)),ob.coords_est(t,2+3*(u-1)),ob.coords_est(t,3+3*(u-1)))-pc))';
    end
    de=sqrt(sum((Pe(I,:)-Pe(J,:)).^2,2)); dt=sqrt(sum((P(I,:)-P(J,:)).^2,2));
    md(t)=median((de-dt)./dt); tk(t)=median(th)*180/pi;
  end
  printf('transmitters 1 km below, R=%5d m | theta_k med %+6.1f deg | signed dilation %+0.4f\n', Rs, median(tk), median(md));
end

disp('--- (6) 2.5 km array: inside vs outside ---')
rng(1001);
oc = GPS_spoofing_sim_20250617(struct('PLOT',false,'grid_step',10,'sigma_pr',0,'max_spoofed_sats',inf,'R_spoofers',2500));
Nt3=size(oc.coords_est,1); scc=reshape(oc.spoofer_coords,3,[])'; R=6371000; cc=oc.center_coords;
din=zeros(Nt3,1); md=zeros(Nt3,1); tk=zeros(Nt3,1); azs=zeros(Nt3,1);
for t=1:Nt3
  tt=oc.coords_grid(t,1); pc=frame(tt,oc.coords(t,2),oc.coords(t,3),oc.coords(t,4)); G=enu(pc);
  c=oc.coords(t,2:3); din(t)=R*norm([(c(1)-cc(1))*pi/180,(c(2)-cc(2))*pi/180*cos(cc(1)*pi/180)]);
  vis=find(oc.elevation(t,:,uc)>0); th=zeros(1,numel(vis)); daz=zeros(1,numel(vis));
  for n=1:numel(vis)
    k=oc.spoof_for_sat(t,vis(n)); q=frame(tt,scc(k,1),scc(k,2),scc(k,3)); v=G'*(q-pc); v=v/norm(v);
    th(n)=asin(v(3)); daz(n)=abs(mod(atan2(v(2),v(1))-oc.azimuth(t,vis(n),uc)+pi,2*pi)-pi);
  end
  P=zeros(Nu,3); Pe=zeros(Nu,3);
  for u=1:Nu
    P(u,:)=(G'*(frame(tt,oc.coords_grid(t,2+3*(u-1)),oc.coords_grid(t,3+3*(u-1)),oc.coords_grid(t,4+3*(u-1)))-pc))';
    Pe(u,:)=(G'*(frame(tt,oc.coords_est(t,1+3*(u-1)),oc.coords_est(t,2+3*(u-1)),oc.coords_est(t,3+3*(u-1)))-pc))';
  end
  de=sqrt(sum((Pe(I,:)-Pe(J,:)).^2,2)); dt=sqrt(sum((P(I,:)-P(J,:)).^2,2));
  md(t)=median((de-dt)./dt); tk(t)=median(th)*180/pi; azs(t)=median(daz)*180/pi;
end
in=din<2500;
printf('inside  (%3d epochs): dilation %+0.4f | theta_k %+6.1f deg | azimuth mismatch %5.1f deg\n', sum(in), median(md(in)), median(tk(in)), median(azs(in)));
printf('outside (%3d epochs): dilation %+0.4f | theta_k %+6.1f deg | azimuth mismatch %5.1f deg\n', sum(~in), median(md(~in)), median(tk(~in)), median(azs(~in)));
disp('VERIFY DONE');

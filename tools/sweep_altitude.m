1;
% Noiseless probe: dilation (measured + first-order model) vs transmitter
% deployment: circle at R=10 km with altitude sweep, and LOS placement
% (WASP-like) at several slant ranges. g = 10 m, seed 1001.
% run from the repository root (needs physconst, e.g. Phased Array / Antenna Toolbox)
function p = frame(t,lat,lon,alt)
  R_E = 6371000; T_E = 86164; th0 = 0.23*2*pi; la = lat*pi/180; lo = lon*pi/180;
  p = (R_E+alt)*[cos(la)*cos(2*pi*t/T_E+th0+lo); cos(la)*sin(2*pi*t/T_E+th0+lo); sin(la)];
end
function G = enu_frame(p)
  u1 = p/norm(p); uN = [0;0;1]; eN = uN - (u1'*uN)*u1; uN1 = eN/norm(eN); uE = cross(uN1,u1); G = [uE, uN1, u1];
end
function [E_meas, dil, E_pred, thk, elmean, offc, satfrac] = evaluate(o)
  Nt = size(o.coords_est,1); Nu = 25; uc = o.u_center; sc = reshape(o.spoofer_coords,3,[])';
  [I,J] = find(triu(true(Nu),1)); E_meas = zeros(Nt,1); dil = E_meas; E_pred = E_meas; thk = E_meas; elmean = E_meas; offc = E_meas;
  for t = 1:Nt
    tt = o.coords_grid(t,1); pc = frame(tt,o.coords(t,2),o.coords(t,3),o.coords(t,4)); G = enu_frame(pc);
    vis = find(o.elevation(t,:,uc) > 0); el = o.elevation(t,vis,uc); az = o.azimuth(t,vis,uc);
    E = [cos(el).*cos(az); cos(el).*sin(az); sin(el)]; V = zeros(3,numel(vis)); th = zeros(1,numel(vis));
    for n = 1:numel(vis)
      k = o.spoof_for_sat(t,vis(n)); q = frame(tt,sc(k,1),sc(k,2),sc(k,3)); v = G'*(q-pc); V(:,n) = v/norm(v); th(n) = asin(V(3,n));
    end
    H = [-E', ones(numel(vis),1)]; M = zeros(3,2);
    for c = 1:2, d = zeros(3,1); d(c) = 1; x = H\((E - V)'*d); M(:,c) = x(1:3); end
    P = zeros(Nu,3); Pe = zeros(Nu,3);
    for u = 1:Nu
      P(u,:) = (G'*(frame(tt,o.coords_grid(t,2+3*(u-1)),o.coords_grid(t,3+3*(u-1)),o.coords_grid(t,4+3*(u-1)))-pc))';
      Pe(u,:) = (G'*(frame(tt,o.coords_est(t,1+3*(u-1)),o.coords_est(t,2+3*(u-1)),o.coords_est(t,3+3*(u-1)))-pc))';
    end
    D = P(J,1:2)-P(I,1:2); L = sqrt(sum(D.^2,2)); S = D./L;
    Dn = [S, zeros(size(S,1),1)] + (M*S')'; stretch = sqrt(sum(Dn.^2,2)) - 1; E_pred(t) = median(abs(stretch));
    de = sqrt(sum((Pe(I,:)-Pe(J,:)).^2,2)); dt = sqrt(sum((P(I,:)-P(J,:)).^2,2));
    E_meas(t) = median(abs(de-dt)./dt); dil(t) = median(de./dt)-1; thk(t) = median(th)*180/pi; elmean(t) = mean(el)*180/pi;
    offc(t) = norm(mean(Pe)-mean(P));
  end
  satfrac = o.theta_sat_frac;
end
fid = fopen('probe_alt2.csv','w');
fprintf(fid,'label,param,theta_k_med,el_mean,E_meas,dil_signed,E_pred,offset_m,margin_min_dB\n');
function m = minmargin(o), d = o.Pr_sp_dBW - o.Pr_sat_dBW; d = d(isfinite(d)); m = min(d); end
% (1) altitude sweep, full mission, free-space channel, 10 mW
for h = [0.1 0.5 1 2 3.5 5 7]*1000
  rng(1001); o = GPS_spoofing_sim_20250617(struct('PLOT',false,'grid_step',10,'sigma_pr',0,'max_spoofed_sats',inf,'R_spoofers',10000,'spoofer_alt',h,'channel','freespace','Pt_spoof',1e-2));
  [E,dl,Ep,thk,elm,off,sf] = evaluate(o);
  fprintf(fid,'circle,%g,%.3f,%.3f,%.4f,%.4f,%.4f,%.0f,%.1f\n',h,median(thk),median(elm),median(E),median(dl),median(Ep),median(off),minmargin(o));
  printf('circle FS h=%5.0f m | theta_k med %5.1f | E_meas %.3f  dil %+.3f  E_pred %.3f | offset %.0f m | min margin %.1f dB\n',h,median(thk),median(E),median(dl),median(Ep),median(off),minmargin(o));
end
% (2) loitering mission (traj_scale 0.1): ground circle (Al-Hourani) and LOS placement (free space)
rng(1001); o = GPS_spoofing_sim_20250617(struct('PLOT',false,'grid_step',10,'sigma_pr',0,'max_spoofed_sats',inf,'R_spoofers',10000,'traj_scale',0.1));
[E,dl,Ep,thk,elm,off,sf] = evaluate(o);
fprintf(fid,'loiter_circle,10000,%.3f,%.3f,%.4f,%.4f,%.4f,%.0f,%.1f\n',median(thk),median(elm),median(E),median(dl),median(Ep),median(off),minmargin(o));
printf('loiter circle ground | theta_k %5.1f | E_meas %.3f dil %+.3f E_pred %.3f | offset %.0f m | margin %.1f dB\n',median(thk),median(E),median(dl),median(Ep),median(off),minmargin(o));
for Rl = [2 5 10]*1000
  rng(1001); o = GPS_spoofing_sim_20250617(struct('PLOT',false,'grid_step',10,'sigma_pr',0,'max_spoofed_sats',inf,'placement','los','R_los',Rl,'channel','freespace','Pt_spoof',1e-2,'traj_scale',0.1));
  [E,dl,Ep,thk,elm,off,sf] = evaluate(o);
  fprintf(fid,'loiter_los,%g,%.3f,%.3f,%.4f,%.4f,%.4f,%.0f,%.1f\n',Rl,median(thk),median(elm),median(E),median(dl),median(Ep),median(off),minmargin(o));
  printf('loiter LOS R=%5.0f m | theta_k %5.1f | E_meas %.4f dil %+.4f E_pred %.4f | offset %.0f m | margin %.1f dB | N_sp %d\n',Rl,median(thk),median(E),median(dl),median(Ep),median(off),minmargin(o),numel(o.spoofer_coords)/3);
end
fclose(fid); disp('PROBE2 DONE');

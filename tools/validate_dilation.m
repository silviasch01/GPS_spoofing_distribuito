1;
% A2: per-epoch predicted dilation (first-order LS model) vs measured
% dilation, on a NOISELESS attacked mission at g = 20 m; plus B5: noiseless
% scenario dump at g = 10 m (true grid, estimated grid, offset).
% run from the repository root (needs physconst, e.g. Phased Array / Antenna Toolbox)
R_E = 6371000; T_E = 86164; th0 = 0.23*2*pi;

function p = frame(t,lat,lon,alt)
  R_E = 6371000; T_E = 86164; th0 = 0.23*2*pi;
  la = lat*pi/180; lo = lon*pi/180;
  p = (R_E+alt)*[cos(la)*cos(2*pi*t/T_E+th0+lo); cos(la)*sin(2*pi*t/T_E+th0+lo); sin(la)];
end
function G = enu_frame(p)   % columns: East, North, Up at point p (same as gps_trace)
  u1 = p/norm(p); uN = [0;0;1]; eN = uN - (u1'*uN)*u1; uN1 = eN/norm(eN); uE = cross(uN1,u1);
  G = [uE, uN1, u1];
end

% ---------------- (a) validation at g = 20 m, sigma_pr = 0 ----------------
rng(1001);
o = GPS_spoofing_sim_20250617(struct('PLOT',false,'grid_step',20,'sigma_pr',0,'max_spoofed_sats',inf,'R_spoofers',10000));
Nt = size(o.coords_est,1); Nu = 25; uc = o.u_center;
sc = reshape(o.spoofer_coords,3,[])';
[I,J] = find(triu(true(Nu),1));
out = zeros(Nt,6);
for t = 1:Nt
  tt = o.coords_grid(t,1);
  pc = frame(tt, o.coords(t,2), o.coords(t,3), o.coords(t,4));
  G = enu_frame(pc);
  vis = find(o.elevation(t,:,uc) > 0);
  el = o.elevation(t,vis,uc); az = o.azimuth(t,vis,uc);
  % satellite LOS (ENU) and assigned-spoofer direction (ENU) from the centre
  E = [cos(el).*cos(az); cos(el).*sin(az); sin(el)];          % 3 x Ns
  V = zeros(3,numel(vis));
  for n = 1:numel(vis)
    k = o.spoof_for_sat(t,vis(n));
    q = frame(tt, sc(k,1), sc(k,2), sc(k,3));
    v = G'*(q-pc); V(:,n) = v/norm(v);
  end
  % pseudorange error for a unit horizontal displacement d: eps_s = (e_s - v_k).d
  % LS: eps = H*[dp; dt'] with H = [-e_s', 1]
  H = [-E', ones(numel(vis),1)];
  M = zeros(3,2);
  for c = 1:2
    d = zeros(3,1); d(c) = 1;
    eps = (E - V)'*d;
    x = H\eps; M(:,c) = x(1:3);
  end
  % pair directions of the true grid (ENU, horizontal)
  P = zeros(Nu,3);
  for u = 1:Nu
    P(u,:) = (G'*(frame(tt,o.coords_grid(t,2+3*(u-1)),o.coords_grid(t,3+3*(u-1)),o.coords_grid(t,4+3*(u-1))) - pc))';
  end
  D = P(J,1:2)-P(I,1:2); L = sqrt(sum(D.^2,2)); S = D./L;
  % exact relative length change of each pair under the linear map
  % delta -> delta + M*delta (3-D, including the induced altitude tilt)
  Dn = [S, zeros(size(S,1),1)] + (M*S')';     % mapped unit pair vectors (3-D)
  stretch = sqrt(sum(Dn.^2,2)) - 1;
  E_pred = median(abs(stretch));
  Mh = M(1:2,1:2); stretch_h = sum((S*Mh').*S,2);   % horizontal first-order only
  E_pred_h = median(abs(stretch_h));
  % measured
  Pe = zeros(Nu,3); Pt = zeros(Nu,3);
  for u = 1:Nu
    Pe(u,:) = frame(tt,o.coords_est(t,1+3*(u-1)),o.coords_est(t,2+3*(u-1)),o.coords_est(t,3+3*(u-1)))';
    Pt(u,:) = frame(tt,o.coords_grid(t,2+3*(u-1)),o.coords_grid(t,3+3*(u-1)),o.coords_grid(t,4+3*(u-1)))';
  end
  de = sqrt(sum((Pe(I,:)-Pe(J,:)).^2,2)); dt = sqrt(sum((Pt(I,:)-Pt(J,:)).^2,2));
  E_meas = median(abs(de-dt)./dt); dil_meas = median(de./dt)-1;
  E_simple = 1/mean(cos(el)) - 1;
  out(t,:) = [tt-o.coords_grid(1,1), E_meas, E_pred, E_simple, dil_meas, E_pred_h];
end
dlmwrite('valid_g20.csv', out);
printf('validation g=20 noiseless: median E_meas %.3f, E_pred(3D) %.3f, E_pred(horiz) %.3f, E_simple %.3f, dil_meas %.3f | corr(E_pred,E_meas) %.3f | max |E_pred-E_meas| %.3f | median ratio pred/meas %.3f\n', ...
  median(out(:,2)), median(out(:,3)), median(out(:,6)), median(out(:,4)), median(out(:,5)), corr(out(:,3),out(:,2)), max(abs(out(:,3)-out(:,2))), median(out(:,3)./out(:,2)));
printf('theta_max_deg (g=20, R=10km): %.2f ; saturated fraction %.3f\n', o.theta_max_deg, o.theta_sat_frac);

% ---------------- (b) noiseless scenario at g = 10 m ----------------
rng(1001);
o = GPS_spoofing_sim_20250617(struct('PLOT',false,'grid_step',10,'sigma_pr',0,'max_spoofed_sats',inf,'R_spoofers',10000));
Nt = size(o.coords_est,1); tm = round(Nt/2);
g = o.coords_grid(tm,2:end); ge = o.coords_est(tm,:); gl = zeros(Nu,2); gle = zeros(Nu,2);
for u = 1:Nu, gl(u,:) = [g(3*u-1), g(3*u-2)]; gle(u,:) = [ge(3*u-1), ge(3*u-2)]; end
dlmwrite('scn4_grid_true.csv', gl); dlmwrite('scn4_grid_est.csv', gle);
dlmwrite('scn4_traj.csv',[o.coords(:,3),o.coords(:,2)]); dlmwrite('scn4_spoof.csv',[o.coords_spoof(:,3),o.coords_spoof(:,2)]);
dlmwrite('scn4_spoofers.csv',[sc(:,2),sc(:,1)]);
% dilation of the estimated grid vs true grid at that epoch
Pe = zeros(Nu,3); Pt = zeros(Nu,3); tt = o.coords_grid(tm,1);
for u = 1:Nu
  Pe(u,:) = frame(tt,ge(3*u-2),ge(3*u-1),ge(3*u))'; Pt(u,:) = frame(tt,g(3*u-2),g(3*u-1),g(3*u))';
end
de = sqrt(sum((Pe(I,:)-Pe(J,:)).^2,2)); dt = sqrt(sum((Pt(I,:)-Pt(J,:)).^2,2));
printf('scenario g=10 noiseless, epoch %d: median ratio %.3f, median |dhat-d|/d %.3f, offset %.0f m, theta_max %.2f deg\n', tm, median(de./dt), median(abs(de-dt)./dt), norm(mean(Pe)-mean(Pt)), o.theta_max_deg);

1;
% run from the repository root (needs physconst, e.g. Phased Array / Antenna Toolbox)
% reproduce mission 1 at g=50 m (attack) with the CURRENT code and compare
% with the saved cell (seed 1000*1+1 = 1001, same options as the driver)
L = load('cell_g50.mat');
rng(1001);
o = GPS_spoofing_sim_20250617(struct('PLOT',false,'grid_step',50,'grid_points',5,'sigma_pr',3,'R_spoofers',10000,'N_spoofers',24,'max_spoofed_sats',inf));
% recompute the statistics exactly as distance_tests (ranging noise drawn after the simulation)
Nt=size(o.coords_est,1); Nu=25; [I,J]=find(triu(true(Nu),1)); R=6371000;
P=@(c)(R+c(3:3:end)').*[cos(c(1:3:end)'*pi/180).*cos(c(2:3:end)'*pi/180), cos(c(1:3:end)'*pi/180).*sin(c(2:3:end)'*pi/180), sin(c(1:3:end)'*pi/180)];
T=zeros(Nt,4);
for t=1:Nt
  pt=P(o.coords_grid(t,2:end)); pe=P(o.coords_est(t,:));
  d0=sqrt(sum((pt(I,:)-pt(J,:)).^2,2)); d0=d0+0.1*randn(size(d0)); dg=sqrt(sum((pe(I,:)-pe(J,:)).^2,2)); r=dg./d0;
  T(t,:)=[min(r) max(r) median(r) max(abs(r-1))];
end
S = L.T{1,2}(L.RID{1,2}==1,1:4);
printf('saved cell: det_names=%d, Nruns=%d, rows mission1=%d\n', numel(L.det_names), L.Nruns, size(S,1));
printf('max |T_now - T_saved| over 4 statistics x %d epochs: %.3e   (median ratio now %.4f, saved %.4f)\n', Nt, max(abs(T(:)-S(:))), median(T(:,3)), median(S(:,3)));

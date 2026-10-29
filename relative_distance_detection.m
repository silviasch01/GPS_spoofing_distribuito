% relative_distance_detection.m
%
% Evaluation of classical MULTI-RECEIVER RELATIVE-DISTANCE spoofing detectors
% against the PHYSICAL distributed ground spoofer of GPS_spoofing_sim_20250617.m
%
% ------------------------------------------------------------------ detectors
% A swarm of m receivers flies in a known formation. Each receiver u computes
% its GNSS fix p_hat(u); the receivers also measure their true mutual
% distances with an independent ranging sensor (UWB), d_meas_ij. All the
% detectors work on the normalized ratios  r_ij = d_hat_ij / d_meas_ij :
%
%   Jansen et al. (ACSAC 2016) -- COLLAPSE test. A single-transmitter spoofer
%     forces all receivers onto the same position, so the GNSS distances
%     shrink (r_ij -> 0). Alarm when f(r) < lambda, with f one of:
%        minimal   : f = min_ij r_ij
%        maximal   : f = max_ij r_ij
%        majority  : f = median_ij r_ij   (>=half of the pairs below lambda,
%                    i.e. the paper's normalized majority voting)
%   Formation CONSISTENCY test (two-sided), same family. Alarm when the
%     distances deviate (either way) from the measured formation:
%        T = max_ij |r_ij - 1| > lambda.
%     "Worst pair" form of the GNSS-vs-ranging check used for UAV swarms by
%     Mykytyn et al. (MECO 2023, |d_GPS - d_UWB| > thr) and, in the
%     measurement domain, by Stenberg et al. (NAVIGATION 2022) and Swaszek,
%     Hartnett & Seals (ION GNSS+ 2017). It is NOT the Borio & Gioia
%     sum-of-squares detector (TAES 2016), a two-receiver carrier-phase test.
%
% SemperFi (Sathaye et al., NDSS 2022), GNSS/INS fusion, is NOT evaluated
% (see the note at the bottom of this file).
%
% --------------------------------------------------------------- attack model
%   H0    : no spoofing (used only to calibrate the thresholds)
%   power : the physical distributed attack of the simulator: N_spoofers
%           fixed ground transmitters on a circle of radius R_spoofers around
%           the operating area, one impersonated satellite per transmitter
%           (routing re-evaluated every epoch from the current swarm centre),
%           capture decided by received power through the Al-Hourani
%           ground-to-UAV channel, replayed delay calibrated at the true swarm
%           centre plus the true 3D slant-path difference seen by each drone.
%           No per-victim pre-compensation and no idealized/oracle attack is
%           simulated: the residual formation distortion the detectors see is
%           a property of the physical geometry.
%
% ------------------------------------------------------------------ statistics
% Thresholds are calibrated on H0 (same GNSS/ranging noise) for a target
% false-alarm probability Pfa (Neyman-Pearson). Pd is the fraction of attacked
% epochs that raise an alarm. Because the epochs of one mission share clocks,
% trajectory and satellites, the Monte Carlo unit is the MISSION: Nruns
% independent missions are simulated, and 95% confidence intervals on Pd are
% obtained by a paired bootstrap over missions (H0 and attack of the same
% mission share the random seed; the threshold is re-calibrated in every
% bootstrap resample).

clear
close all

% ---------------------------------------------------------------- parameters
% formation spacing [m]: the 5-20 m range covers the typical inter-UAV
% separation of outdoor GPS-based formation flight (5 m in Cheok et al.,
% SAE Int. J. Aerospace 2022, and Grundler & Ryan, Cal Poly 2025; Guo et
% al., IEEE TAES 2025, derive the minimum spacing from navigation+control
% error, metres with GPS); sub-metre spacings require RTK (light shows)
grid_steps = [5 8 10 12 15 20];
grid_points = 5;                 % 5x5 receivers (m = 25)
Nruns = 20;                      % Monte Carlo missions (random start time, clocks, noise)
sigma_pr = 3;                    % pseudorange noise std [m]
sigma_range = 0.1;               % inter-receiver ranging noise std [m]
Pfa = 0.01;                      % target false alarm probability
roc_step = 10;                   % grid_step used for the ROC figure
Nboot = 500;                     % bootstrap resamples (over missions)
% attack deployment: the ground array must enclose the operating area of the
% swarm (the mission trajectory spans ~7 km from its start point)
R_spoofers = 10000;              % [m]
N_spoofers = 24;

scen(1) = struct('name','H0',    'opts',struct('Pt_spoof',0));
scen(2) = struct('name','power', 'opts',struct('max_spoofed_sats',inf));
det_names = {'Jansen min','Jansen max','Jansen majority','Formation consistency'};
% collapse detectors alarm for SMALL statistics (f < lambda), the consistency
% detector for LARGE statistics (T > lambda)
alarm_small = [true true true false];

Nscen = numel(scen);
Ndet = numel(det_names);
Ngs = numel(grid_steps);

% ------------------------------------------------------------- simulations
T = cell(Ngs,Nscen);    % test statistics, one row per epoch
E = cell(Ngs,Nscen);    % median relative distance error per epoch
RID = cell(Ngs,Nscen);  % mission id of every epoch (for the bootstrap)
theta_sat = zeros(Ngs,Nruns);
for g = 1:Ngs
	for s = 1:Nscen
		T{g,s} = []; E{g,s} = []; RID{g,s} = [];
		for r = 1:Nruns
			% same random realization (start time, clocks, noise) for H0 and
			% attack of a given mission
			seed = 1000*r + g;
			rng(seed);
			opts = scen(s).opts;
			opts.PLOT = false;
			opts.grid_step = grid_steps(g);
			opts.grid_points = grid_points;
			opts.sigma_pr = sigma_pr;
			opts.R_spoofers = R_spoofers;
			opts.N_spoofers = N_spoofers;
			out = GPS_spoofing_sim_20250617(opts);
			if isfield(out,'theta_sat_frac') && s == 2
				theta_sat(g,r) = out.theta_sat_frac;
			end

			[Tr,Er] = distance_tests(out,sigma_range);
			T{g,s} = [T{g,s}; Tr];
			E{g,s} = [E{g,s}; Er];
			RID{g,s} = [RID{g,s}; r*ones(size(Er))];
			fprintf('grid_step %4g m | %-6s | mission %2d/%d | median rel. dist. error %6.2f %%\n',...
				grid_steps(g),scen(s).name,r,Nruns,100*median(Er));
		end
	end
end

% -------------------------------------------- Pd at given Pfa, bootstrap CI
Pd = zeros(Ngs,Ndet);          % point estimate (all missions)
Pd_ci = zeros(Ngs,Ndet,2);     % 95% CI from the paired mission bootstrap
for g = 1:Ngs
	for k = 1:Ndet
		Pd(g,k) = pd_at_pfa(T{g,1}(:,k),T{g,2}(:,k),Pfa,alarm_small(k));
		pb = zeros(Nboot,1);
		for b = 1:Nboot
			rs = randi(Nruns,Nruns,1);
			i0 = ismember_runs(RID{g,1},rs);
			i1 = ismember_runs(RID{g,2},rs);
			pb(b) = pd_at_pfa(T{g,1}(i0,k),T{g,2}(i1,k),Pfa,alarm_small(k));
		end
		Pd_ci(g,k,:) = quantile_(pb,[0.025 0.975]);
	end
end

fprintf('\nDetection probability at Pfa = %g (sigma_pr = %g m, sigma_range = %g m)\n',...
	Pfa,sigma_pr,sigma_range);
fprintf('%d missions per spacing, 95%% CI from paired bootstrap over missions\n',Nruns);
for k = 1:Ndet
	fprintf('\n%s\n  grid_step [m]      Pd      [95%% CI]\n',det_names{k});
	for g = 1:Ngs
		fprintf('  %13g   %6.3f   [%6.3f, %6.3f]\n',grid_steps(g),Pd(g,k),Pd_ci(g,k,1),Pd_ci(g,k,2));
	end
end

fprintf('\nMedian relative distance error |d_hat-d|/d [%%] (median over all epochs)\n  grid_step [m] ');
fprintf('%10s',scen.name);
fprintf('\n');
for g = 1:Ngs
	fprintf('  %13g ',grid_steps(g));
	for s = 1:Nscen
		fprintf('%10.2f',100*median(E{g,s}));
	end
	fprintf('\n');
end
fprintf('\nFraction of spoofer links with elevation angle > 10 deg (Al-Hourani fit saturated): %.3f (mean over missions)\n', mean(theta_sat(:)));

% ------------------------------------------------------------------ figures
colors = lines(4);
figure(1)
for k = 1:Ndet
	subplot(2,2,k)
	lo = squeeze(Pd_ci(:,k,1)); hi = squeeze(Pd_ci(:,k,2));
	fill([grid_steps fliplr(grid_steps)],[lo' fliplr(hi')],colors(2,:),'FaceAlpha',0.2,'EdgeColor','none');
	hold on
	semilogx(grid_steps,Pd(:,k),'-o','Color',colors(2,:),'LineWidth',1.5);
	semilogx(grid_steps,Pfa*ones(1,Ngs),'k--');
	set(gca,'XScale','log')
	grid on; ylim([0 1.05])
	xlabel('formation spacing [m]'); ylabel('P_D')
	title(sprintf('%s (P_{FA}=%g)',det_names{k},Pfa))
	if k==1, legend({'95% CI','physical attack','P_{FA}'},'Location','best'); end
end

g = find(grid_steps==roc_step,1);
if ~isempty(g)
	figure(2)
	for k = 1:Ndet
		subplot(2,2,k)
		t0 = sort(T{g,1}(~isnan(T{g,1}(:,k)),k));
		thr_c = [-inf;t0;inf];
		tk = T{g,2}(:,k); tk = tk(~isnan(tk));
		pfa_c = zeros(numel(thr_c),1); pd_c = pfa_c;
		for i = 1:numel(thr_c)
			if alarm_small(k)
				pfa_c(i) = mean(t0 < thr_c(i)); pd_c(i) = mean(tk < thr_c(i));
			else
				pfa_c(i) = mean(t0 > thr_c(i)); pd_c(i) = mean(tk > thr_c(i));
			end
		end
		plot(pfa_c,pd_c,'Color',colors(2,:),'LineWidth',1.5); hold on
		plot([0 1],[0 1],'k--'); grid on; axis square
		xlabel('P_{FA}'); ylabel('P_D')
		title(sprintf('ROC %s, spacing %g m',det_names{k},roc_step))
		if k==1, legend({'physical attack','chance'},'Location','southeast'); end
	end
end

figure(3)
semilogx(grid_steps,cellfun(@(e) 100*median(e),E(:,1)),'--s','Color',[.5 .5 .5],'LineWidth',1.5); hold on
semilogx(grid_steps,cellfun(@(e) 100*median(e),E(:,2)),'-o','Color',colors(2,:),'LineWidth',1.5);
grid on
xlabel('formation spacing [m]'); ylabel('median |d_{hat}-d|/d [%]')
legend({'H0 (noise floor)','physical attack'},'Location','best')
title('Distortion of the relative distances')


% ======================================================== local functions

% Jansen collapse variants + formation consistency test, per epoch.
% Columns of Tr: [min(r), max(r), median(r), max|r-1|].
function [Tr,Er] = distance_tests(out,sigma_range)
	Nt = size(out.coords_est,1);
	Nusers = size(out.coords_est,2)/3;
	[I,J] = find(triu(true(Nusers),1));
	Tr = zeros(Nt,4);
	Er = zeros(Nt,1);
	for t = 1:Nt
		p_true = coords_to_ecef(out.coords_grid(t,2:end));
		p_est = coords_to_ecef(out.coords_est(t,:));
		d0_clean = sqrt(sum((p_true(I,:)-p_true(J,:)).^2,2));  % true distances
		d0 = d0_clean + sigma_range*randn(size(d0_clean));     % UWB measurement
		d_gnss = sqrt(sum((p_est(I,:)-p_est(J,:)).^2,2));      % GNSS distances
		r = d_gnss ./ d0;                                      % normalized ratio
		Tr(t,1) = min(r);
		Tr(t,2) = max(r);
		Tr(t,3) = median(r);
		Tr(t,4) = max(abs(r-1));
		Er(t) = median(abs(d_gnss-d0_clean)./d0_clean);
	end
end

% Pd at the threshold giving Pfa on the H0 statistics t0 (empirical quantile)
function pd = pd_at_pfa(t0,t1,Pfa,alarm_small)
	t0 = t0(~isnan(t0)); t1 = t1(~isnan(t1));
	if alarm_small
		thr = lower_quantile(t0,Pfa); pd = mean(t1 < thr);
	else
		thr = upper_quantile(t0,Pfa); pd = mean(t1 > thr);
	end
end

% logical index of the epochs whose mission id is in the resample rs
% (a mission drawn twice contributes twice)
function idx = ismember_runs(rid,rs)
	idx = [];
	for r = rs(:)'
		idx = [idx; find(rid == r)]; %#ok<AGROW>
	end
end

% [lat1,lon1,alt1,lat2,...] -> (#users,3) Cartesian positions (spherical Earth,
% same model as the simulator; the common Earth-rotation term does not change
% the mutual distances).
function p = coords_to_ecef(c)
	R_Earth = physconst('EarthRadius');
	lat = c(1:3:end).'*pi/180;
	lon = c(2:3:end).'*pi/180;
	alt = c(3:3:end).';
	p = (R_Earth+alt).*[cos(lat).*cos(lon), cos(lat).*sin(lon), sin(lat)];
end

% threshold with a fraction Pfa of x ABOVE it (consistency test)
function thr = upper_quantile(x,Pfa)
	x = sort(x(:)); n = numel(x);
	thr = x(min(n,max(1,ceil((1-Pfa)*n))));
end

% threshold with a fraction Pfa of x BELOW it (collapse test)
function thr = lower_quantile(x,Pfa)
	x = sort(x(:)); n = numel(x);
	thr = x(min(n,max(1,floor(Pfa*n)+1)));
end

% empirical quantiles without the statistics toolbox
function q = quantile_(x,p)
	x = sort(x(:)); n = numel(x);
	q = x(min(n,max(1,round(p*n))));
end

% ---------------------------------------------------------------------------
% NOTE on SemperFi (GNSS/INS fusion, Sathaye et al., NDSS 2022) -- not tested.
%
% SemperFi commands a pseudorandom maneuver and compares the trajectory tracked
% by the IMU (dead-reckoning) against the one reported by GPS, raising an alarm
% when the averaged Euclidean deviation exceeds ~5 m. Its statistic is built on
% displacements from the manoeuvre start, so it responds to the time-varying
% part of the induced position error. Whether the physical distributed attack
% simulated here would be seen by it depends on how the induced error varies
% along the manoeuvre, which this script does not model; it is therefore left
% out of the evaluation rather than claimed either way.
% ---------------------------------------------------------------------------

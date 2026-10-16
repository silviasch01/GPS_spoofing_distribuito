% relative_distance_detection.m
%
% Experimental take-down of classical MULTI-RECEIVER RELATIVE-DISTANCE
% spoofing detection, applied to the distributed-spoofing simulator in
% GPS_spoofing_sim_20250617.m.
%
% ------------------------------------------------------------------ detectors
% A swarm of m receivers flies in a known rigid formation. Each receiver u
% computes its GNSS fix p_hat(u); the controller knows the commanded formation,
% hence the nominal mutual distances d0_ij. All the classical detectors work on
% the normalized ratios  r_ij = d_hat_ij / d0_ij :
%
%   Jansen et al. (ACSAC 2016) -- COLLAPSE test. A single-transmitter spoofer
%     forces all receivers onto the same position, so the GNSS distances
%     shrink (r_ij -> 0). Alarm when f(r) < lambda, with f one of:
%        minimal   : f = min_ij r_ij
%        maximal   : f = max_ij r_ij
%        majority  : f = median_ij r_ij   (>=half of the pairs below lambda)
%     Jansen reports the (normalized) majority voting as the best choice; the
%     three variants are all implemented here.
%   Swaszek et al. (2013) / Borio & Gioia (2022) -- two-sided CONSISTENCY test,
%     same family. Alarm when the distances deviate (either way) from the known
%     formation:  T = max_ij |r_ij - 1| > lambda.
%
% SemperFi (Sathaye et al., NDSS 2022), GNSS/INS fusion, is NOT evaluated here.
% It is cited as blind for the SAME structural reason: a rigid translation of
% the whole swarm preserves the SHAPE of the trajectory flown during the
% anti-spoofing maneuver, so the IMU-vs-GPS track comparison (which only sees
% relative motion, not absolute position) falls in the same translational
% null-space. See the note at the bottom of this file.
%
% --------------------------------------------------------------- attack models
%   H0     : no spoofing                 (used only to calibrate the thresholds)
%   rigid  : null-space attack. A coordinated distributed spoofer imposes the
%            SAME position offset Delta on every receiver (spoofMode='rigid').
%            The absolute position is wrong, but every r_ij ~ 1, so the whole
%            detector family is BLIND. --> this is the result we want to show.
%   power  : naive distributed attack with the physical ground-spoofer power
%            model. It distorts the formation (r_ij spread), so a consistency
%            test DOES fire. Shown as contrast: the attacker must be coordinated
%            (rigid) to defeat the detectors.
%
% The thresholds are calibrated on H0 (same GNSS/ranging noise) for a target
% false alarm probability Pfa; the detection probability Pd is the fraction of
% spoofed epochs that raise an alarm. The formation spacing (grid_step) is
% swept, because the formation distortion of the NAIVE attack scales with the
% baselines while the GNSS noise does not.

clear
close all

% ---------------------------------------------------------------- parameters
grid_steps = [5 10 20 50 100];   % formation spacing [m]
grid_points = 5;                 % 5x5 receivers (m = 25)
Nruns = 2;                       % Monte Carlo runs (random start time)
sigma_pr = 3;                    % pseudorange noise std [m]
sigma_range = 0.1;               % inter-receiver ranging noise std [m]
Pfa = 0.01;                      % target false alarm probability
roc_step = 20;                   % grid_step used for the ROC figure

scen(1) = struct('name','H0',    'opts',struct('Pt_spoof',0));
scen(2) = struct('name','rigid', 'opts',struct('spoofMode','rigid'));
scen(3) = struct('name','power', 'opts',struct('spoofMode','power','max_spoofed_sats',inf));
det_names = {'Jansen min','Jansen max','Jansen majority','Swaszek/Borio consistency'};

Nscen = numel(scen);
Ndet = numel(det_names);
Ngs = numel(grid_steps);

% ------------------------------------------------------------- simulations
T = cell(Ngs,Nscen);   % test statistics, one row per epoch
E = cell(Ngs,Nscen);   % median relative distance error per epoch
for g = 1:Ngs
	for s = 1:Nscen
		T{g,s} = [];
		E{g,s} = [];
		for r = 1:Nruns
			% same random realization (start time, clocks, noise) across the
			% scenarios of a given run/spacing
			seed = 1000*r + g;
			rng(seed);
			opts = scen(s).opts;
			opts.PLOT = false;
			opts.grid_step = grid_steps(g);
			opts.grid_points = grid_points;
			opts.sigma_pr = sigma_pr;
			out = GPS_spoofing_sim_20250617(opts);

			[Tr,Er] = distance_tests(out,sigma_range);
			T{g,s} = [T{g,s}; Tr];
			E{g,s} = [E{g,s}; Er];
			fprintf('grid_step %4g m | %-6s | run %d | median rel. dist. error %6.2f %%\n',...
				grid_steps(g),scen(s).name,r,100*median(Er));
		end
	end
end

% -------------------------------------------- thresholds and Pd at given Pfa
% Collapse detectors (Jansen) raise the alarm for SMALL statistics (f < lambda);
% the consistency detector raises it for LARGE statistics (T > lambda).
alarm_small = [true true true false];
Pd = zeros(Ngs,Nscen,Ndet);
for g = 1:Ngs
	for k = 1:Ndet
		t0 = T{g,1}(:,k); t0 = t0(~isnan(t0));
		if alarm_small(k)
			thr = lower_quantile(t0,Pfa);
			for s = 1:Nscen
				tk = T{g,s}(:,k); Pd(g,s,k) = mean(tk < thr);
			end
		else
			thr = upper_quantile(t0,Pfa);
			for s = 1:Nscen
				tk = T{g,s}(:,k); Pd(g,s,k) = mean(tk > thr);
			end
		end
	end
end

fprintf('\nDetection probability at Pfa = %g (sigma_pr = %g m, sigma_range = %g m)\n',...
	Pfa,sigma_pr,sigma_range);
for k = 1:Ndet
	fprintf('\n%s\n  grid_step [m] ',det_names{k});
	fprintf('%10s',scen.name);
	fprintf('\n');
	for g = 1:Ngs
		fprintf('  %13g ',grid_steps(g));
		fprintf('%10.3f',Pd(g,:,k));
		fprintf('\n');
	end
end

fprintf('\nMedian relative distance error |d_hat-d|/d [%%]\n  grid_step [m] ');
fprintf('%10s',scen.name);
fprintf('\n');
for g = 1:Ngs
	fprintf('  %13g ',grid_steps(g));
	for s = 1:Nscen
		fprintf('%10.2f',100*median(E{g,s}));
	end
	fprintf('\n');
end

% ------------------------------------------------------------------ figures
colors = lines(Nscen);
figure(1)
for k = 1:Ndet
	subplot(2,2,k)
	for s = 2:Nscen
		semilogx(grid_steps,Pd(:,s,k),'-o','Color',colors(s,:),'LineWidth',1.5);
		hold on
	end
	semilogx(grid_steps,Pfa*ones(1,Ngs),'k--');
	grid on; ylim([0 1.05])
	xlabel('formation spacing [m]'); ylabel('P_D')
	title(sprintf('%s (P_{FA}=%g)',det_names{k},Pfa))
	if k==1, legend({scen(2:end).name,'P_{FA}'},'Location','best'); end
end

g = find(grid_steps==roc_step,1);
if ~isempty(g)
	figure(2)
	for k = 1:Ndet
		subplot(2,2,k)
		t0 = sort(T{g,1}(~isnan(T{g,1}(:,k)),k));
		thr_c = [-inf;t0;inf];
		for s = 2:Nscen
			tk = T{g,s}(:,k); tk = tk(~isnan(tk));
			pfa_c = zeros(numel(thr_c),1); pd_c = pfa_c;
			for i = 1:numel(thr_c)
				if alarm_small(k)
					pfa_c(i) = mean(t0 < thr_c(i)); pd_c(i) = mean(tk < thr_c(i));
				else
					pfa_c(i) = mean(t0 > thr_c(i)); pd_c(i) = mean(tk > thr_c(i));
				end
			end
			plot(pfa_c,pd_c,'Color',colors(s,:),'LineWidth',1.5); hold on
		end
		plot([0 1],[0 1],'k--'); grid on; axis square
		xlabel('P_{FA}'); ylabel('P_D')
		title(sprintf('ROC %s, spacing %g m',det_names{k},roc_step))
		if k==1, legend({scen(2:end).name,'chance'},'Location','southeast'); end
	end
end

figure(3)
for s = 2:Nscen
	semilogx(grid_steps,cellfun(@(e) 100*median(e),E(:,s)),'-o',...
		'Color',colors(s,:),'LineWidth',1.5); hold on
end
grid on
xlabel('formation spacing [m]'); ylabel('median |d_{hat}-d|/d [%]')
legend({scen(2:end).name},'Location','best')
title('Distortion of the relative distances')


% ======================================================== local functions

% Jansen collapse variants + Swaszek/Borio consistency, per epoch.
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
		d0_clean = sqrt(sum((p_true(I,:)-p_true(J,:)).^2,2));  % known formation
		d0 = d0_clean + sigma_range*randn(size(d0_clean));     % ranging noise
		d_gnss = sqrt(sum((p_est(I,:)-p_est(J,:)).^2,2));      % GNSS distances
		r = d_gnss ./ d0;                                      % normalized ratio
		Tr(t,1) = min(r);
		Tr(t,2) = max(r);
		Tr(t,3) = median(r);
		Tr(t,4) = max(abs(r-1));
		Er(t) = median(abs(d_gnss-d0_clean)./d0_clean);
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

% ---------------------------------------------------------------------------
% NOTE on SemperFi (GNSS/INS fusion, Sathaye et al., NDSS 2022) -- not tested.
%
% SemperFi commands a pseudorandom maneuver and compares the trajectory tracked
% by the IMU (dead-reckoning) against the one reported by GPS, raising an alarm
% when the averaged Euclidean deviation exceeds ~5 m. It detects attacks that
% make the GPS track DISAGREE with the inertial track in shape/relative motion.
%
% The rigid-translation attack adds the SAME constant offset Delta to every GPS
% fix at every epoch. The GPS-derived displacement between consecutive samples,
% (p+Delta)[t+1] - (p+Delta)[t] = p[t+1]-p[t], is therefore identical to the
% inertial one: the maneuver shape is preserved and the IMU-vs-GPS deviation
% stays ~0. The attack lives in the same translational null-space as the
% relative-distance detectors, so SemperFi is cited as blind for the same
% structural reason, without a separate experiment.
% ---------------------------------------------------------------------------

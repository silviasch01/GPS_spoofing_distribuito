% relative_distance_detection.m
%
% Classical GNSS spoofing detection based on RELATIVE DISTANCES within a
% swarm/formation of receivers, applied to the simulator in
% GPS_spoofing_sim_20250617.m.
%
% Every receiver u computes its GNSS fix p_hat(u). Independently, the
% receivers measure their true mutual distances d_ij (e.g. UWB/radio ranging
% or a known rigid formation). Classical detectors compare the two:
%
%   D1 "collapse" (Tippenhauer et al., CCS 2011; Jansen et al., Crowd-GPS-Sec
%      2018): a single spoofer makes all receivers compute (almost) the same
%      position, so the GNSS distances shrink:
%         T1 = 1 - sum(d_hat_ij)/sum(d_ij)          (alarm if T1 > gamma)
%   D2 per-epoch distance consistency (chi-square / RSS residual test):
%         T2 = sum((d_hat_ij - d_ij).^2)
%   D3 worst-pair test:
%         T3 = max|d_hat_ij - d_ij|
%   D4 time-integrated version of D2 (moving average over Wint epochs).
%
% Thresholds are set empirically from the no-spoofing scenario (H0, with the
% same GNSS/ranging noise) to obtain a target false alarm probability Pfa;
% the detection probability Pd is the fraction of spoofed epochs above the
% threshold.
%
% Scenarios (same trajectory/satellites/clock errors for each run):
%   H0          : no spoofing (used only to set the thresholds for Pfa)
%   distributed : distributed spoofing (24 spoofers on a circle), all the
%                 satellites spoofed
%
% The formation spacing (grid_step) is swept: the formation distortion caused
% by the distributed attack is proportional to the baselines, while the GNSS
% noise is not, so the result depends on the ratio between the two.

clear
close all

% ---------------------------------------------------------------- parameters
grid_steps = [5 10 20 50 100];   % formation spacing [m]
grid_points = 5;                 % 5x5 receivers
Nruns = 2;                       % Monte Carlo runs (random start time)
sigma_pr = 3;                    % pseudorange noise std [m] (uncorrelated part)
sigma_range = 0.1;               % inter-receiver ranging noise std [m]
Pfa = 0.01;                      % target false alarm probability
Wint = 12;                       % D4 integration window [epochs] (12*5 s = 1 min)
roc_step = 20;                   % grid_step used for the ROC figure

scen(1) = struct('name','H0',          'opts',struct('Pt_spoof',0));
scen(2) = struct('name','distributed', 'opts',struct('N_spoofers',24,'max_spoofed_sats',inf));
det_names = {'D1 collapse','D2 RSS','D3 max pair','D4 RSS integrated'};

Nscen = numel(scen);
Ndet = numel(det_names);
Ngs = numel(grid_steps);

% ------------------------------------------------------------- simulations
% T{g,s}: (#epochs, Ndet) test statistics; E{g,s}: relative distance error
T = cell(Ngs,Nscen);
E = cell(Ngs,Nscen);
for g = 1:Ngs
	for s = 1:Nscen
		T{g,s} = [];
		E{g,s} = [];
		for r = 1:Nruns
			% same random realization (start time, clocks, noise) for all
			% the scenarios of a given run
			seed = 1000*r + g;
			rng(seed);
			opts = scen(s).opts;
			opts.PLOT = false;
			opts.grid_step = grid_steps(g);
			opts.grid_points = grid_points;
			opts.sigma_pr = sigma_pr;
			out = GPS_spoofing_sim_20250617(opts);

			[Tr,Er] = distance_tests(out,sigma_range,Wint);
			T{g,s} = [T{g,s}; Tr];
			E{g,s} = [E{g,s}; Er];
			fprintf('grid_step %4g m | %-10s | run %d | median rel. dist. error %6.2f %%\n',...
				grid_steps(g),scen(s).name,r,100*median(Er));
		end
	end
end

% -------------------------------------------- thresholds and Pd at given Pfa
Pd = zeros(Ngs,Nscen,Ndet);
for g = 1:Ngs
	for k = 1:Ndet
		t0 = T{g,1}(:,k);
		t0 = t0(~isnan(t0));
		thr = upper_quantile(t0,Pfa);
		for s = 1:Nscen
			tk = T{g,s}(:,k);
			tk = tk(~isnan(tk));
			Pd(g,s,k) = mean(tk > thr);
		end
	end
end

fprintf('\nDetection probability at Pfa = %g (sigma_pr = %g m, sigma_range = %g m)\n',...
	Pfa,sigma_pr,sigma_range);
for k = 1:Ndet
	fprintf('\n%s\n  grid_step [m] ',det_names{k});
	fprintf('%12s',scen.name);
	fprintf('\n');
	for g = 1:Ngs
		fprintf('  %13g ',grid_steps(g));
		fprintf('%12.3f',Pd(g,:,k));
		fprintf('\n');
	end
end

fprintf('\nMedian relative distance error |d_hat-d|/d [%%]\n  grid_step [m] ');
fprintf('%12s',scen.name);
fprintf('\n');
for g = 1:Ngs
	fprintf('  %13g ',grid_steps(g));
	for s = 1:Nscen
		fprintf('%12.2f',100*median(E{g,s}));
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
	grid on
	ylim([0 1.05])
	xlabel('formation spacing [m]')
	ylabel('P_D')
	title(sprintf('%s (P_{FA}=%g)',det_names{k},Pfa))
	if k==1
		legend({scen(2:end).name,'P_{FA}'},'Location','best')
	end
end

g = find(grid_steps==roc_step,1);
if ~isempty(g)
	figure(2)
	for k = 1:Ndet
		subplot(2,2,k)
		t0 = T{g,1}(:,k);
		t0 = sort(t0(~isnan(t0)));
		for s = 2:Nscen
			tk = T{g,s}(:,k);
			tk = tk(~isnan(tk));
			pfa_c = zeros(numel(t0)+1,1);
			pd_c = zeros(numel(t0)+1,1);
			thr_c = [t0;inf];
			for i = 1:numel(thr_c)
				pfa_c(i) = mean(t0 >= thr_c(i));
				pd_c(i) = mean(tk >= thr_c(i));
			end
			plot(pfa_c,pd_c,'Color',colors(s,:),'LineWidth',1.5);
			hold on
		end
		plot([0 1],[0 1],'k--')
		grid on
		axis square
		xlabel('P_{FA}')
		ylabel('P_D')
		title(sprintf('ROC %s, spacing %g m',det_names{k},roc_step))
		if k==1
			legend({scen(2:end).name,'chance'},'Location','southeast')
		end
	end
end

figure(3)
for s = 1:Nscen
	semilogx(grid_steps,cellfun(@(e) 100*median(e),E(:,s)),'-o',...
		'Color',colors(s,:),'LineWidth',1.5);
	hold on
end
grid on
xlabel('formation spacing [m]')
ylabel('median |d_{hat}-d|/d [%]')
legend({scen.name},'Location','best')
title('Distortion of the relative distances')


% ======================================================== local functions

% test statistics for every epoch of one simulation
function [Tr,Er] = distance_tests(out,sigma_range,Wint)
	Nt = size(out.coords_est,1);
	Nusers = size(out.coords_est,2)/3;
	[I,J] = find(triu(true(Nusers),1));
	Tr = zeros(Nt,4);
	Er = zeros(Nt,1);
	for t = 1:Nt
		p_true = coords_to_ecef(out.coords_grid(t,2:end));
		p_est = coords_to_ecef(out.coords_est(t,:));
		d_true = sqrt(sum((p_true(I,:)-p_true(J,:)).^2,2));
		d_meas = d_true + sigma_range*randn(size(d_true)); % ranging
		d_gnss = sqrt(sum((p_est(I,:)-p_est(J,:)).^2,2));
		res = d_gnss - d_meas;
		Tr(t,1) = 1 - sum(d_gnss)/sum(d_meas);
		Tr(t,2) = sum(res.^2);
		Tr(t,3) = max(abs(res));
		Er(t) = median(abs(d_gnss-d_true)./d_true);
	end
	% D4: moving average of D2 over the last Wint epochs
	Tr(:,4) = NaN;
	for t = Wint:Nt
		Tr(t,4) = mean(Tr(t-Wint+1:t,2));
	end
end

% [lat1,lon1,alt1,lat2,...] -> (#users,3) Cartesian positions. The same
% spherical Earth model of the simulator is used; the Earth rotation term is
% common to all the receivers and does not change their distances.
function p = coords_to_ecef(c)
	R_Earth = physconst('EarthRadius');
	lat = c(1:3:end).'*pi/180;
	lon = c(2:3:end).'*pi/180;
	alt = c(3:3:end).';
	p = (R_Earth+alt).*[cos(lat).*cos(lon), cos(lat).*sin(lon), sin(lat)];
end

% threshold such that a fraction Pfa of x is above it (no toolbox needed)
function thr = upper_quantile(x,Pfa)
	x = sort(x(:));
	n = numel(x);
	thr = x(min(n,max(1,ceil((1-Pfa)*n))));
end

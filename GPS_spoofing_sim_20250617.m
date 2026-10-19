function out = GPS_mobility_v4_prove(opts)
	% opts (optional struct) overrides the default parameters below, e.g.
	%   out = GPS_spoofing_sim_20250617(struct('PLOT',false,'N_spoofers',1));
	% out returns the true/estimated positions so that other scripts (e.g.
	% relative_distance_detection.m) can run detection tests on them.
	if nargin < 1
		opts = struct();
	end
	PLOT = getopt(opts,'PLOT',true);

	% 1: Time Difference of Arrival
	% 2: Linearized Gradient Descent
	fixMethod = getopt(opts,'fixMethod',2);
	KALMAN_FILTER = getopt(opts,'KALMAN_FILTER',false); % smoothing?
	CONSENSUS = getopt(opts,'CONSENSUS',false);

	% list of [time(s), lat(°), lon(°), alt(m)];
	WP = [0.0,  38.1, 13.3, 100.0; ...
		300.0,  38.15, 13.3, 500.0; ...
		600.0,  38.15, 13.35, 300.0; ...
		900.0,  38.1, 13.35, 700.0; ...
		1200.0, 38.1, 13.3, 100.0];
	Ts = 5; % s
	% WP(:,1) = WP(:,1)*86400/1200;
	% Ts = Ts*86400/1200;
	WP(:,1) = WP(:,1)+86400*log(1/rand);

	% returns the list of [time, lat, lon, alt] for every point along the path,
	% excluding the last
	coords = WP_to_coords(WP,Ts);
	Nt = size(coords,1);

	% expand the list of positions to a square (on the sphere) grid of positions
	grid_step = getopt(opts,'grid_step',100);
	% grid_step = 2000;

	grid_points = getopt(opts,'grid_points',5);
	coords_grid = grid_expand(coords,grid_points,grid_step);
	Nusers = grid_points^2;
	u_center = ceil(grid_points/2)+(ceil(grid_points/2)-1)*grid_points;
	u_color = jet(grid_points^2);

	% generate the true delays and elevations for every point in the coords
	% list/grid
	% returns matrices with size (#times, #satellites, #users)
	[delay,azimuth,elevation,pt_sat] = gps_trace(coords_grid);
	Nsats = size(delay,2);

	% waypoint model for the spoofed positions
	% - the center of the spoofed area (center_coords) will measure a set of
	%   delays resulting in the sequence (coords_spoof) of spoofed positions
	% - other points in the affected area will measure different delays (due to
	%   the distributed spoofing method) and will estimate a position close to
	%   it, but different
	WP_spoof = [...
		0.0,    38.120, 13.360, 100.0; ...
		1200.0, 38.120, 13.360, 100.0];
	WP_spoof(:,1) = WP_spoof(:,1) + WP(1,1);
	coords_spoof = WP_to_coords(WP_spoof,Ts);
	[delay_spoof,azimuth_spoof,elevation_spoof] = gps_trace(coords_spoof);

	% generate the (static) spoofer's positions (N_spoofers on a
	% circumference centered around the spoofed area center_coords, with
	% radius R_spoofers)
	% N_spoofers = 24;
	N_spoofers = getopt(opts,'N_spoofers',24);
	center_coords = [38.1, 13.3, 100];
	if N_spoofers > 1
		R_spoofers = getopt(opts,'R_spoofers',2500);
	else
		R_spoofers = 0;
	end
	spoofer_coords = points_on_circle(N_spoofers,center_coords,R_spoofers);

	% determine the best (closest) spoofer for every satellite
	spoof_for_sat = spoofer_assign(center_coords,spoofer_coords,azimuth_spoof);

	% returns a matrix of size (#times, #satellites, #users) with the measured
	% delays, affected by timing error delta_t(user).
	% - out-of-sight satellites have a infinite delay
	Gt_spoof = 10^(1/10);
	if CONSENSUS
		spoofed_ratio = 0.5;
		spoofed_sats_rel_power = [1]; % 1: the strongest, 2: the second strongest, ... -1: the last, -2: second to last
		[measured_delay,delta_t,Pt_spoof_dBW] = spoofing_sim_consensus(...
			spoofed_ratio, spoofed_sats_rel_power, Gt_spoof, spoofer_coords,...
			coords_grid, elevation, delay, delay_spoof, R_spoofers);
		currentfigure(2)
		hold off
		plot(Pt_spoof_dBW(:,:,1)+30)
		hold on
		plot(Pt_spoof_dBW(:,:,2)+30,'--')
		plot(Pt_spoof_dBW(:,:,3)+30,'--')
		xlabel('time index')
		ylabel('Pt_spoof [dBm]')
		grid on
	elseif strcmpi(getopt(opts,'spoofMode','power'),'rigid')
		% Null-space (rigid-translation) attack: an idealized coordinated
		% distributed spoofer that knows the swarm geometry and imposes the
		% SAME position offset Delta on every receiver, for every satellite.
		% The absolute position is wrong by Delta, but all the relative
		% distances between receivers are preserved exactly (up to pseudorange
		% noise). This is the attack that lies in the null-space of any
		% relative-distance / rigid-formation detector (e.g. Jansen 2016).
		[measured_delay,delta_t,is_spoofed] = spoofing_sim_rigid(...
			coords,coords_spoof,coords_grid,elevation,pt_sat);
	else
		Pt_spoof = getopt(opts,'Pt_spoof',1e-3);
		% maximum number of satellites replaced by the spoofer for every user
		max_spoofed_sats = getopt(opts,'max_spoofed_sats',inf);
		[measured_delay,delta_t,is_spoofed] = spoofing_sim(Pt_spoof,Gt_spoof,spoofer_coords,...
			spoof_for_sat,coords_grid,elevation,delay,delay_spoof,R_spoofers,max_spoofed_sats);
	end

	% pseudorange noise (std in meters), independent for every measurement
	sigma_pr = getopt(opts,'sigma_pr',0);
	if sigma_pr > 0
		measured_delay = measured_delay + ...
			sigma_pr/physconst('LightSpeed')*randn(size(measured_delay));
	end

	% estimate the positions
	% - fixMethod: TDoA or Gradient Descent
	% - KALMAN_FILTER: KF smoothing
	coords_est = coords_estimate(fixMethod,KALMAN_FILTER,measured_delay,pt_sat);

	if nargout > 0
		out.coords = coords;
		out.coords_grid = coords_grid; % [time, lat1, lon1, alt1, lat2, ...]
		out.coords_est = coords_est;   % [lat1, lon1, alt1, lat2, ...]
		out.coords_spoof = coords_spoof;
		out.center_coords = center_coords;
		out.spoofer_coords = spoofer_coords;
		out.delta_t = delta_t;
		if exist('is_spoofed','var')
			out.is_spoofed = is_spoofed; % (#times, #satellites, #users)
		end
		out.elevation = elevation;
	end
	if ~PLOT
		return
	end

	% plot
	dxlim = 2*max(max(coords_grid(:,3:3:end)-coords(:,3))); % lon
	dylim = 2*max(max(coords_grid(:,2:3:end)-coords(:,2))); % lat
	Nt = size(coords_grid,1);
	currentfigure(1)
	for t = 1:Nt
		hold off
		for u = 1:Nusers
			plot(coords_grid(t,3+3*(u-1)),... % lon
				coords_grid(t,2+3*(u-1)), ... % lat
				'o','Color',u_color(u,:),'MarkerFaceColor',u_color(u,:)),
			if u==1
				hold on
			end
			plot(coords_est(t,2+3*(u-1)),... % lon
				coords_est(t,1+3*(u-1)), ... % lat
				'square','MarkerFaceColor','w','MarkerEdgeColor',u_color(u,:)),
		end
		% xlim(coords(t,3)+dxlim*[-1,1]); % lon
		% ylim(coords(t,2)+dylim*[-1,1]); % lat
		ind = find(elevation(t,:,u_center)>0);
		plot(center_coords(2),center_coords(1),'rx','MarkerSize',10);
		plot(coords_spoof(t,3),coords_spoof(t,2),'ro','MarkerSize',10);
		for i = 1:N_spoofers
			plot(spoofer_coords(2+3*(i-1)),spoofer_coords(1+3*(i-1)),'^r')
		end
		plot(coords(t,3)+dxlim*[0.8;1]*cos(azimuth(t,ind,u_center)),...
			coords(t,2)+dylim*[0.8;1]*sin(azimuth(t,ind,u_center)),...
			'-k','LineWidth',2);
		grid on
		drawnow
	end
end


function coords_est = coords_estimate(fixMethod,KALMAN_FILTER,measured_delay,pt_sat)
	c = physconst('LightSpeed');
	R_Earth = physconst('EarthRadius');
	T_Earth = 86164; % [s], rotation period with respect to stars
	theta0_Earth = 0.23*2*pi; % initial phase (rotation)

	Nt = size(measured_delay,1);
	Nsats = size(measured_delay,2);
	Nusers = size(measured_delay,3);
	Ts = diff(pt_sat(4,1,1:2));

	As = [1,Ts;0,1];
	q_lat  = 1e-6; % piccola se il sistema è lento, 0.1/Ts^2?
	q_lon  = 1e-6;
	q_alt  = 1e-6;
	q_time = 1e-6;
	% Rw = [Ts^2/4,Ts/2;Ts/2,1]; % da moltiplicare per la q appropriata
	Rw = [Ts^3/3,Ts^2/2;Ts^2/2,Ts];
	% modello dell'osservazione:
	% x_meas = C*[x;x'] + v
	C = [1 0];
	sigma_v_lat = 0.01/R_Earth*180/pi;
	sigma_v_lon = 0.01/R_Earth*180/pi;
	sigma_v_alt = 100;
	sigma_v_time = 1e-4;

	coords_est = zeros(Nt,3*Nusers);

	for u = 1:Nusers
		lat_est = randn(2, 1);
		Plat_est = eye(2);
		lon_est = randn(2, 1);
		Plon_est = eye(2);
		alt_est = randn(2, 1);
		Palt_est = eye(2);
		time_est = randn(2, 1);
		Ptime_est = eye(2);

		for t = 1:Nt
			visible_sats = find(~isinf(measured_delay(t,:,u)));
			Nvisible_sats = numel(visible_sats);
			if fixMethod==1
				% we choose the satellite with minimum measured delay as
				% reference satellite
				[~,sat_ref] = min(measured_delay(t,visible_sats,u));

				A = zeros(Nvisible_sats - 1, 4);
				b = zeros(Nvisible_sats - 1, 1);
				k = 0;
				for sat = 1:Nvisible_sats
					if sat ~= sat_ref
						k = k + 1;
						xi = pt_sat(1, visible_sats(sat), t);
						xj = pt_sat(1, visible_sats(sat_ref), t);
						yi = pt_sat(2, visible_sats(sat), t);
						yj = pt_sat(2, visible_sats(sat_ref), t);
						zi = pt_sat(3, visible_sats(sat), t);
						zj = pt_sat(3, visible_sats(sat_ref), t);
						ti = measured_delay(t, visible_sats(sat), u);
						tj = measured_delay(t, visible_sats(sat_ref), u);
						A(k, :) = [-2 * (xi - xj), -2 * (yi - yj), -2 * (zi - zj), 2 * c^2 * (ti - tj)];
						b(k) = xj^2 - xi^2 + yj^2 - yi^2 + zj^2 - zi^2 + c^2 * (ti^2 - tj^2);
					end
				end

				% Least squares estimate
				pt_est = A \ b;

			elseif fixMethod==2

				% Linearized Gradient Descent
				X_total = 0;
				Y_total = 0;
				Z_total = 0;
				Delta_t = 0;
				R_total = zeros(Nvisible_sats,1);
				PSR = zeros(Nvisible_sats,1);
				matrix = zeros(Nvisible_sats,4);
				iteration = 5;
				for idx = 1:iteration
					for sat = 1:Nvisible_sats
						sat_x = pt_sat(1, visible_sats(sat), t);
						sat_y = pt_sat(2, visible_sats(sat), t);
						sat_z = pt_sat(3, visible_sats(sat), t);
						R_total(sat) = sqrt((sat_x - X_total)^2 + ...
							(sat_y - Y_total)^2 + (sat_z - Z_total)^2);
						PSR(sat) = (measured_delay(t,visible_sats(sat),u) - Delta_t)*c;
						matrix(sat,:)=[(X_total - sat_x)/R_total(sat),...
							(Y_total - sat_y)/R_total(sat), ...
							(Z_total - sat_z)/R_total(sat), c];
					end

					Delta = matrix \ (PSR-R_total);
					X_total = X_total + Delta(1);
					Y_total = Y_total + Delta(2);
					Z_total = Z_total + Delta(3);
					Delta_t = Delta(4);
				end
				pt_est = [X_total;Y_total;Z_total;Delta_t];
			end

			% convert to lat, lon, alt
			lat_meas = 180/pi*asin(pt_est(3) / norm(pt_est(1:3)));
			lon_meas = 180/pi*(mod(atan2(pt_est(2), pt_est(1)) - 2*pi*pt_sat(4,1,t)/T_Earth-theta0_Earth+pi,2*pi)-pi);
			alt_meas = norm(pt_est(1:3)) - R_Earth;
			time_meas = pt_est(4);

			if KALMAN_FILTER

				% prediction
				lat_est_pri = As*lat_est;
				Plat_est_pri = As*Plat_est*As'+q_lat*Rw;
				lon_est_pri = As*lon_est;
				Plon_est_pri = As*Plon_est*As'+q_lon*Rw;
				alt_est_pri = As*alt_est;
				Palt_est_pri = As*Palt_est*As'+q_alt*Rw;
				time_est_pri = As*time_est;
				Ptime_est_pri = As*Ptime_est*As'+q_time*Rw;

				% update
				Ilat = lat_meas-C*lat_est_pri;
				Slat = C*Plat_est_pri*C'+sigma_v_lat^2;
				Klat = Plat_est_pri*C'/Slat;
				lat_est = lat_est_pri + Klat*Ilat;
				Plat_est = (eye(2)-Klat*C)*Plat_est_pri;

				Ilon = lon_meas-C*lon_est_pri;
				Slon = C*Plon_est_pri*C'+sigma_v_lon^2;
				Klon = Plon_est_pri*C'/Slon;
				lon_est = lon_est_pri + Klon*Ilon;
				Plon_est = (eye(2)-Klon*C)*Plon_est_pri;

				Ialt = alt_meas-C*alt_est_pri;
				Salt = C*Palt_est_pri*C'+sigma_v_alt^2;
				Kalt = Palt_est_pri*C'/Salt;
				alt_est = alt_est_pri + Kalt*Ialt;
				Palt_est = (eye(2)-Kalt*C)*Palt_est_pri;

				Itime = time_meas-C*time_est_pri;
				Stime = C*Ptime_est_pri*C'+sigma_v_time^2;
				Ktime = Ptime_est_pri*C'/Stime;
				time_est = time_est_pri + Ktime*Itime;
				Ptime_est = (eye(2)-Ktime*C)*Ptime_est_pri;

				% errore di posizione
				coords_est(t,(u-1)*3+(1:3)) = [lat_est(1),lon_est(1),alt_est(1)];
			else
				coords_est(t,(u-1)*3+(1:3)) = [lat_meas,lon_meas,alt_meas];
			end
		end
	end
end

% Null-space / rigid-translation attack model.
%
% Every receiver u is pushed by the SAME offset vector Delta(i) (in the
% inertial frame used by gps_trace), equal to the displacement that moves the
% true trajectory center onto the spoofed center. For every visible satellite
% the measured delay is the geometric delay from that satellite to the
% receiver's TRUE position shifted by Delta(i), minus the receiver clock
% error. The spoofers are NOT moved: this only models the delays that an
% idealized coordinated ground-based distributed spoofer would have to impose.
%
% Result: coords_est(u) = true_position(u) + Delta, so all pairwise distances
% between receivers are preserved exactly (up to pseudorange noise), while the
% absolute positions are wrong by Delta. This is invisible to any detector
% that only checks the rigid formation geometry.
function [measured_delay,delta_t,is_spoofed] = spoofing_sim_rigid(...
		coords,coords_spoof,coords_grid,elevation,pt_sat)
	c = physconst('LightSpeed');
	Nt = size(coords_grid,1);
	Nusers = (size(coords_grid,2)-1)/3;
	Nsats = size(elevation,2);

	delta_t = 2*rand(1,Nusers)-1;
	delta_t = repmat(delta_t,[Nt,1]);

	measured_delay = zeros(Nt,Nsats,Nusers);
	is_spoofed = false(Nt,Nsats,Nusers);
	for i = 1:Nt
		t = coords(i,1);
		% common rigid offset in the inertial frame
		p_center_true = coords_to_frame(t,coords(i,2),coords(i,3),coords(i,4));
		p_center_spoof = coords_to_frame(t,coords_spoof(i,2),coords_spoof(i,3),coords_spoof(i,4));
		Delta = p_center_spoof - p_center_true;
		sats = pt_sat(1:3,:,i); % satellite positions, same frame
		for u = 1:Nusers
			p_u = coords_to_frame(t, ...
				coords_grid(i,2+3*(u-1)), ...
				coords_grid(i,3+3*(u-1)), ...
				coords_grid(i,4+3*(u-1)));
			p_shift = p_u + Delta;
			for sat = 1:Nsats
				if elevation(i,sat,u) > 0
					d = norm(sats(:,sat)-p_shift)/c;
					measured_delay(i,sat,u) = d - delta_t(i,u);
					is_spoofed(i,sat,u) = true;
				else
					measured_delay(i,sat,u) = inf;
				end
			end
		end
	end
end

% [lat,lon,alt] at time t -> position in the inertial frame used by gps_trace
function p = coords_to_frame(t,lat,lon,alt)
	R_Earth = physconst('EarthRadius');
	T_Earth = 86164;
	theta0_Earth = 0.23*2*pi;
	la = lat*pi/180;
	lo = lon*pi/180;
	u1 = [cos(la)*cos(2*pi*t/T_Earth+theta0_Earth+lo); ...
		cos(la)*sin(2*pi*t/T_Earth+theta0_Earth+lo); ...
		sin(la)];
	p = (R_Earth+alt)*u1;
end

function [measured_delay,delta_t] = spoofing_sim_freespace(Pt_spoof,Gt_spoof,spoofer_coords,...
		spoof_for_sat,coords_grid,elevation,delay,delay_spoof,R_spoofers)
	Nt = size(coords_grid,1);
	Nusers = (size(coords_grid,2)-1)/3;
	Nsats = size(spoof_for_sat,2);

	% generate the set of receiver clock errors delta_t (#users), independent
	% for every user (positive if the receiver's clock is early)
	delta_t = 2*rand(1,Nusers)-1;
	% no drift, constant for all times
	delta_t = repmat(delta_t,[Nt,1]);

	% for every time, for every visible satellite, for every user:
	% - determine the received power from the spoofer and compare it to the
	%   received power from the actual satellite
	% - if the actual satellite received power is stronger than the spoofed
	%   signal, the measured delay is the actual delay - delta_t
	% - otherwise, if the spoofer power is higher than the power from real
	%   satellite, the measured delay in the area affected by spoofing is:
	%      delay_spoof - delta_t - (R_spoof - spoofer_to_rx_distance)/c
	fc = 1575e6;
	Pt_sat = 25; % W
	Gt_sat = 10^(13.5/10);
	Gr = 10^(5/10);
	c = physconst('LightSpeed');
	measured_delay = zeros(Nt,Nsats,Nusers);
	for i = 1:Nt
		for u = 1:Nusers
			for sat = 1:Nsats
				if elevation(i,sat,u) > 0
					spoofer_to_rx_distance = coords_distance(...
						spoofer_coords((spoof_for_sat(i,sat)-1)*3+(1:3)),...
						coords_grid(i,1+(u-1)*3+(1:3)));
					Pr_sat = Pt_sat*Gt_sat*Gr*(1/(4*pi*fc*delay(i,sat,u)))^2;
					Pr_spoofer = Pt_spoof*Gt_spoof*Gr*(1/(4*pi*fc*spoofer_to_rx_distance/c))^2;
					if Pr_sat > Pr_spoofer
						measured_delay(i,sat,u) = delay(i,sat,u) - delta_t(i,u);
					else
						measured_delay(i,sat,u) = delay_spoof(i,sat)-delta_t(i,u)...
							- (R_spoofers - spoofer_to_rx_distance)/c;
					end
				else
					measured_delay(i,sat,u) = inf;
				end
			end
		end
	end

	% add noise?
	% Tchip = 1/1.023e6;
	% Tr = 300;
	% Fr = 10^(3/10); % cifra di rumore?
	% Pw = kBol*Tr*Fr;
	% CA_SF = 1023;
	% SNR = CA_SF*Pr*Tchip./Pw;
	% sigma_delay = Tchip*sqrt(0.001./SNR);
	% randn_delay = sigma_delay.*randn(Nsat,Nt,Nusers);
	% delay(sigma_delay<inf) = delay(sigma_delay<inf) + randn_delay(sigma_delay<inf);
end

function [measured_delay,delta_t,is_spoofed] = spoofing_sim(Pt_spoof,Gt_spoof,spoofer_coords,...
		spoof_for_sat,coords_grid,elevation,delay,delay_spoof,R_spoofers,max_spoofed_sats)
	if nargin < 10
		max_spoofed_sats = inf;
	end

	Nt = size(coords_grid,1);
	Nusers = (size(coords_grid,2)-1)/3;
	Nsats = size(spoof_for_sat,2);

	% Generate the set of receiver clock errors delta_t (#users), independent
	% for every user (positive if the receiver's clock is early)
	delta_t = 2*rand(1,Nusers)-1;
	% No drift, constant for all times
	delta_t = repmat(delta_t,[Nt,1]);

	alpha = 3.04;
	A = -23.29;
	B = 4.14;
	theta_0 = -3.61;
	eta_0 = 20.70;
	a = -0.41;
	sigma_0 = 5.86;

	fc = 1575e6;
	Pt_sat = 25; % W
	Gt_sat = 10^(13.5/10);
	Gr = 10^(5/10);
	c = physconst('LightSpeed');
	R_Earth = physconst('EarthRadius');

	elevation_drones= drone_elevation_from_spoofers(coords_grid, spoofer_coords,Nusers);
	measured_delay = zeros(Nt,Nsats,Nusers);
	is_spoofed = false(Nt,Nsats,Nusers);

	for i = 1:Nt
		for u = 1:Nusers
			N_spoofed_sats = 0;
			for sat = 1:Nsats
				if elevation(i,sat,u) > 0

					spoofer_to_rx_distance = coords_distance(...
						spoofer_coords((spoof_for_sat(i,sat)-1)*3+(1:3)), ...
						coords_grid(i, 1+(u-1)*3+(1:3)));

					Pr_sat = Pt_sat*Gt_sat*Gr*(1/(4*pi*fc*delay(i,sat,u)))^2;

					theta_rad = elevation_drones(i, spoof_for_sat(i,sat), u);
					theta_deg = min(rad2deg(theta_rad),10); % between -2 and 10°
					% if theta_deg < -2 || theta_deg > 10
					% 	fprintf('!');
					% end

					% sigma_db = a*theta_deg + sigma_0;
					% X_db = sigma_db*randn();

					terrestrial_PL_db = 10*alpha *log10(spoofer_to_rx_distance); % -100@1000m
					excess_PL_db = A*(theta_deg - theta_0)*exp(-(theta_deg - theta_0)/B)+eta_0;

					% PL_db = terrestrial_PL_db + excess_PL_db + X_db;
					PL_db = terrestrial_PL_db + excess_PL_db;
					PL_linear = 10^(PL_db/10);

					Pr_spoofer = Pt_spoof*Gt_spoof*Gr/PL_linear;

					if Pr_sat < Pr_spoofer && N_spoofed_sats < max_spoofed_sats
						measured_delay(i,sat,u) = delay_spoof(i,sat) - delta_t(i,u) ...
							- (R_spoofers - spoofer_to_rx_distance)/c;
						N_spoofed_sats = N_spoofed_sats + 1;
						is_spoofed(i,sat,u) = true;
					else
						measured_delay(i,sat,u) = delay(i,sat,u) - delta_t(i,u);
					end
				else
					measured_delay(i,sat,u) = inf;
				end
			end
		end
	end

	% Optional: add noise
	% Tchip = 1/1.023e6;
	% Tr = 300;
	% Fr = 10^(3/10); % noise figure
	% Pw = kBol*Tr*Fr;
	% CA_SF = 1023;
	% SNR = CA_SF*Pr*Tchip./Pw;
	% sigma_delay = Tchip*sqrt(0.001./SNR);
	% randn_delay = sigma_delay.*randn(Nsat,Nt,Nusers);
	% delay(sigma_delay<inf) = delay(sigma_delay<inf) + randn_delay(sigma_delay<inf);
end

function [measured_delay, delta_t, Pt_spoof_dBW] = spoofing_sim_consensus(...
		spoof_ratio, spoofed_sats_rel_power, Gt_spoof, spoofer_coords, ...
		coords_grid, elevation, delay, delay_spoof, R_spoofers)

	Nt = size(coords_grid,1);
	Nusers = (size(coords_grid,2)-1)/3;
	Nsats = size(elevation,2);

	% assumes the users are on a rectangular grid
	grid_points = sqrt(Nusers);
	u_center = ceil(grid_points/2)+(ceil(grid_points/2)-1)*grid_points;

	% Generate the set of receiver clock errors delta_t (#users), independent
	% for every user (positive if the receiver's clock is early)
	delta_t = 2*rand(1,Nusers)-1;
	% No drift, constant for all times
	delta_t = repmat(delta_t,[Nt,1]);

	% ground-to-drone channel model parameters
	% alpha = 3.04;
	% A = -23.29;
	% B = 4.14;
	% theta_0 = -3.61;
	% eta_0 = 20.70;
	% shadowing is disabled, should be spatially and temporally correlated
	% between drones
	% a = -0.41;
	% sigma_0 = 5.86;

	fc = 1575e6;
	Pt_sat = 25; % W
	Gt_sat = 10^(13.5/10);
	Gr = 10^(5/10);
	c = physconst('LightSpeed');
	% R_Earth = physconst('EarthRadius');

	Pt_sat_dBW = 10*log10(Pt_sat);
	Gt_sat_dB = 10*log10(Gt_sat);
	Gr_dB = 10*log10(Gr);
	Gt_spoof_dB = 10*log10(Gt_spoof);

	% elevation_drones = drone_elevation_from_spoofers(coords_grid,spoofer_coords,Nusers); % rad

	N_spoofers = numel(spoofer_coords)/3;
	if N_spoofers > 1
		fprintf('error: spoofing_sim_consensus expects only one spoofer\n')
		return
	end
	N_spoofed_sats = numel(spoofed_sats_rel_power);

	Pt_spoof_dBW = zeros(Nt, N_spoofed_sats, 3);
	measured_delay = zeros(Nt, Nsats, Nusers);
	for i = 1:Nt

		% determine the set spoofed_sats of spoofed satellites based on the
		% satellites' received power at the center of the receivers' grid and
		% the selection spoofed_sats_rel_power
		Pr_sat_dBW_center = zeros(1,Nsats);
		N_visible_sats = 0;
		for sat = 1:Nsats
			if elevation(i,sat,u_center) > 0
				N_visible_sats = N_visible_sats + 1;
				Pr_sat_dBW_center(sat) = Pt_sat_dBW + Gt_sat_dB + Gr_dB - 20*log10(4*pi*fc*delay(i,sat,u_center));
			else
				Pr_sat_dBW_center(sat) = -inf;
			end
		end
		[~,sat_rel_power_idx]=sort(Pr_sat_dBW_center,'descend');
		spoofed_sats = zeros(1,N_spoofed_sats);
		for sat = 1:N_spoofed_sats
			if spoofed_sats_rel_power(sat) > 0
				spoofed_sats(sat) = sat_rel_power_idx(spoofed_sats_rel_power(sat));
			else
				spoofed_sats(sat) = sat_rel_power_idx(N_visible_sats+1+spoofed_sats_rel_power(sat));
			end
		end

		% determine the PathLoss from the spoofer to every drone
		total_PL_dB = zeros(1,Nusers);
		Pr_sat_dBW = zeros(N_spoofed_sats,Nusers);
		spoofer_to_rx_distance = zeros(1,Nusers);
		for u = 1:Nusers
			spoofer_to_rx_distance(u) = coords_distance(...
				spoofer_coords, coords_grid(i, 1+(u-1)*3+(1:3)));
			% if 0
			% 	theta_rad = elevation_drones(i, 1, u);
			% 	theta_deg = min(rad2deg(theta_rad), 10);
			% 	terrestrial_PL_dB = 10*alpha * log10(spoofer_to_rx_distance(u));
			% 	excess_PL_dB = A*(theta_deg - theta_0)*exp(-(theta_deg - theta_0)/B) + eta_0;
			% 	total_PL_dB(u) = terrestrial_PL_dB + excess_PL_dB - Gt_spoof_dB - Gr_dB;
			% else
			total_PL_dB(u) = 20*log10(4*pi*fc*spoofer_to_rx_distance(u)/c) - Gt_spoof_dB - Gr_dB;
			% end
		end

		% determine the exactly needed spoofer transmitted power to spoof a
		% fraction spoof_ratio of the users
		%
		% Pr(Pt_spoof_dBW - total_PL_dB > Pr_sat_dBW) = spoof_ratio
		for sat = 1:N_spoofed_sats
			for u = 1:Nusers
				Pr_sat_dBW(sat,u) = Pt_sat_dBW + Gt_sat_dB + Gr_dB - 20*log10(4*pi*fc*delay(i,spoofed_sats(sat),u));
			end
			tmp_val = sort(Pr_sat_dBW(sat,:) + total_PL_dB);
			Pt_spoof_dBW(i,sat,1) = mean(tmp_val(floor(Nusers*spoof_ratio)+(0:1)));
			Pt_spoof_dBW(i,sat,2) = tmp_val(1);
			Pt_spoof_dBW(i,sat,3) = tmp_val(end)+0.001;
		end
		is_spoofed = Pt_spoof_dBW(i,:,1).'-total_PL_dB > Pr_sat_dBW;

		for u = 1:Nusers
			for sat = 1:Nsats
				if elevation(i,sat,u) > 0
					if any(sat == spoofed_sats) && is_spoofed(u)
						measured_delay(i,sat,u) = delay_spoof(i,sat) - delta_t(i,u) ...
							- (R_spoofers - spoofer_to_rx_distance(u))/c;
					else
						measured_delay(i,sat,u) = delay(i,sat,u) - delta_t(i,u);
					end
				else
					measured_delay(i,sat,u) = inf;
				end
			end
		end
	end
end

function elevation_drones = drone_elevation_from_spoofers(coords_grid, spoofer_coords,Nusers)
	R_Earth = physconst('EarthRadius');

	T_Earth = 86164;
	theta0_Earth = 0.23*2*pi;
	Nt = size(coords_grid, 1);

	N_spoofers = size(spoofer_coords, 2)/3;
	elevation_drones = zeros(Nt, N_spoofers, Nusers);

	for i = 1:Nt
		t = coords_grid(i, 1);
		for d = 1:Nusers

			lat_drone = pi/180*coords_grid(i, 2 + 3*(d-1));
			lon_drone = pi/180*coords_grid(i, 3 + 3*(d-1));
			alt_drone = coords_grid(i, 4 + 3*(d-1));

			u1_drone = [cos(lat_drone)*cos(2*pi*t/T_Earth + theta0_Earth + lon_drone); ...
				cos(lat_drone)*sin(2*pi*t/T_Earth + theta0_Earth + lon_drone); ...
				sin(lat_drone)];
			%p = (R_Earth+z)*u1
			drone = (R_Earth+alt_drone)*u1_drone;

			for s = 1:N_spoofers
				lat_spoofer = pi/180 * spoofer_coords(1 + 3*(s-1));
				lon_spoofer = pi/180 * spoofer_coords(2 + 3*(s-1));
				alt_spoofer = spoofer_coords(3 + 3*(s-1));


				u1_spoofer = [cos(lat_spoofer)*cos(2*pi*t/T_Earth + theta0_Earth + lon_spoofer); ...
					cos(lat_spoofer)*sin(2*pi*t/T_Earth + theta0_Earth + lon_spoofer); ...
					sin(lat_spoofer)];
				%
				% u1 = [cos(lat)*cos(2*pi*t/T_Earth+theta0_Earth+lon); ...
				%     cos(lat)*sin(2*pi*t/T_Earth+theta0_Earth+lon);sin(lat)];


				uN = [0; 0; 1];
				eN_1 = uN - (u1_spoofer'*uN)*u1_spoofer;
				uN_1 = eN_1/norm(eN_1);
				uE = cross(uN_1, u1_spoofer);
				G_spoofer = [uE, uN_1, u1_spoofer];

				spoofer = (R_Earth+alt_spoofer)*u1_spoofer;
				% distance_p_to_s = sqrt(sum(abs(s-p).^2));
				% p_to_s_versor = (s-p)./distance_p_to_s;
				% p_to_s_local = G'*p_to_s_versor;

				distance = sqrt(sum(abs(drone-spoofer).^2));
				spoofer_to_drone_versor = (drone-spoofer)./distance;
				spoofer_to_drone_local = G_spoofer' * spoofer_to_drone_versor;
				elevation_drones(i, s, d) = asin(spoofer_to_drone_local(3));

			end
		end
	end
end


function spoof_for_sat = spoofer_assign(center_coords,spoofer_coords,azimuth_spoof)
	u0 = coords_to_unitv(center_coords);
	uN = [0; 0; 1];
	eN_0 = uN - (u0'*uN)*u0;
	uN_0 = eN_0/norm(eN_0);
	uE = cross(uN_0,u0);
	G = [uE, uN_0, u0]; % x->East, y->North, z->Space

	% azimuth for every spoofer, from the point of view of center_coords
	N_spoof = numel(spoofer_coords)/3;
	azimuth_spoof_local = zeros(1,N_spoof);
	for spoof = 1:N_spoof
		u_spoof = coords_to_unitv(spoofer_coords(3*(spoof-1)+(1:3)));
		spoof_versor = (u_spoof-u0)/norm(u_spoof-u0);
		spoof_versor_local = G'*spoof_versor;
		azimuth_spoof_local(spoof) = atan2(spoof_versor_local(2),spoof_versor_local(1));
	end

	Nt = size(azimuth_spoof,1);
	Nsat = size(azimuth_spoof,2);
	spoof_for_sat = zeros(Nt,Nsat);
	for i = 1:Nt
		for sat = 1:Nsat
			[~,ind]=min(abs(mod(azimuth_spoof(i,sat)-azimuth_spoof_local+pi,2*pi)-pi));
			spoof_for_sat(i,sat) = ind;
		end
	end
end

function spoofer_coords = points_on_circle(N_spoofers,center_coords,R_spoofers)
	R_Earth = physconst('EarthRadius');
	u0 = coords_to_unitv(center_coords);
	uN = [0; 0; 1];
	eN_0 = uN-(u0'*uN)*u0;
	uN_0 = eN_0/norm(eN_0);
	uE = cross(uN_0,u0);
	alpha = R_spoofers/R_Earth;
	spoofer_coords = zeros(1,3*N_spoofers);
	for i = 1:N_spoofers
		beta = 2*pi*(i-1)/N_spoofers;
		u = cos(alpha)*u0 + ...
			sin(alpha)*(uE*cos(beta) + uN_0*sin(beta));
		spoofer_coords(3*(i-1)+(1:3)) = [unitv_to_coords(u),center_coords(3)];
	end
end

function [delay,azimuth,elevation,pt_sat] = gps_trace(coords,plot_visible_sats)
	if nargin < 2
		plot_visible_sats = false;
	end
	R_Earth = physconst('EarthRadius');
	c = physconst('lightspeed');

	Nt = size(coords,1);
	Nusers = (size(coords,2)-1)/3;
	if size(coords,2) ~= 3*Nusers+1
		fprintf('gps_trace: the coords parameter has a wrong format.\n');
		return
	end

	% Earth parameters
	T_Earth = 86164; % [s], rotation period with respect to stars
	theta0_Earth = 0.23*2*pi; % initial phase (rotation)
	% tilt_Earth_deg = 23.439281; % ° wrt ecliptic
	% tilt_Earth = pi*tilt_Earth_deg/180;
	% Rtilt_Earth = [...
	%     1                0               0; ...
	%     0  cos(tilt_Earth) sin(tilt_Earth); ...
	%     0 -sin(tilt_Earth) cos(tilt_Earth) ];

	% GPS parameters (approximate)
	Rh = 20200e3; % satellite altitudes [m]
	tilt_Sat_deg = 55; % [°] orbital plane inclination
	Norbits = 6; % number of orbits
	d_phi_orbit = 2*pi/Norbits;
	phi0_orbit = 0;
	Nsats_per_orbit = 4; % number of satellites per orbit
	theta0_orbit = [0 30 135 255]*pi/180;
	theta0_sat_orbit = 2*pi/Nsats_per_orbit*(mod(0:Norbits-1,3)/3);

	% standard gravitational parameter of the Earth: Me*G
	% G = 6.674e-11; Newton's gravitational constant [m^3/kg s^2]
	% Me = 5.972e24; Earth's mass [kg]
	mu = 3.986e14; % [m^3/s^2]
	Rsat = R_Earth + Rh;
	Tsat = 2*pi*sqrt(Rsat^3/mu);
	tilt_Sat = tilt_Sat_deg * pi/180;
	theta_sat = zeros(Nsats_per_orbit,Norbits);
	Nsats = Nsats_per_orbit*Norbits;
	for k = 1:Norbits
		for sat = 1:Nsats_per_orbit
			theta_sat(sat,k) = theta0_sat_orbit(k) + theta0_orbit(sat);
		end
	end

	% allocate the output variables
	delay = zeros(Nt,Nsats,Nusers);
	azimuth = zeros(Nt,Nsats,Nusers);
	elevation = zeros(Nt,Nsats,Nusers);
	if nargout > 3
		pt_sat = zeros(4,Nsats,Nt);
	end

	for i = 1:Nt
		t = coords(i,1);
		% relative positions of every satellite, in a frame of reference
		% where z is along the rotation axis
		s = zeros(3,Nsats);
		for k = 1:Norbits
			phi_orbit = (k-1)*d_phi_orbit + phi0_orbit;
			for sat = 1:Nsats_per_orbit
				s(:,sat+(k-1)*Nsats_per_orbit) = ...
					[cos(phi_orbit) sin(phi_orbit) 0; -sin(phi_orbit) cos(phi_orbit) 0; 0 0 1] * ...
					[1 0 0; 0 cos(tilt_Sat) -sin(tilt_Sat); 0 sin(tilt_Sat) cos(tilt_Sat)] * ...
					[Rsat*cos(2*pi*t/Tsat+theta_sat(sat,k)); Rsat*sin(2*pi*t/Tsat+theta_sat(sat,k)); 0];
				% s_log(:,i,j+(k-1)*Nsats_per_orbit) = s(:,j+(k-1)*Nsats_per_orbit);
				% lat_sat = asin(s(3, j+(k-1)*Nsats_per_orbit) / Rsat);
				% lon_sat = atan2(s(2, j+(k-1)*Nsats_per_orbit), s(1, j+(k-1)*Nsats_per_orbit));
				% alt_sat = Rh;
				% Salvataggio delle coordinate geografiche
				% s_geo(:,i,j+(k-1)*Nsats_per_orbit) = [rad2deg(lat_sat); rad2deg(lon_sat); alt_sat];
			end
		end
		if nargout > 3
			pt_sat(1:3,:,i) = s;
			pt_sat(4,:,i) = t;
		end

		% s = Rtilt_Earth * s;
		% s1(:,i,j+(k-1)*Nsats_per_orbit) = s_geo(:,i,j+(k-1)*Nsats_per_orbit);

		for u = 1:Nusers
			lat = pi/180*coords(i,2+3*(u-1));
			lon = pi/180*coords(i,3+3*(u-1));
			z = coords(i,4+3*(u-1));

			u1 = [cos(lat)*cos(2*pi*t/T_Earth+theta0_Earth+lon); ...
				cos(lat)*sin(2*pi*t/T_Earth+theta0_Earth+lon);sin(lat)];
			uN = [0; 0; 1];
			eN_1 = uN - (u1'*uN)*u1;
			uN_1 = eN_1/norm(eN_1);
			uE = cross(uN_1,u1);
			G = [uE, uN_1, u1]; % x->East, y->North, z->Space
			p = (R_Earth+z)*u1;
			% p_log(:,i,u) = p;

			% lat = asin(p(3) / norm(p));
			% lon = mod(atan2(p(2), p(1)) - 2*pi*t(i)/T_Earth-theta0_Earth + pi,2*pi)-pi;
			% alt = norm(p) - R_Earth;
			% coord_list_estimate(u, :,i) = [rad2deg(lat), rad2deg(lon), alt];
			% p1(:,i,u) =  [rad2deg(lat), rad2deg(lon), alt];
			% p = Rtilt_Earth * p;
			% p = [rad2deg(lat), rad2deg(lon), alt];

			% theta_p, elevazione apparente dei satelliti
			distance_p_to_s = sqrt(sum(abs(s-p).^2));
			p_to_s_versor = (s-p)./distance_p_to_s;
			p_to_s_local = G'*p_to_s_versor;
			delay(i,:,u) = distance_p_to_s/c;
			azimuth(i,:,u) = atan2(p_to_s_local(2,:),p_to_s_local(1,:));
			elevation(i,:,u) = asin(p_to_s_local(3,:));
		end
	end
	if plot_visible_sats
		for i=1:Nt
			ind=find(elevation(i,:,1)>0);
			hold off
			plot(cos(azimuth(i,ind,1)).*cos(elevation(i,ind,1)),...
				sin(azimuth(i,ind,1)).*cos(elevation(i,ind,1)),'x');
			hold on
			plot(cos(0:pi/64:2*pi),sin(0:pi/64:2*pi),'k--')
			plot(0.5*cos(0:pi/64:2*pi),0.5*sin(0:pi/64:2*pi),'k--')
			plot(0,0,'ok')
			axis square
			drawnow
		end
	end
end

function coords_grid = grid_expand(coords,grid_points,grid_step,print_d)
	if nargin < 4
		print_d = false;
	end

	R_Earth = physconst('EarthRadius');
	Nt = size(coords,1);
	coords_grid = zeros(Nt,1+3*grid_points^2);
	for t = 1:Nt
		% copy time value and altitude
		coords_grid(t,1) = coords(t,1);
		coords_grid(t,4:3:end) = coords(t,4);

		u1 = coords_to_unitv(coords(t,2:3));
		uN = [0;0;1];
		eN_1 = uN-(u1'*uN)*u1;
		uN_1 = eN_1/norm(eN_1);
		uE = cross(u1,uN_1);
		alpha_step = grid_step/R_Earth;
		for i = 1:grid_points
			alpha = alpha_step*(i-1-(grid_points-1)/2);
			ui = u1*cos(alpha)+uN_1*sin(alpha);
			beta_step = grid_step/R_Earth;
			for j = 1:grid_points
				beta = beta_step*(j-1-(grid_points-1)/2);
				u = ui*cos(beta)+uE*sin(beta);
				coords_grid(t,((j-1)+(i-1)*grid_points)*3+(2:3)) = ...
					unitv_to_coords(u);
			end
		end

		if print_d
			for i = 1:grid_points-1
				for j = 1:grid_points-1
					ofs = ((j-1)+(i-1)*grid_points)*3;
					ofs_i = ((j-1)+i*grid_points)*3;
					d = coords_distance(coords_grid(1,ofs+(2:3)),coords_grid(1,ofs_i+(2:3)));
					fprintf('d(%d,%d)->d(%d,%d): %g\n',i,j,i+1,j,d);
				end
			end
			for i = 1:grid_points-1
				for j = 1:grid_points-1
					ofs = ((j-1)+(i-1)*grid_points)*3;
					ofs_j = (j+(i-1)*grid_points)*3;
					d = coords_distance(coords_grid(1,ofs+(2:3)),coords_grid(1,ofs_j+(2:3)));
					fprintf('d(%d,%d)->d(%d,%d): %g\n',i,j,i,j+1,d);
				end
			end
		end
	end
end

% returns the list of [time, lat, lon, alt] for every point along the path,
% excluding the last
%
% time >= WP(1,1) && time < WP(end,1) && mod(time,Ts)==0
%
function coords = WP_to_coords(WP,Ts,print_v)
	if nargin < 3
		print_v = false;
	end
	if print_v
		R_Earth = physconst('EarthRadius');
	end

	% determine the number of points and allocate the list
	t_ini = min(WP(:,1));
	t_fin = max(WP(:,1));
	i_ini = ceil(t_ini/Ts);
	i_fin = ceil(t_fin/Ts)-1;
	Nt = i_fin-i_ini+1;
	coords = zeros(Nt,4);

	% we use j to index the time interval, for interpolation between WP(j,:) and
	% WP(j+1,:)
	j = 1;
	t1 = WP(j,1);a1 = WP(j,4);
	u1 = coords_to_unitv(WP(j,2:3));
	t2 = WP(j+1,1);a2 = WP(j+1,4);
	u2 = coords_to_unitv(WP(j+1,2:3));
	% angle between the starting and ending unit vector for this segment
	theta12 = acos(u1'*u2);
	if print_v
		v = norm((R_Earth+a1)*u1-(R_Earth+a2)*u2)/(t2-t1);
		fprintf('v: %g\n',v);
	end
	% projection error and unit vector for interpolation between u1 and u2
	if theta12
		e2_1 = u2-cos(theta12)*u1;
		u2_1 = e2_1/norm(e2_1);
	else
		u2_1 = zeros(3,1);
	end

	for i = i_ini:i_fin
		t = i*Ts;

		% if the current time is beyon the end of the interval, skip to the next
		% until WP(j,1) <= t < WP(j+1,1)
		while t >= WP(j+1,1)
			j = j + 1;
			t1 = t2;a1 = a2;u1 = u2;
			t2 = WP(j+1,1);a2 = WP(j+1,4);
			u2 = coords_to_unitv(WP(j+1,2:3));
			theta12 = acos(u1'*u2);
			if print_v
				v = norm((R_Earth+a1)*u1-(R_Earth+a2)*u2)/(t2-t1);
				fprintf('v: %g\n',v);
			end
			if theta12
				e2_1 = u2-cos(theta12)*u1;
				u2_1 = e2_1/norm(e2_1);
			else
				u2_1 = zeros(3,1);
			end
		end
		% interpolate along the geodetic
		coords(i-i_ini+1,1) = t;
		alpha = theta12*(t-t1)/(t2-t1);
		coords(i-i_ini+1,2:3) = unitv_to_coords(u1*cos(alpha)+u2_1*sin(alpha));
		coords(i-i_ini+1,4) = a1+(t-t1)*(a2-a1)/(t2-t1);
	end
end

% convert [lat,lon] to unit vector
function u = coords_to_unitv(coords)
	lat = coords(1)*pi/180;
	lon = coords(2)*pi/180;
	u = [cos(lat)*cos(lon);cos(lat)*sin(lon);sin(lat)];
end

% convert unit vector to [lat,lon]
function coords = unitv_to_coords(u)
	coords=180/pi*[asin(u(3)),atan2(u(2),u(1))];
end

% distance between two points with coordinates coords1 and coords2
function distance = coords_distance(coords1,coords2)
	R_Earth = physconst('EarthRadius');
	distance = R_Earth*acos(coords_to_unitv(coords1)'*coords_to_unitv(coords2));
end

function v = getopt(opts,name,default)
	if isfield(opts,name)
		v = opts.(name);
	else
		v = default;
	end
end

function currentfigure(nf)
	if ishghandle(nf)
		set(0,'CurrentFigure',nf)
	else
		figure(nf)
	end
	set(gcf,'WindowStyle','docked');
end

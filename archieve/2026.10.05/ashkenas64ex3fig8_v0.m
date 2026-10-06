clear; clc; close all; clear functions
%% 1. Linear system G(s)
Kp = 1;             % the pilot gain representing the gradient of the sidestick of the pilot model for the X-15 Soft Glide Landing PIO Incident in 1959
Kp_Yp = 13.68;      % controller gain (here pilot gain has been inherited as an internal parameter of the X-15 Soft Glide Landing PIO Incident in 1959)
M_del_e = 0.537;    % Input signal entering the rate limiter
omega_n = 2.3;      % Natural frequency
zeta_sp = 1.42 / omega_n / 2;   % Damping ratio
num = Kp_Yp*M_del_e*[1 0.82];
den = [1 1.42 omega_n^2 0];
Gs = tf(num,den);

%% 2. Frequency range for G(jw)
w_G = linspace(0.3,3.5,1000);

%% 3. Rate Limiter Element parameters
R  = 15;           % Rate limit, deg/s
Ai = 13.68;        % Input amplitude, deg
w_N = linspace(0.01, 10, 1001);   % Frequencies, rad/s
w_onset = R / Ai;

%% 4. Describing function N(Ai,w)
N = zeros(size(w_N));
for k = 1:length(w_N)
    alpha = w_N(k) / w_onset;
    % Region I
    if alpha <= 1
        M = 1;
        phi = 0;
    % Region II
    elseif alpha < 1.862
        M = 0.2908*alpha^3 - 1.4396*alpha^2 + 1.9232*alpha + 0.223;
        phi = 0.528*alpha^3 - 2.6213*alpha^2 + 3.5056*alpha - 1.4171;
    % Region III
    else
        varpi = w_onset / w_N(k);
        M = (4/pi)*varpi;
        phi = -acos((pi/2)*varpi);
    end
    N(k) = M*exp(1j*phi);
end

%% 5. Calculate -1/N(Ai,w)
minus_inv_N = -1 ./ N;

%% 6. Convert -1/N to an FRD model
response = reshape(minus_inv_N,1,1,[]);
sys_N = frd(response,w_N);

%% 7. Nichols Chart
figure(Name='Nichols Chart',NumberTitle='off');
p1 = nicholsplot(Gs,w_G);
hold on
p2 = nicholsplot(sys_N);

%% 8. Shift -1/N from +180 deg to -180 deg
p2.PhaseMatchingEnabled = 'on';
p2.PhaseMatchingFrequency = w_N(1);
phase_first = rad2deg(angle(minus_inv_N(1)));
p2.PhaseMatchingValue = phase_first - 360;

%% 9. Figure settings
grid on
xlim([-180 -50])
ylim([-5 20])
yline(0,'--','DisplayName','0 dB');
title('Nichols Chart of G(j\omega) and -1/N(A_i,\omega)')
legend('G(j\omega)','-1/N(A_i,\omega)','0 dB')

%% 10. Find Intersections (Simple Interpolation Method)
% Get numerical data from the linear system plot
[mag_G, phase_G, w_G_out] = nichols(Gs, w_G);
phase_G = squeeze(phase_G)';
mag_G_dB = 20*log10(squeeze(mag_G))';
w_G_out = squeeze(w_G_out)';

% Get numerical data for the describing function (shifted by -360 to match plot)
phase_N = rad2deg(unwrap(angle(minus_inv_N))) - 360; 
mag_N_dB = 20*log10(abs(minus_inv_N));

% Step A: Create a common X-axis (Phase) where both curves exist
min_phase = max(min(phase_G), min(phase_N));
max_phase = min(max(phase_G), max(phase_N));
common_phase = linspace(min_phase, max_phase, 2000);

% Step B: Interpolate both curves to share this common X-axis
% (unique is used to avoid errors if there are duplicate phase points)
[u_phase_G, idx_G] = unique(phase_G);
mag_G_interp = interp1(u_phase_G, mag_G_dB(idx_G), common_phase);
w_G_interp = interp1(u_phase_G, w_G_out(idx_G), common_phase);

[u_phase_N, idx_N] = unique(phase_N);
mag_N_interp = interp1(u_phase_N, mag_N_dB(idx_N), common_phase);

% Step C: Subtract the curves. An intersection occurs where the difference is zero.
mag_diff = mag_G_interp - mag_N_interp;

% Find zero-crossings by looking for where adjacent points have different signs
crossings = find(mag_diff(1:end-1) .* mag_diff(2:end) <= 0);

% Step D: Print the results and plot a marker
fprintf('\n--- Limit Cycle Intersections ---\n');
if isempty(crossings)
    fprintf('No intersections found.\n');
else
    for i = 1:length(crossings)
        idx = crossings(i);
        
        % Extract values at the crossing point
        int_w = w_G_interp(idx);
        int_phase = common_phase(idx);
        int_mag = mag_G_interp(idx);
        
        fprintf('Intersection %d:\n', i);
        fprintf('  Frequency (w) : %.4f rad/s\n', int_w);
        fprintf('  Phase         : %.2f deg\n', int_phase);
        fprintf('  Magnitude     : %.2f dB\n\n', int_mag);
        
        % Plot a red circle on the graph to show the intersection
        plot(int_phase, int_mag, 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'HandleVisibility', 'off');
    end
end
clear; clc; close all; clear functions
%% 1. Linear system G(s)


% CTRL+R / CTRL+T to comment and uncomment to switch between two incidents:

% 1)

% -------------------------------------------------------------------------

% THE PARAMETERS OF THE X-15 SOFT GLIDE LANDING PIO INCIDENT in 1959

Kp = 1;  % The pilot gain representing the gradient of the sidestick of the pilot model for the X-15 Soft Glide Landing PIO Incident in 1959
R  = 15; % Slew Rate limit of the X-15 Soft Glide Landing PIO Incident in 1959, deg/s

% The parameters of the transfer function of the aircraft for the X-15 Soft Glide Landing PIO Incident in 1959 
qco = 1;                  % pilot command amplitude of the X-15 Soft Glide Landing PIO Incident in 1959
Kp_Yp   = 13.68;            % controller gain (here pilot gain has been inherited as an internal parameter of the X-15 Soft Glide Landing PIO Incident in 1959)
M_del_e = 0.537;   % Input signal entering the rate limiter
omega_n = 2.3;     % Natural frequency
zeta_sp = 1.42 / omega_n / 2;   % Damping ratio
num = Kp_Yp*M_del_e*[1 0.82];
den = [1 1.42 omega_n^2 0];

Gac = qco*Kp*tf(num,den); % the transfer function of the aircraft for the X-15 Soft Glide Landing PIO Incident in 1959

% -------------------------------------------------------------------------

% 2)

% -------------------------------------------------------------------------

% THE PARAMETERS OF THE X-15 Flight 3-65-97, 1967 FATAL CRASH PIO INCIDENT

% Kp = 13.68;                 % the pilot gain representing the gradient of the sidestick of the pilot model for X-15 Flight 3-65-97, 1967 Fatal Crash PIO Incident's Transfer function of the Aircraft
% R   = 60; % actuator rate limit of X-15 Flight 3-65-97, 1967 Fatal Crash PIO Incident's Transfer function of the Aircraft

% qco = 1; of X-15 Flight 3-65-97, 1967 Fatal Crash PIO Incident's Transfer function of the Aircraft
% X-15 Flight 3-65-97, 1967 Fatal Crash PIO Incident's Transfer function of the Controller
% num_Gc = 5.21 * conv([1 -57.36], conv([1 4.26], [1 0.55])); 
% den_Gc = conv([1 2*0.442*22.85 22.85^2], conv([1 0], [1 1.16]));
% Gc  = tf(num_Gc, den_Gc);
% X-15 Flight 3-65-97, 1967 Fatal Crash PIO Incident's Transfer function of the Aircraft
% num_Gac = -10.524 * conv([1 1.562], conv([1 0.038], [1 0]));
% den_Gac = conv([1 2*0.212*0.088 0.088^2], conv([1 3.75], [1 -1.44]));

% Gac = qco*Kptf(num_Gac, den_Gac);

% -------------------------------------------------------------------------

% The Linearized Transfer Function
Gs = qco*Kp*Gac; % Linearized Transfer Function of the convolution of the controller and the aircraft


%% 2. Resolution and the Frequency Range
n = 1000; % Resolution
w = logspace(-2, 2, n); % Common Frequency Rangerad/s 

%% 3. Rate Limiter Element parameters (For the Category II PIO Detection with the OLOP Method)

Ai = qco*Kp;        % Input amplitude, in degs for using the X-15 Flight 3-65-97, 1967 Fatal Crash PIO Incident's Transfer function 
w_onset = R / Ai;

%% 4. Describing function N(Ai,w)
N = zeros(size(w));
for k = 1:length(w)
    alpha = w(k) / w_onset;
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
        varpi = w_onset / w(k);
        M = (4/pi)*varpi;
        phi = -acos((pi/2)*varpi);
    end
    N(k) = M*exp(1j*phi);
end

%% 5. Calculate -1/N(Ai,w)
minus_inv_N = -1 ./ N;

%% 6. Convert -1/N to an FRD model
response = reshape(minus_inv_N,1,1,[]);
sys_N = frd(response, w);

%% 7. Nichols Chart
figure(Name='Nichols Chart',NumberTitle='off');
p1 = nicholsplot(Gs, w);
hold on
p2 = nicholsplot(sys_N);

%% 8. Shift -1/N from +180 deg to -180 deg
p2.PhaseMatchingEnabled = 'on';
p2.PhaseMatchingFrequency = w(1);
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
[mag_G, phase_G, w_G_out] = nichols(Gs, w);
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
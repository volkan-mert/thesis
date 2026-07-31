%% PIO Limit Cycle Analysis using the Describing Function Method
% Simple version (uses nicholsplot)
% Plot the aircraft G(jw) and -1/N on a Nichols chart, find where they
% cross, and print the limit cycle (PIO) frequency in the command window.

clear; clc; close all;

%% 1. Aircraft transfer function
% Control law
num_claw = [5.21, -273.7855, -1425.456, -700.224];
den_claw = [1, 21.36, 545.6, 605.7, 0];
G_claw = tf(num_claw, den_claw);

% Lateral dynamics
num_dyn = [-10.524, -16.8384, -0.6209, 0];
den_dyn = [1, 2.35, -5.31, 0.184, -0.041];
G_dyn = tf(num_dyn, den_dyn);

% Full open-loop aircraft (q / q_command)
G = G_claw * G_dyn;

%% 2. Rate limiter describing function
R = 60;              % max actuator rate (deg/s)
A = 13.68;           % input amplitude (deg)
w_onset = R / A;     % onset frequency (rad/s)

w = logspace(-1, 2, 500);   % frequency vector (rad/s)

mag = zeros(size(w));   % gain of N
phi = zeros(size(w));   % phase of N (radians)

for k = 1:length(w)
    alpha = w(k) / w_onset;

    if alpha < 1
        % Region I: no saturation
        mag(k) = 1;
        phi(k) = 0;

    elseif alpha < 1.862
        % Region II: transition (cubic fit)
        mag(k) = 0.2908*alpha^3 - 1.4396*alpha^2 + 1.9232*alpha + 0.223;
        phi(k) = 0.5280*alpha^3 - 2.6213*alpha^2 + 3.5056*alpha - 1.4171;

    else
        % Region III: fully developed saturation
        wbar = 1 / alpha;
        mag(k) = (4 * wbar) / pi;
        phi(k) = -acos(pi * wbar / 2);
    end
end

% Build -1/N as a frequency response (FRD) model for nicholsplot
N = mag .* exp(1i*phi);
indf = frd(-1 ./ N, w);
indf.Name = '-1/N';

%% 3. Nichols chart with nicholsplot
figure;
opts = nicholsoptions;
opts.Grid = 'on';
opts.PhaseMatching = 'on';        % put both curves on the same branch
opts.PhaseMatchingFreq = 15;      % align near the onset frequency
opts.PhaseMatchingValue = 540;    % on the 540 deg branch

nicholsplot(G, indf, opts);
title('PIO Limit Cycle - Describing Function Method');
hold on;
axis([360 720 -40 50]);

%% 4. Get numeric data to find the crossing
[m_ac, p_ac] = nichols(G, w);
mag_ac   = 20*log10(squeeze(m_ac));   % dB
phase_ac = squeeze(p_ac);             % deg

mag_df   = 20*log10(1 ./ mag);        % dB
phase_df = 540 - phi*180/pi;          % deg (540 branch)

%% 5. Find where the two curves cross (the limit cycle point)
w_pio = NaN;

for i = 1:length(w)-1
    for j = 1:length(w)-1
        [xi, yi, t, hit] = seg_cross(phase_ac(i), mag_ac(i), ...
                                     phase_ac(i+1), mag_ac(i+1), ...
                                     phase_df(j), mag_df(j), ...
                                     phase_df(j+1), mag_df(j+1));
        if hit
            % Frequency at the crossing (interpolated along aircraft curve)
            w_pio = w(i) + t*(w(i+1) - w(i));

            % Mark it on the chart
            plot(xi, yi, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
        end
    end
end

%% 6. Print the intersection frequency in the command window
if isnan(w_pio)
    fprintf('No crossing found - no limit cycle predicted.\n');
else
    fprintf('--- Predicted Limit Cycle (PIO) ---\n');
    fprintf('Intersection frequency = %.4f rad/s (%.4f Hz)\n', w_pio, w_pio/(2*pi));
end


%% Helper function: do two line segments cross?
function [xi, yi, t, hit] = seg_cross(x1,y1, x2,y2, x3,y3, x4,y4)
    % Segment 1: (x1,y1)->(x2,y2)   Segment 2: (x3,y3)->(x4,y4)
    hit = false; 
    xi = 0; 
    yi = 0; 
    t = 0;

    d = (x2-x1)*(y4-y3) - (y2-y1)*(x4-x3);
    if abs(d) < 1e-12
        return;   % lines are parallel
    end

    t = ((x3-x1)*(y4-y3) - (y3-y1)*(x4-x3)) / d;   % along segment 1
    u = ((x3-x1)*(y2-y1) - (y3-y1)*(x2-x1)) / d;   % along segment 2

    if t >= 0 && t <= 1 && u >= 0 && u <= 1
        hit = true;
        xi = x1 + t*(x2-x1);
        yi = y1 + t*(y2-y1);
    end
end

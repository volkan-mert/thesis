clear; clc; close all

%% Define Linear System Components G(s)
s = tf('s');

% Control Law Transfer Function
num_claw = 5.21 * [1, -52.55, -273.6, -134.4];
den_claw = [1, 21.36, 545.6, 605.7, 0];
G_claw = tf(num_claw, den_claw);

% Longitudinal Dynamics Transfer Function
num_ac = -10.524 * [1, 1.6, 0.059, 0];
den_ac = [1, 2.35, -5.31, 0.184, -0.041];
G_ac = tf(num_ac, den_ac);

Kp = 1; % Pilot gain, with the pilot model taken inside the loop

% Combined open-loop transfer function G(s)
G = Kp * G_claw * G_ac;

%% Calculate -1/N Describing-Function Loci
n = 1000;
w = logspace(-1, 2, n);
R = 15; % Rate limit, deg/s
% u_rle = 0.3;% Alternative values for R = 15 deg/s
u_rle = [0.2, 0.3, 0.31, 0.5, 1, 2, 5, 10, 50]; % Alternative values for R = 15 deg/s
% u_rle = [1, 1.25, 1.5, 2, 5, 10, 50]; % Alternative values for R = 60 deg/s

sys_inv_list = cell(1, numel(u_rle));
inv_N_list = cell(1, numel(u_rle));

for k = 1:numel(u_rle)
    A_i = u_rle(k); % Input amplitude of the rate-limiting element

    % Describing function N for the fully developed region
    varpi = R ./ (A_i * w);
    N = (4/pi) * varpi .* exp(-1i * acos((pi/2) * varpi));
    N((A_i * w) <= R) = 1; % Linear region, as used in the original script

    % Inverse describing function (-1/N)
    inv_N = -1 ./ N;
    inv_N_list{k} = inv_N;

    % Reshape to 3-D for the FRD object (SISO format)
    inv_N_3D = reshape(inv_N, 1, 1, []);
    sys_inv_list{k} = frd(inv_N_3D, w);
end

%% Plot G(jw) and -1/N on the Nyquist Plot
figure('Color', 'w');

% Plot G(s) together with all -1/N curves
nyquist(G, sys_inv_list{:}, w);
grid;

legend_labels = [{'G(j\omega)'}, arrayfun(@(u) sprintf('-1/N (u_{rle} = %.4f)', u), u_rle, 'UniformOutput', false)];
legend(legend_labels, 'Location', 'northeast');

title_slew_rate = sprintf('Slew Rate, R = %d deg/s', R);
title({'The Nyquist Plot of a Rate Limiting Element', 'Open-Loop Plant G(j\omega) vs. -1/N Loci', ['\color{red} (' title_slew_rate ')']});

%% Find Geometrical Intersections in Nichols Coordinates
% Important distinction:
%   w_G is the frequency parameter of G(jw) at the graphical crossing.
%   w_N is the frequency parameter of the plotted -1/N locus at that crossing.
% For a frequency-dependent describing function, a physical limit-cycle
% solution requires the same frequency on both sides.

target_u = 0.31; % Curve to annotate, as in the screenshot

Gjw = squeeze(freqresp(G, w));
Gjw = Gjw(:).';
phase_G = rad2deg(unwrap(angle(Gjw)));
gain_G_dB = 20 * log10(abs(Gjw));

result_rows = [];
annotated = false;

for k = 1:numel(u_rle)
    inv_N = inv_N_list{k};
    phase_N = rad2deg(unwrap(angle(inv_N)));
    gain_N_dB = 20 * log10(abs(inv_N));

    [phase_int, gain_int, seg_G, seg_N, tau_G, tau_N] = ...
        polylineIntersections(phase_G, gain_G_dB, phase_N, gain_N_dB);

    for m = 1:numel(phase_int)
        % Logarithmic interpolation is appropriate because w is log-spaced.
        w_G = logInterpolate(w(seg_G(m)), w(seg_G(m)+1), tau_G(m));
        w_N = logInterpolate(w(seg_N(m)), w(seg_N(m)+1), tau_N(m));

        % Since the RLE locus depends on A*w, this is the amplitude that
        % makes the graphical intersection frequency-consistent.
        u_required_same_w = u_rle(k) * w_N / w_G;

        result_rows = [result_rows; ...
            u_rle(k), phase_int(m), gain_int(m), ...
            w_G, w_N, u_required_same_w]; %#ok<AGROW>

        % Add one screenshot-style annotation for the selected curve.
        if ~annotated && abs(u_rle(k) - target_u) < 1e-12
            ax = gca;
            hold(ax, 'on');
            plot(ax, phase_int(m), gain_int(m), 'ro', 'MarkerSize', 13, 'LineWidth', 1.8);
            text(ax, phase_int(m) + 7, gain_int(m) - 2, sprintf(['Intersection point\n', '\omega_G = %.4f rad/s\n', '\omega_N = %.4f rad/s'], w_G, w_N), 'Color', 'r', 'FontWeight', 'bold');

            annotated = true;
        end
    end
end

% If the selected -1/N curve does not intersect G(jw), display a clear
% no-limit-cycle message directly on the Nichols plot.
if ~annotated
    ax = gca;
    hold(ax, 'on');
    text(ax, 0.5, 0.5, 'NO LIMIT CYCLE!', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Color', 'r', ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'Margin', 6);
end

if isempty(result_rows)
    warning('No graphical intersection was found over the selected frequency range.');
else
    intersection_results = array2table(result_rows, 'VariableNames', {'u_rle', 'Phase_deg', 'Gain_dB', 'w_G_rad_s', 'w_N_rad_s', 'u_required_for_same_w'});

    fprintf('\nGraphical intersections in Nichols coordinates:\n');
    disp(intersection_results);

    target_rows = intersection_results(abs(intersection_results.u_rle-target_u) < 1e-12, :);
    if ~isempty(target_rows)
        fprintf(['For the plotted u_rle = %.4f curve, the visual crossing occurs at\n', ...
                 'phase = %.4f deg and gain = %.4f dB.\n', ...
                 'The plant frequency is w_G = %.6f rad/s, whereas the point on\n', ...
                 'the plotted -1/N curve corresponds to w_N = %.6f rad/s.\n', ...
                 'Therefore, this is not a same-frequency solution for u_rle = %.4f.\n', ...
                 'The frequency-consistent amplitude is approximately %.6f deg,\n', ...
                 'with limit-cycle frequency w = %.6f rad/s.\n\n'], ...
                 target_u, target_rows.Phase_deg(1), target_rows.Gain_dB(1), ...
                 target_rows.w_G_rad_s(1), target_rows.w_N_rad_s(1), ...
                 target_u, target_rows.u_required_for_same_w(1), ...
                 target_rows.w_G_rad_s(1));
    end
end

%% Local Functions
function wq = logInterpolate(w1, w2, tau)
%LOGINTERPOLATE Interpolate between two points on a logarithmic frequency grid.
    wq = exp(log(w1) + tau * (log(w2) - log(w1)));
end

function [xI, yI, seg1, seg2, t1, t2] = ...
    polylineIntersections(x1, y1, x2, y2)
%POLYLINEINTERSECTIONS Find intersections between two 2-D polylines.
% This implementation does not require Mapping Toolbox or polyxpoly.

    x1 = x1(:); y1 = y1(:);
    x2 = x2(:); y2 = y2(:);

    xI = []; yI = [];
    seg1 = []; seg2 = [];
    t1 = []; t2 = [];

    parallel_tol = 1e-12;
    bound_tol = 1e-10;
    duplicate_tol = 1e-6;

    for i = 1:(numel(x1)-1)
        p = [x1(i); y1(i)];
        r = [x1(i+1)-x1(i); y1(i+1)-y1(i)];

        if norm(r) < parallel_tol
            continue
        end

        p_min = min(p, p+r) - bound_tol;
        p_max = max(p, p+r) + bound_tol;

        for j = 1:(numel(x2)-1)
            q = [x2(j); y2(j)];
            s = [x2(j+1)-x2(j); y2(j+1)-y2(j)];

            if norm(s) < parallel_tol
                continue
            end

            q_min = min(q, q+s) - bound_tol;
            q_max = max(q, q+s) + bound_tol;

            % Fast bounding-box rejection.
            if any(p_max < q_min) || any(q_max < p_min)
                continue
            end

            denominator = cross2D(r, s);
            if abs(denominator) < parallel_tol
                continue
            end

            q_minus_p = q - p;
            tau1 = cross2D(q_minus_p, s) / denominator;
            tau2 = cross2D(q_minus_p, r) / denominator;

            if tau1 >= -bound_tol && tau1 <= 1+bound_tol && ...
               tau2 >= -bound_tol && tau2 <= 1+bound_tol

                tau1 = min(max(tau1, 0), 1);
                tau2 = min(max(tau2, 0), 1);
                point = p + tau1*r;

                % Suppress duplicate detections at adjacent segments.
                if isempty(xI) || all(hypot(xI-point(1), yI-point(2)) > duplicate_tol)
                    xI(end+1,1) = point(1); %#ok<AGROW>
                    yI(end+1,1) = point(2); %#ok<AGROW>
                    seg1(end+1,1) = i; %#ok<AGROW>
                    seg2(end+1,1) = j; %#ok<AGROW>
                    t1(end+1,1) = tau1; %#ok<AGROW>
                    t2(end+1,1) = tau2; %#ok<AGROW>
                end
            end
        end
    end
end

function value = cross2D(a, b)
%CROSS2D Scalar 2-D cross product.
    value = a(1)*b(2) - a(2)*b(1);
end

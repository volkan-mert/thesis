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
% u_rle = 0.3; % Testing a no Limit Cycle Value for R = 15 deg/s
u_rle = [0.2, 0.3, 0.31, 0.5, 1, 5]; % Alternative values for R = 15 deg/s
% u_rle = [1, 1.25, 1.5, 2]; % Alternative values for R = 60 deg/s

target_u = 0.31; % Selected -1/N curve to annotate

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

opt = nyquistoptions;
opt.Grid = 'on';
opt.ShowFullContour = 'off'; % Display the positive-frequency branch only

% Plot the open-loop system together with all inverse-DF loci
h = nyquistplot(G, sys_inv_list{:}, w);
setoptions(h, opt);

drawnow;
ax = gca;
hold(ax, 'on');
axis(ax, 'equal');

legend_labels = [{'G(j\omega)'}, ...
    arrayfun(@(u) sprintf('-1/N (u_{rle} = %.4f)', u), ...
    u_rle, 'UniformOutput', false)];
legend(ax, legend_labels, 'Location', 'best');

xlabel(ax, 'Real Axis');
ylabel(ax, 'Imaginary Axis');

title_slew_rate = sprintf('Slew Rate, R = %d deg/s', R);
title(ax, {'The Nyquist Plot of a Rate Limiting Element', ...
           'Open-Loop Plant G(j\omega) vs. -1/N Loci', ...
           ['\color{red}(' title_slew_rate ')']});

%% Find Geometrical Intersections in the Complex Plane
% Important distinction:
%   w_G is the frequency parameter of G(jw) at the graphical crossing.
%   w_N is the frequency parameter of the plotted -1/N locus at that crossing.
% For a frequency-dependent describing function, a physical limit-cycle
% solution requires the same frequency on both sides.

Gjw = squeeze(freqresp(G, w));
Gjw = Gjw(:).';

real_G = real(Gjw);
imag_G = imag(Gjw);

result_rows = [];
annotated = false;

for k = 1:numel(u_rle)
    inv_N = inv_N_list{k};
    real_N = real(inv_N);
    imag_N = imag(inv_N);

    [real_int, imag_int, seg_G, seg_N, tau_G, tau_N] = ...
        polylineIntersections(real_G, imag_G, real_N, imag_N);

    for m = 1:numel(real_int)
        % Logarithmic interpolation is appropriate because w is log-spaced.
        w_G = logInterpolate(w(seg_G(m)), w(seg_G(m)+1), tau_G(m));
        w_N = logInterpolate(w(seg_N(m)), w(seg_N(m)+1), tau_N(m));

        % Since the RLE locus depends on A*w, this is the amplitude that
        % makes the graphical intersection frequency-consistent.
        u_required_same_w = u_rle(k) * w_N / w_G;

        result_rows = [result_rows; ...
            u_rle(k), real_int(m), imag_int(m), ...
            w_G, w_N, u_required_same_w]; %#ok<AGROW>

        % Annotate and automatically zoom around the first intersection
        % of the selected -1/N curve.
        if ~annotated && abs(u_rle(k) - target_u) < 1e-12
            x_int = real_int(m);
            y_int = imag_int(m);

            plot(ax, x_int, y_int, 'ro', ...
                'MarkerSize', 13, 'LineWidth', 1.8);

            %% Automatically arrange xlim and ylim around the intersection
            % Use neighbouring samples from both curves so that the local
            % crossing geometry remains visible rather than zooming only
            % to the single intersection point.
            zoom_samples = 30;

            idx_G = max(1, seg_G(m)-zoom_samples): ...
                    min(numel(w), seg_G(m)+1+zoom_samples);
            idx_N = max(1, seg_N(m)-zoom_samples): ...
                    min(numel(w), seg_N(m)+1+zoom_samples);

            x_local = [real_G(idx_G), real_N(idx_N)];
            y_local = [imag_G(idx_G), imag_N(idx_N)];

            finite_local = isfinite(x_local) & isfinite(y_local);
            x_local = x_local(finite_local);
            y_local = y_local(finite_local);

            % Maximum local distance from the intersection. A minimum
            % radius prevents an excessively tight view for dense grids.
            local_radius = max(hypot(x_local-x_int, y_local-y_int));
            minimum_radius = 0.05 * max([abs(x_int), abs(y_int), 1]);
            view_radius = 1.35 * max(local_radius, minimum_radius);

            % Equal x- and y-ranges preserve the Nyquist-plane geometry.
            xlim(ax, x_int + [-view_radius, view_radius]);
            ylim(ax, y_int + [-view_radius, view_radius]);
            axis(ax, 'equal');

            % Determine the label offset after applying the new limits.
            xl = xlim(ax);
            yl = ylim(ax);
            x_offset = 0.04 * diff(xl);
            y_offset = 0.04 * diff(yl);

            % Write the limit-cycle frequency next to the intersection.
            text(ax, x_int + x_offset, y_int + y_offset, ...
                sprintf('\\omega_{LC} = %.4f rad/s', w_G), ...
                'Color', 'r', ...
                'FontSize', 11, ...
                'FontWeight', 'bold', ...
                'Interpreter', 'tex', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom', ...
                'BackgroundColor', 'w', ...
                'Margin', 3, ...
                'Clipping', 'on');

            annotated = true;
        end
    end
end

% If the selected -1/N curve does not intersect G(jw), display a clear
% no-limit-cycle message directly on the Nyquist plot.
if ~annotated
    text(ax, 0.50, 0.50, 'NO LIMIT CYCLE!', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Color', 'r', ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'Margin', 6);
end

%% Display Numerical Intersection Results
if isempty(result_rows)
    warning('No graphical intersection was found over the selected frequency range.');
else
    intersection_results = array2table(result_rows, ...
        'VariableNames', {'u_rle', 'RealPart', 'ImagPart', ...
                          'w_G_rad_s', 'w_N_rad_s', ...
                          'u_required_for_same_w'});

    fprintf('\nGraphical intersections in Nyquist coordinates:\n');
    disp(intersection_results);

    target_rows = intersection_results( ...
        abs(intersection_results.u_rle-target_u) < 1e-12, :);

    if ~isempty(target_rows)
        fprintf(['For the plotted u_rle = %.4f curve, the visual crossing occurs at\n', ...
                 'Re = %.6f and Im = %.6f.\n', ...
                 'The plant frequency is w_G = %.6f rad/s, whereas the point on\n', ...
                 'the plotted -1/N curve corresponds to w_N = %.6f rad/s.\n', ...
                 'Therefore, this is not necessarily a same-frequency solution.\n', ...
                 'The frequency-consistent amplitude is approximately %.6f deg,\n', ...
                 'with limit-cycle frequency w = %.6f rad/s.\n\n'], ...
                 target_u, target_rows.RealPart(1), target_rows.ImagPart(1), ...
                 target_rows.w_G_rad_s(1), target_rows.w_N_rad_s(1), ...
                 target_rows.u_required_for_same_w(1), ...
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
    duplicate_tol = 1e-8;

    for i = 1:(numel(x1)-1)
        p = [x1(i); y1(i)];
        r = [x1(i+1)-x1(i); y1(i+1)-y1(i)];

        if any(~isfinite([p; r])) || norm(r) < parallel_tol
            continue
        end

        p_min = min(p, p+r) - bound_tol;
        p_max = max(p, p+r) + bound_tol;

        for j = 1:(numel(x2)-1)
            q = [x2(j); y2(j)];
            s = [x2(j+1)-x2(j); y2(j+1)-y2(j)];

            if any(~isfinite([q; s])) || norm(s) < parallel_tol
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
                if isempty(xI) || ...
                   all(hypot(xI-point(1), yI-point(2)) > duplicate_tol)
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

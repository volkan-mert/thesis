%% APPLIED NONLINEAR CONTROL: DESCRIBING FUNCTION ANALYSIS
% Frequency-Dependent Describing Function (FDDF) Limit Cycle Detection
% Plant Model: Gs_ac = Gs_ldynac

clear; clc; close all;

%% 1. Linear Aircraft Plant Definition (Gs_ldynac only)
num_ldynac = [-10.524, -16.8384, -0.6209, 0];
den_ldynac = [1, 2.35, -5.31, 0.184, -0.041];
Gs_ldynac  = tf(num_ldynac, den_ldynac);

Gs_ac      = Gs_ldynac;               % Linear OLTF set to Gs_ldynac

%% 2. Rate Limiter Element (RLE) Describing Function Evaluation
R = 60;                       % Slew rate / Max actuator rate (deg/s)
A = 13.68;                    % Input amplitude (deg)
w = logspace(-1, 2, 2000);    % Frequency vector (rad/s)

% Calculate onset and critical saturation frequencies
w_onset = R / A;              
w_crit  = 1.862 * w_onset;    

mag = zeros(1, length(w));
phi = zeros(1, length(w));

for k = 1:length(w)
    alpha = w(k) / w_onset;   
    if alpha < 1
        % Region I: No saturation (Eq. 15)
        mag(k) = 1; phi(k) = 0;
    elseif alpha < 1.862
        % Region II: Transition via cubic spline interpolation (Eq. 16)
        mag(k) = 0.2908*alpha^3 - 1.4396*alpha^2 + 1.9232*alpha + 0.223;
        phi(k) = 0.5280*alpha^3 - 2.6213*alpha^2 + 3.5056*alpha - 1.4171; % rad
    else
        % Region III: Fully developed saturation (Eq. 17)
        wbar = 1 / alpha;    
        mag(k) = (4 * wbar) / pi;
        phi(k) = -acos(pi * wbar / 2);                                  % rad
    end
end

%% 3. FRD Model Creation & Plotting
N_jw = mag .* exp(1j * phi);
sys_rle = frd(N_jw, w);
indf_rle = -1/sys_rle; 

figure;
nyquistplot(Gs_ac, indf_rle);
title('FREQUENCY-DEPENDENT DESCRIBING FUNCTION (Gs\_ldynac)', 'FontWeight', 'bold');
grid on;

%% 4. Mark the critical point (-1, 0)
hold on;
h_crit = plot(-1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
text(-1.05, 0.2, ' (-1, 0i)', 'Color', 'r', 'FontWeight', 'bold');
ylim([-5, 5]);
xlim([-5, 5]);

%% 5. Find Intersection Points

% Extract frequency response vector for Gs_ac
G_jw = squeeze(freqresp(Gs_ac, w)).'; % 1 x N vector
invN_jw = -1 ./ N_jw;                 % 1 x N vector

% --- Inline FDDF Evaluators (Prevents function scoping errors) ---
get_G = @(wv) squeeze(freqresp(Gs_ac, wv));
get_invN = @(wv) -1 / (...
    ( (wv/w_onset < 1) * 1 + ...
      (wv/w_onset >= 1 & wv/w_onset < 1.862) * (0.2908*(wv/w_onset)^3 - 1.4396*(wv/w_onset)^2 + 1.9232*(wv/w_onset) + 0.223) + ...
      (wv/w_onset >= 1.862) * (4 / (pi * (wv/w_onset))) ) * ...
    exp(1j * ( (wv/w_onset >= 1 & wv/w_onset < 1.862) * (0.5280*(wv/w_onset)^3 - 2.6213*(wv/w_onset)^2 + 3.5056*(wv/w_onset) - 1.4171) + ...
               (wv/w_onset >= 1.862) * (-acos(pi / (2 * (wv/w_onset)))) )) ...
);

%% A. Frequency-Matched Search (True FDDF Limit Cycle)
dist_vec = abs(G_jw - invN_jw);
[min_dist, coarse_idx] = min(dist_vec);

% Objective function using inline handles
obj_fun = @(w_val) abs(get_G(w_val) - get_invN(w_val));
w_search_min = w(max(1, coarse_idx - 10));
w_search_max = w(min(length(w), coarse_idx + 10));

[w_lc, fval] = fminbnd(obj_fun, w_search_min, w_search_max);
G_lc = get_G(w_lc);

fprintf('\n=======================================================\n');
fprintf('  1. FREQUENCY-MATCHED INTERSECTION (LIMIT CYCLE)\n');
fprintf('=======================================================\n');
fprintf('  Frequency (w_lc)   : %.4f rad/s\n', w_lc);
fprintf('  Complex Coordinate : %.4f + %.4fi\n', real(G_lc), imag(G_lc));
fprintf('  Magnitude Error    : %.4e\n', fval);

%% B. Geometric 2D Curve Intersection (Spatial Crossing)
x1 = real(G_jw);    y1 = imag(G_jw);
x2 = real(invN_jw); y2 = imag(invN_jw);

[x_geom, y_geom, idx1, idx2] = find_2d_intersections(x1, y1, x2, y2);

fprintf('\n=======================================================\n');
fprintf('  2. GEOMETRIC 2D CROSSING POINTS (%d Found)\n', length(x_geom));
fprintf('=======================================================\n');
for i = 1:length(x_geom)
    fprintf('  Point %d: (%.4f, %.4fi) | w_G = %.4f rad/s, w_N = %.4f rad/s\n', ...
        i, x_geom(i), y_geom(i), w(idx1(i)), w(idx2(i)));
end

%% 6. Plot Intersections on Nyquist Diagram
figure(1); hold on;

% Proxy handles for legend binding (prevents controllib chart warnings)
h_G    = plot(NaN, NaN, 'b-', 'LineWidth', 1.5);
h_invN = plot(NaN, NaN, 'r--', 'LineWidth', 1.5);

leg_handles = [h_G, h_invN, h_crit];
leg_labels  = {'G(s)', '-1/N(j\omega)', 'Critical Pt (-1,0)'};

% Plot Frequency-Matched (Limit Cycle) Point
if fval < 0.05
    h_lc = plot(real(G_lc), imag(G_lc), 'ko', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'g');
    text(real(G_lc) + 0.2, imag(G_lc) + 0.3, ...
        sprintf('  Limit Cycle\n  \\omega = %.2f rad/s', w_lc), ...
        'Color', 'k', 'FontWeight', 'bold', 'FontSize', 9);
    leg_handles(end+1) = h_lc;
    leg_labels{end+1}  = 'True LC Point';
end

% Plot Geometric Crossings
if ~isempty(x_geom)
    h_geom = plot(x_geom(1), y_geom(1), 'mx', 'MarkerSize', 12, 'LineWidth', 2);
    for i = 2:length(x_geom)
        plot(x_geom(i), y_geom(i), 'mx', 'MarkerSize', 12, 'LineWidth', 2);
    end
    for i = 1:length(x_geom)
        text(x_geom(i) + 0.2, y_geom(i) - 0.3, ...
            sprintf('  Geom X-ing %d', i), ...
            'Color', 'm', 'FontWeight', 'bold', 'FontSize', 9);
    end
    leg_handles(end+1) = h_geom;
    leg_labels{end+1}  = 'Geometric Crossing';
end

legend(leg_handles, leg_labels, 'Location', 'northeast');


%% LOCAL HELPER FUNCTION (Placed at the very end of the file)


function [x_int, y_int, idx1, idx2] = find_2d_intersections(x1, y1, x2, y2)
    x_int = []; y_int = []; idx1 = []; idx2 = [];
    N1 = length(x1) - 1;
    N2 = length(x2) - 1;
    
    for i = 1:N1
        p1 = [x1(i); y1(i)]; p2 = [x1(i+1); y1(i+1)];
        dp = p2 - p1;
        
        min_x1 = min(p1(1), p2(1)); max_x1 = max(p1(1), p2(1));
        min_y1 = min(p1(2), p2(2)); max_y1 = max(p1(2), p2(2));
        
        for j = 1:N2
            q1 = [x2(j); y2(j)]; q2 = [x2(j+1); y2(j+1)];
            
            if max(q1(1), q2(1)) < min_x1 || min(q1(1), q2(1)) > max_x1 || ...
               max(q1(2), q2(2)) < min_y1 || min(q1(2), q2(2)) > max_y1
                continue;
            end
            
            dq = q2 - q1;
            A_mat = [dp, -dq];
            if abs(det(A_mat)) < 1e-12, continue; end
            
            tu = A_mat \ (q1 - p1);
            t = tu(1); u = tu(2);
            
            if t >= 0 && t <= 1 && u >= 0 && u <= 1
                pt = p1 + t * dp;
                x_int(end+1) = pt(1); 
                y_int(end+1) = pt(2); 
                idx1(end+1) = i;      
                idx2(end+1) = j;      
            end
        end
    end
end
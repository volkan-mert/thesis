clear; clc; close all

%% 1. Define the open-loop system
s = tf('s');

% Control law
num_claw = 5.21 * [1, -52.55, -273.6, -134.4];
den_claw = [1, 21.36, 545.6, 605.7, 0];
G_claw = tf(num_claw, den_claw);

% Aircraft longitudinal dynamics
num_ac = -10.524 * [1, 1.6, 0.059, 0];
den_ac = [1, 2.35, -5.31, 0.184, -0.041];
G_ac = tf(num_ac, den_ac);

Kp = 1;
G = Kp * G_claw * G_ac;

%% 2. Calculate the inverse describing-function curves
w = logspace(-1, 2, 1000);       % Frequency vector, rad/s
R = 15;                           % Rate limit, deg/s
u_rle = 0.3; % Testing a no Limit Cycle Value for R = 15 deg/s
% u_rle = [0.2, 0.3, 0.31, 0.5, 1, 5]; % Alternative Values for R = 0.15
% u_rle = [1, 1.25, 1.5, 2]; % Alternative values for R = 60 deg/s
target_u = 0.31;                  % Curve used for intersection detection

inv_N_all = cell(1, length(u_rle));
frd_all = cell(1, length(u_rle));

for k = 1:length(u_rle)
    A = u_rle(k);
    varpi = R ./ (A*w);

    % Describing function
    N = (4/pi) * varpi .* exp(-1i*acos((pi/2)*varpi));
    N(A*w <= R) = 1;

    % Inverse describing function
    inv_N_all{k} = -1 ./ N;
    frd_all{k} = frd(reshape(inv_N_all{k}, 1, 1, []), w);
end

%% 3. Draw the Nyquist plot
figure('Color', 'w');

opt = nyquistoptions;
opt.Grid = 'on';
opt.ShowFullContour = 'off';

h = nyquistplot(G, frd_all{:}, w);
setoptions(h, opt);

drawnow;
ax = gca;
hold(ax, 'on');
axis(ax, 'equal');

legend_text = [{'G(j\omega)'}, ...
    arrayfun(@(A) sprintf('-1/N, u_{rle}=%.2f', A), ...
    u_rle, 'UniformOutput', false)];
legend(ax, legend_text, 'Location', 'best');

xlabel(ax, 'Real Axis');
ylabel(ax, 'Imaginary Axis');
title(ax, sprintf('Nyquist Plot for R = %.0f deg/s', R));

%% 4. Find the intersection for the selected curve
Gjw = squeeze(freqresp(G, w));
Gjw = Gjw(:).';

xG = real(Gjw);
yG = imag(Gjw);

k_target = find(abs(u_rle-target_u) < 1e-12, 1);
xN = real(inv_N_all{k_target});
yN = imag(inv_N_all{k_target});

intersection_found = false;

% Compare every segment of G(jw) with every segment of -1/N
for i = 1:length(w)-1
    p1 = [xG(i),   yG(i)];
    p2 = [xG(i+1), yG(i+1)];

    for j = 1:length(w)-1
        q1 = [xN(j),   yN(j)];
        q2 = [xN(j+1), yN(j+1)];

        [is_crossing, t] = lineIntersection(p1, p2, q1, q2);

        if is_crossing
            % Intersection coordinates
            x_int = p1(1) + t*(p2(1)-p1(1));
            y_int = p1(2) + t*(p2(2)-p1(2));

            % Frequency at the intersection
            w_LC = exp(log(w(i)) + t*(log(w(i+1))-log(w(i))));

            intersection_found = true;
            i_int = i;
            j_int = j;
            break
        end
    end

    if intersection_found
        break
    end
end

%% 5. Mark the result and arrange the plot limits
if intersection_found
    plot(ax, x_int, y_int, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName','Limit Cycle Frequency');

    % Use nearby points to define automatic x and y limits
    n_zoom = 25;
    i1 = max(1, i_int-n_zoom);
    i2 = min(length(w), i_int+n_zoom);
    j1 = max(1, j_int-n_zoom);
    j2 = min(length(w), j_int+n_zoom);

    x_near = [xG(i1:i2), xN(j1:j2)];
    y_near = [yG(i1:i2), yN(j1:j2)];

    x_range = max(x_near) - min(x_near);
    y_range = max(y_near) - min(y_near);

    % Avoid zero or very small plotting ranges
    x_range = max(x_range, 0.1);
    y_range = max(y_range, 0.1);

    xlim(ax, [min(x_near)-0.25*x_range, max(x_near)+0.25*x_range]);
    ylim(ax, [min(y_near)-0.25*y_range, max(y_near)+0.25*y_range]);

    % Write the frequency near the intersection
    xl = xlim(ax);
    yl = ylim(ax);

    text(ax, x_int + 0.03*diff(xl), y_int + 0.03*diff(yl), ...
        sprintf('\\omega_{LC} = %.4f rad/s', w_LC), ...
        'Color', 'r', ...
        'FontSize', 11, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'Margin', 2);

    fprintf('Limit-cycle frequency = %.4f rad/s\n', w_LC);
else
    text(ax, 0.5, 0.5, 'NO LIMIT CYCLE!', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'Color', 'r', ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w');

    fprintf('NO LIMIT CYCLE!\n');
end

%% Local function
function [crossing, t] = lineIntersection(p1, p2, q1, q2)
% Find the intersection of two line segments.

    r = p2 - p1;
    s = q2 - q1;

    denominator = r(1)*s(2) - r(2)*s(1);
    crossing = false;
    t = 0;

    % Parallel lines do not have one unique intersection
    if abs(denominator) < 1e-12
        return
    end

    qp = q1 - p1;
    t_temp = (qp(1)*s(2) - qp(2)*s(1)) / denominator;
    u_temp = (qp(1)*r(2) - qp(2)*r(1)) / denominator;

    % The segments intersect when both parameters are between 0 and 1
    if t_temp >= 0 && t_temp <= 1 && u_temp >= 0 && u_temp <= 1
        crossing = true;
        t = t_temp;
    end
end

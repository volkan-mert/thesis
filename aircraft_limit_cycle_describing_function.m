%% AIRCRAFT LIMIT-CYCLE ANALYSIS USING DESCRIBING FUNCTION
% Adapted from the Van der Pol describing-function script.
%
% Linear system:
% G(s) = Kf * sys_ctrl * sys_dyn
%
% The describing function is kept the same as in the attached script:
% N(A,w) = j*mu*w*A^2/4
%
% Limit-cycle condition:
% G(jw)*N(A,w) = -1

clear; clc; close all;

%% 1. Forward Gain
Kf = 13.68;

%% 2. Control Law
num_ctrl = 5.21 * [1, -52.55, -273.6, -134.4];
den_ctrl = [1, 21.36, 545.6, 605.7, 0];

sys_ctrl = tf(num_ctrl, den_ctrl);

%% 3. Longitudinal Dynamics of the Aircraft
num_dyn = -10.524 * [1, 1.6, 0.059, 0];
den_dyn = [1, 2.35, -5.31, 0.184, -0.041];

sys_dyn = tf(num_dyn, den_dyn);

%% 4. Complete Linear System
sys = Kf * sys_ctrl * sys_dyn;

%% 5. Describing Function Parameter
% Same describing function used in the attached Van der Pol script.
mu = 1;

% N(A,w) = j*mu*w*A^2/4

%% 6. Frequency Response G(jw)
% A wide logarithmic range is used because the aircraft system has
% possible crossings at both low and high frequencies.
w = logspace(-3, 3, 30000);       % rad/s
G = squeeze(freqresp(sys, w));
G = G(:).';

%% 7. Find Frequencies Where real(G(jw)) = 0
% For N(A,w) = j*mu*w*A^2/4, a limit cycle requires G(jw) to be
% on the positive imaginary axis.

realG = real(G);

% Find intervals where real(G) changes sign.
k_cross = find(realG(1:end-1).*realG(2:end) <= 0);

w_cross = [];

for k = 1:length(k_cross)
    w1 = w(k_cross(k));
    w2 = w(k_cross(k)+1);

    % More accurate zero of real(G(jw)).
    wc = fzero(@(ww) real(squeeze(freqresp(sys, ww))), [w1 w2]);

    % Avoid saving the same crossing twice.
    if isempty(w_cross) || all(abs(w_cross - wc) > 1e-5)
        w_cross = [w_cross; wc];
    end
end

%% 8. Calculate Limit-Cycle Amplitude at Each Valid Crossing
% At real(G)=0:
%
% G(jw) = j*b
% N(A,w) = j*mu*w*A^2/4
%
% Therefore
% G*N = -b*mu*w*A^2/4 = -1
%
% and
% A = sqrt(4/(mu*w*b))
%
% A real positive solution exists only when imag(G) > 0.

result = [];

for k = 1:length(w_cross)

    wc = w_cross(k);
    Gc = squeeze(freqresp(sys, wc));

    if imag(Gc) > 0

        A_limit = sqrt(4/(mu*wc*imag(Gc)));

        N_limit = 1i*mu*wc*A_limit^2/4;
        L_limit = Gc*N_limit;

        % Columns:
        % 1 = amplitude A
        % 2 = frequency w
        % 3 = real(G*N)
        % 4 = imag(G*N)
        % 5 = real(G)
        % 6 = imag(G)
        result = [result;
                  A_limit, wc, real(L_limit), imag(L_limit), ...
                  real(Gc), imag(Gc)];
    end
end

%% 9. Display Results
fprintf('\nLIMIT-CYCLE RESULTS\n');
fprintf('-------------------\n');

if isempty(result)
    fprintf('NO LIMIT CYCLE FOUND!\n\n');
else
    for k = 1:size(result,1)
        fprintf('Candidate %d:\n', k);
        fprintf('A = %.6f\n', result(k,1));
        fprintf('w = %.6f rad/s\n', result(k,2));
        fprintf('f = %.6f Hz\n', result(k,2)/(2*pi));
        fprintf('G(jw) = %.6f %+.6fj\n', result(k,5), result(k,6));
        fprintf('G(jw)N(A,w) = %.6f %+.6fj\n\n', ...
                result(k,3), result(k,4));
    end
end

%% 10. Nyquist Plot: G(jw) and -1/N(A,w)
figure(Name='G(jw) and -1/N(A,w)', NumberTitle='off');

plot(real(G), imag(G), 'LineWidth', 1.5);
hold on;

% Plot -1/N(A,w) using the detected limit-cycle amplitudes.
if ~isempty(result)
    for k = 1:size(result,1)
        A_plot = result(k,1);
        N = 1i*mu*w*A_plot^2/4;
        minus_inv_N = -1./N;

        plot(real(minus_inv_N), imag(minus_inv_N), '--', 'LineWidth', 1.2);

        % Mark the intersection on G(jw).
        Gc = result(k,5) + 1i*result(k,6);
        plot(real(Gc), imag(Gc), 'ko', 'MarkerFaceColor', 'k');
        text(real(Gc), imag(Gc), ...
             sprintf('  A = %.3f, \\omega = %.3f', result(k,1), result(k,2)));
    end
end

grid on;
xlabel('Real');
ylabel('Imaginary');
title('Describing-Function Limit-Cycle Condition');
subtitle('G(j\omega) = K_f G_{ctrl}(j\omega)G_{dyn}(j\omega)');

if isempty(result)
    legend('G(j\omega)', 'Location', 'best');
else
    legend_entries = {'G(j\omega)'};
    for k = 1:size(result,1)
        legend_entries{end+1} = sprintf('-1/N, A = %.3f', result(k,1));
        legend_entries{end+1} = sprintf('Limit cycle %d', k);
    end
    legend(legend_entries, 'Location', 'best');
end

%% 11. Plot G(jw)N(A,w) Around the Critical Point (-1,0)
figure(Name='G(jw)N(A,w)', NumberTitle='off');
hold on;

if isempty(result)
    % Use some amplitudes only for visualization if no solution is found.
    A_values = [0.1 0.2 0.3 0.4 0.5];
else
    % Plot the amplitudes that satisfy the limit-cycle condition.
    A_values = result(:,1).';
end

for k = 1:length(A_values)
    A = A_values(k);
    N = 1i*mu*w*A^2/4;
    L = G.*N;

    plot(real(L), imag(L), 'LineWidth', 1.5);
end

% Critical point.
plot(-1, 0, 'ko', 'MarkerFaceColor', 'k');
text(-1, 0, '  (-1,0)');

% Mark each calculated limit-cycle point.
if ~isempty(result)
    for k = 1:size(result,1)
        plot(result(k,3), result(k,4), 'ro', 'MarkerFaceColor', 'r');
        text(result(k,3), result(k,4), ...
             sprintf('  A = %.3f, \\omega = %.3f', result(k,1), result(k,2)));
    end
end

grid on;
xlabel('Real');
ylabel('Imaginary');
title('Nyquist Plot of G(j\omega)N(A,\omega)');
subtitle('Limit cycle occurs at (-1,0)');

legend_entries = cell(1,length(A_values)+1);
for k = 1:length(A_values)
    legend_entries{k} = sprintf('A = %.3f', A_values(k));
end
legend_entries{end} = '(-1,0)';
legend(legend_entries, 'Location', 'best');

% Zoom around the critical point when solutions exist.
if ~isempty(result)
    xlim([-2 0.5]);
    ylim([-1.5 1.5]);
end

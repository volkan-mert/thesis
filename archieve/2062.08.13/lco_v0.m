
%% LIMIT CYCLE ANALYSIS WITH RLE DESCRIBING FUNCTION

clear;
clc;
close all;

%% 1. Aircraft dynamics

numGac = [-10.5240, -16.8384, -0.6247, 0];
denGac = [1, 2.3473, -5.3061, -0.1836, -0.0418];

Gac = tf(numGac, denGac);


%% 2. Controller

numGc = [5.21, -273.7855, -1425.2, -700.1952];
denGc = [1, 21.3594, 545.5538, 605.6621, 0];

Gc = tf(numGc, denGc);


%% 3. Linear systems

Kp = 1;

% Open-loop linear part
L = Kp * Gc * Gac;

% Linear closed-loop system
sysCL = feedback(L, 1);


%% 4. Rate limiter

R  = 60;          % deg/s
Ai = 13.68;       % deg

w_onset = R / Ai;

w = logspace(log10(0.1*w_onset), ...
             log10(10*w_onset), 2000);


%% 5. RLE describing function for Ai = 13.68 deg

N = rle_df(Ai, w, R);

minus_inv_N = -1 ./ N;


%% 6. Frequency responses

G_OL = squeeze(freqresp(L, w));
G_OL = G_OL(:).';

G_CL = squeeze(freqresp(sysCL, w));
G_CL = G_CL(:).';


%% ==============================================================
%  PART 1
%  Check sysCL*N
% ==============================================================

% If N is placed around the already closed-loop sysCL,
% the condition would be:
%
%        1 + sysCL(jw)*N(A,w) = 0
%
% or
%
%        sysCL(jw) = -1/N(A,w)

error_CL = abs(1 + G_CL .* N);

[min_error_CL, kCL] = min(error_CL);

w_best_CL = w(kCL);

fprintf('\n--------------------------------------------\n');
fprintf('CHECK OF sysCL * N\n');
fprintf('--------------------------------------------\n');
fprintf('Closest frequency = %.4f rad/s\n', w_best_CL);
fprintf('Minimum error     = %.6f\n', min_error_CL);

if min_error_CL < 1e-3

    fprintf('LIMIT CYCLE FOUND for sysCL*N\n');

else

    fprintf('NO LIMIT CYCLE for sysCL*N at Ai = %.2f deg\n', Ai);

end


%% Plot sysCL against -1/N

figure;

plot(real(G_CL), imag(G_CL), ...
    'LineWidth', 1.5);

hold on;

plot(real(minus_inv_N), imag(minus_inv_N), ...
    '--', 'LineWidth', 1.5);

grid on;
axis equal;

xlabel('Real');
ylabel('Imaginary');

title('sysCL(j\omega) and -1/N(A_i,\omega)');

legend('sysCL(j\omega)', ...
       '-1/N(A_i,\omega)', ...
       'Location', 'best');


%% ==============================================================
%  PART 2
%  Draw Gc*Gac against -1/N
% ==============================================================

figure;

plot(real(G_OL), imag(G_OL), ...
    'LineWidth', 1.5);

hold on;

plot(real(minus_inv_N), imag(minus_inv_N), ...
    '--', 'LineWidth', 1.5);

grid on;
axis equal;

xlabel('Real');
ylabel('Imaginary');

title('G_c(j\omega)G_{ac}(j\omega) and -1/N(A_i,\omega)');

legend('G_cG_{ac}', ...
       '-1/N(A_i,\omega)', ...
       'Location', 'best');


%% ==============================================================
%  PART 3
%  Find actual limit-cycle amplitude and frequency
%
%  Gc(jw)*Gac(jw) = -1/N(A,w)
% ==============================================================

% Search ranges

A_search = linspace(1, 60, 300);       % deg
w_search = logspace(-1, 2, 2000);      % rad/s

G_search = squeeze(freqresp(L, w_search));
G_search = G_search(:).';


% Initial values

best_error = inf;

A0 = NaN;
w0 = NaN;


% Search over amplitude and frequency

for i = 1:length(A_search)

    A_test = A_search(i);

    N_test = rle_df(A_test, w_search, R);

    error_test = abs(1 + G_search .* N_test);

    [error_now, k] = min(error_test);

    if error_now < best_error

        best_error = error_now;

        A0 = A_test;
        w0 = w_search(k);

    end

end


%% Refine the result with fminsearch

x0 = log([w0, A0]);

options = optimset('Display', 'off', ...
                   'TolX', 1e-10, ...
                   'TolFun', 1e-12);

x = fminsearch(@(x) LC_error(x, L, R), ...
               x0, options);


% Final result

w_LC = exp(x(1));
A_LC = exp(x(2));

N_LC = rle_df(A_LC, w_LC, R);

G_LC = squeeze(freqresp(L, w_LC));

LC_error_final = abs(1 + G_LC*N_LC);

alpha_LC = w_LC*A_LC/R;


fprintf('\n--------------------------------------------\n');
fprintf('DESCRIBING FUNCTION LIMIT CYCLE\n');
fprintf('--------------------------------------------\n');

fprintf('Amplitude A     = %.4f deg\n', A_LC);
fprintf('Frequency w     = %.4f rad/s\n', w_LC);
fprintf('Frequency       = %.4f Hz\n', w_LC/(2*pi));
fprintf('alpha           = %.4f\n', alpha_LC);
fprintf('Final error     = %.8e\n', LC_error_final);


if LC_error_final < 1e-4

    fprintf('\nLIMIT CYCLE FOUND!\n');

else

    fprintf('\nNO LIMIT CYCLE FOUND!\n');

end


%% 7. Plot at the limit-cycle amplitude

w_plot = logspace(-1, 2, 3000);

G_plot = squeeze(freqresp(L, w_plot));
G_plot = G_plot(:).';

N_plot = rle_df(A_LC, w_plot, R);

minus_inv_N_plot = -1 ./ N_plot;


figure;

plot(real(G_plot), imag(G_plot), ...
    'LineWidth', 1.5);

hold on;

plot(real(minus_inv_N_plot), ...
     imag(minus_inv_N_plot), ...
     '--', 'LineWidth', 1.5);


% Limit-cycle point

plot(real(G_LC), imag(G_LC), ...
    'ko', ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 7);


text(real(G_LC), imag(G_LC), ...
    sprintf('  A = %.2f deg,  \\omega = %.3f rad/s', ...
    A_LC, w_LC));


grid on;
axis equal;

xlabel('Real');
ylabel('Imaginary');

title('Describing Function Limit Cycle');

legend('G_cG_{ac}', ...
       '-1/N(A,\omega)', ...
       'Limit Cycle', ...
       'Location', 'best');


%% ==============================================================
%  LOCAL FUNCTION: RLE DESCRIBING FUNCTION
% ==============================================================

function N = rle_df(A, w, R)

% Onset frequency
w_onset = R/A;

% Normalized frequency
alpha = w/w_onset;

% Initialize
M   = ones(size(w));
phi = zeros(size(w));


% Three regions

region1 = alpha <= 1;

region2 = alpha > 1 & alpha < 1.862;

region3 = alpha >= 1.862;


%% Region I

M(region1)   = 1;
phi(region1) = 0;


%% Region II

a = alpha(region2);

M(region2) = ...
    0.2908*a.^3 ...
    - 1.4396*a.^2 ...
    + 1.9232*a ...
    + 0.223;

phi(region2) = ...
    0.528*a.^3 ...
    - 2.6213*a.^2 ...
    + 3.5056*a ...
    - 1.4171;


%% Region III

varpi = w_onset ./ w(region3);

M(region3) = (4/pi).*varpi;

phi(region3) = ...
    -acos((pi/2).*varpi);


%% Complex describing function

N = M .* exp(1j*phi);

end


%% ==============================================================
%  LOCAL FUNCTION: LIMIT-CYCLE ERROR
% ==============================================================

function error = LC_error(x, G, R)

% Log variables are used so A and w stay positive

w = exp(x(1));
A = exp(x(2));

N = rle_df(A, w, R);

Gjw = squeeze(freqresp(G, w));

% Limit-cycle condition:
%
%       1 + G(jw)*N(A,w) = 0

error = abs(1 + Gjw*N)^2;

end
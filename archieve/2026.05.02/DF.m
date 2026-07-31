%% Rate Limiter Describing Function Implementation
clear; clc; close all;

%% 1. Define Example System Parameters
R     = 10;                     % Rate limit value (units/s)
u_rle = 5;                      % Input amplitude
omega = logspace(-1, 2, 500);   % Frequency vector (rad/s) from 0.1 to 100

%% 2. Compute Onset Frequency and Describing Function
% Equation 1: Onset Frequency
omega_onset = R / u_rle;

% Equation 2: Describing Function N_rle(j*omega, u_rle)
% Note: Using element-wise operators (./ and .*) for vector compatibility
arg_acos = (pi / 2) * (omega_onset ./ omega);

% Evaluate N_rle
N_rle = (4 / pi) * (omega_onset ./ omega) .* exp(-1j * acos(arg_acos));

%% 3. Handle Linear Region (Optional Physical Saturation Check)
% In physical rate limiters, if the frequency is below the rate-limiting 
% threshold (arg_acos > 1), acos() produces complex numbers. In reality, 
% the system operates linearly (N_rle = 1 + 0j) before rate limiting occurs.
linear_region = (arg_acos > 1);
N_rle(linear_region) = 1.0 + 0j;

%% 4. Plot Magnitude and Phase
figure('Name', 'Rate Limiter Describing Function', 'Color', 'w');

% Magnitude Plot
subplot(2, 1, 1);
semilogx(omega, abs(N_rle), 'b-', 'LineWidth', 2);
grid on;
ylabel('Magnitude |N_{rle}|');
title(sprintf('Rate Limiter Describing Function (R = %.1f, u_{rle} = %.1f, \\omega_{onset} = %.2f rad/s)', ...
    R, u_rle, omega_onset));
xlim([min(omega) max(omega)]);
ylim([0 1.1]);

% Phase Plot
subplot(2, 1, 2);
semilogx(omega, rad2deg(angle(N_rle)), 'r-', 'LineWidth', 2);
grid on;
xlabel('Frequency \omega (rad/s)');
ylabel('Phase (deg)');
xlim([min(omega) max(omega)]);

%% 5. Reusable Function Block
% You can copy this function to a separate file named 'calc_Nrle.m'
function [N_rle, omega_onset] = calc_Nrle(omega, u_rle, R)
% CALC_NRLE Computes the rate limiter describing function
%
% Inputs:
%   omega - Frequency (scalar or vector in rad/s)
%   u_rle - Input signal amplitude (scalar or vector)
%   R     - Rate limit value (scalar)
%
% Outputs:
%   N_rle       - Complex describing function value
%   omega_onset - Onset frequency of rate limiting

omega_onset = R ./ u_rle;

arg = (pi / 2) .* (omega_onset ./ omega);
N_rle = (4 / pi) .* (omega_onset ./ omega) .* exp(-1j * acos(arg));

% Ensure physical fidelity in the linear (un-limited) operating regime
N_rle(arg > 1) = 1.0;
end
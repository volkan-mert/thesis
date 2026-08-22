%% REGION III LIMIT CYCLE ANALYSIS FOR SELECTED Kp VALUES
%
% Fully developed rate saturation only
%
% Limit-cycle condition:
%
%       1 + G(jw)*N(A,w) = 0
%
% where
%
%       G(s) = Kp*Gc(s)*Gac(s)
%
% Equivalent conditions:
%
%       G(jw) = -1/N(A,w)
%
%       G(jw)*N(A,w) = -1
%
% Region III:
%
%       w >= w_crit
%
%       w_onset = R/A
%
%       w_crit = 1.862*w_onset

clear;
clc;
close all;

%% 1. Controller transfer function

numGc = [5.21, -273.7855, -1425.2, -700.1952];
denGc = [1, 21.3594, 545.5538, 605.6621, 0];

Gc = tf(numGc,denGc);

%% 2. Aircraft longitudinal dynamics

numGac = [-10.5240, -16.8384, -0.6247, 0];
denGac = [1, 2.3473, -5.3061, -0.1836, -0.0418];

Gac = tf(numGac,denGac);

%% 3. Selected Kp values

Kp_values = [0.1 0.5 1 5 10 15 20 50 100];

%% 4. Rate limiter parameter

R = 60;                 % Maximum actuator rate, deg/s

%% 5. Search ranges

A_min = 1;              % Minimum amplitude, deg
A_max = 500;            % Maximum amplitude, deg

w_min = 0.1;            % Minimum frequency, rad/s
w_max = 100;            % Maximum frequency, rad/s

A_values = linspace(A_min,A_max,500);
w_values = logspace(log10(w_min),log10(w_max),4000);

%% 6. Limit-cycle tolerance

error_limit = 1e-2;

%% 7. Initialize result variables

nKp = length(Kp_values);

Region = strings(nKp,1);
Kp_result = zeros(nKp,1);
A_result = zeros(nKp,1);
w_result = zeros(nKp,1);
w_onset_result = zeros(nKp,1);
w_crit_result = zeros(nKp,1);
alpha_result = zeros(nKp,1);
error_result = zeros(nKp,1);
Status = strings(nKp,1);

%% 8. Calculate Region III solution for each Kp

for m = 1:nKp

    Kp = Kp_values(m);

    % Linear system
    G = Kp*Gc*Gac;

    % Frequency response
    Gjw = squeeze(freqresp(G,w_values));
    Gjw = Gjw(:).';

    % Initial values
    best_error = inf;
    A0 = NaN;
    w0 = NaN;

    % Coarse search
    for i = 1:length(A_values)

        A_test = A_values(i);

        % Onset frequency
        w_onset_test = R/A_test;

        % Critical frequency
        w_crit_test = 1.862*w_onset_test;

        % Region III frequencies only
        region3 = w_values >= w_crit_test;

        if ~any(region3)
            continue;
        end

        w_test = w_values(region3);
        G_test = Gjw(region3);

        % Region III describing function
        varpi = w_onset_test./w_test;
        M = (4/pi).*varpi;
        phi = -acos((pi/2).*varpi);
        N_test = M.*exp(1j*phi);

        % Limit-cycle error
        error_test = abs(1 + G_test.*N_test);

        [error_now,index_now] = min(error_test);

        % Keep best candidate
        if error_now < best_error

            best_error = error_now;
            A0 = A_test;
            w0 = w_test(index_now);

        end

    end

    % Refine the result
    x0 = [log(w0) log(A0)];

    options = optimset('Display','off','TolX',1e-10,'TolFun',1e-12);

    x = fminsearch(@(x) region3_error(x,G,R,A_min,A_max,w_min,w_max),x0,options);

    % Final amplitude and frequency
    w_LC = exp(x(1));
    A_LC = exp(x(2));

    % Onset frequency
    w_onset_LC = R/A_LC;

    % Critical frequency
    w_crit_LC = 1.862*w_onset_LC;

    % Normalized frequency
    alpha_LC = w_LC*A_LC/R;

    % Region III describing function
    varpi_LC = w_onset_LC/w_LC;
    M_LC = (4/pi)*varpi_LC;
    phi_LC = -acos((pi/2)*varpi_LC);
    N_LC = M_LC*exp(1j*phi_LC);

    % Linear system at calculated frequency
    G_LC = squeeze(freqresp(G,w_LC));

    % Limit-cycle equation
    GN_LC = G_LC*N_LC;
    final_error = abs(1 + GN_LC);

    % Region III check
    region3_ok = w_LC >= w_crit_LC;

    % Decide if limit cycle is found
    if final_error < error_limit && region3_ok

        Status(m) = "LIMIT CYCLE FOUND";

    else

        Status(m) = "NO LIMIT CYCLE";

    end

    % Save results
    Region(m) = "Region III";
    Kp_result(m) = Kp;
    A_result(m) = A_LC;
    w_result(m) = w_LC;
    w_onset_result(m) = w_onset_LC;
    w_crit_result(m) = w_crit_LC;
    alpha_result(m) = alpha_LC;
    error_result(m) = final_error;

end

%% 9. Create result table

LimitCycleTable = table(Region,Kp_result,A_result,w_result,w_onset_result,w_crit_result,alpha_result,error_result,Status);

LimitCycleTable.Properties.VariableNames = {'Region','Kp','A_deg','w_rad_s','w_onset','w_crit','alpha','Error','Status'};

fprintf('\n');
fprintf('========================================================================================================\n');
fprintf('                                    REGION III RESULTS\n');
fprintf('========================================================================================================\n\n');

format short g

disp(LimitCycleTable);

%% 10. Create table containing only detected limit cycles

FoundLimitCycles = LimitCycleTable(LimitCycleTable.Status == "LIMIT CYCLE FOUND",:);

fprintf('\n');
fprintf('========================================================================================================\n');
fprintf('                                REGION III LIMIT CYCLES FOUND\n');
fprintf('========================================================================================================\n\n');

if isempty(FoundLimitCycles)

    fprintf('NO REGION III LIMIT CYCLES WERE FOUND.\n');

else

    disp(FoundLimitCycles);

end

%% 11. Nyquist and Nichols plots for detected limit cycles

for m = 1:height(FoundLimitCycles)

    % Current limit-cycle values
    Kp = FoundLimitCycles.Kp(m);
    A_LC = FoundLimitCycles.A_deg(m);
    w_LC = FoundLimitCycles.w_rad_s(m);
    w_onset_LC = FoundLimitCycles.w_onset(m);
    w_crit_LC = FoundLimitCycles.w_crit(m);

    % Linear system
    G = Kp*Gc*Gac;

    % Region III frequency range
    w_plot = logspace(log10(w_crit_LC),log10(w_max),4000);

    % Frequency response of G(jw)
    G_plot = squeeze(freqresp(G,w_plot));
    G_plot = G_plot(:).';

    % Region III describing function
    varpi_plot = w_onset_LC./w_plot;
    M_plot = (4/pi).*varpi_plot;
    phi_plot = -acos((pi/2).*varpi_plot);
    N_plot = M_plot.*exp(1j*phi_plot);

    % -1/N(A,w)
    minus_inv_N = -1./N_plot;

    % G(jw)*N(A,w)
    GN_plot = G_plot.*N_plot;

    % Convert responses to FRD models
    G_frd = frd(reshape(G_plot,1,1,length(w_plot)),w_plot);
    minus_inv_N_frd = frd(reshape(minus_inv_N,1,1,length(w_plot)),w_plot);
    GN_frd = frd(reshape(GN_plot,1,1,length(w_plot)),w_plot);

    % Values exactly at calculated limit-cycle frequency
    G_LC = squeeze(freqresp(G,w_LC));

    varpi_LC = w_onset_LC/w_LC;
    M_LC = (4/pi)*varpi_LC;
    phi_LC = -acos((pi/2)*varpi_LC);
    N_LC = M_LC*exp(1j*phi_LC);

    minus_inv_N_LC = -1/N_LC;
    GN_LC = G_LC*N_LC;

    %% Figure 1: Nyquist plot of G(jw) versus -1/N(A,w)

    figure;

    nyquist(G_frd,minus_inv_N_frd);

    grid on;
    hold on;

    % Mark the two points at w = w_LC
    plot(real(G_LC),imag(G_LC),'ko','MarkerFaceColor','k','MarkerSize',7);
    plot(real(minus_inv_N_LC),imag(minus_inv_N_LC),'rx','MarkerSize',10,'LineWidth',2);

    title(sprintf('Nyquist: G(j\\omega) vs -1/N(A,\\omega), K_p = %.1f, A = %.2f deg',Kp,A_LC));

    text(real(G_LC),imag(G_LC),sprintf('  \\omega = %.3f rad/s',w_LC),'FontWeight','bold');

    legend('G(j\omega)','-1/N(A,\omega)','G(j\omega_{LC})','-1/N(A,\omega_{LC})','Location','best');

    %% Figure 2: Nyquist plot of G(jw)*N(A,w)

    figure;

    nyquist(GN_frd);

    grid on;
    hold on;

    % Critical point
    plot(-1,0,'rx','MarkerSize',12,'LineWidth',2);

    % Calculated limit-cycle point
    plot(real(GN_LC),imag(GN_LC),'ko','MarkerFaceColor','k','MarkerSize',7);

    title(sprintf('Nyquist: G(j\\omega)N(A,\\omega), K_p = %.1f, A = %.2f deg',Kp,A_LC));

    text(real(GN_LC),imag(GN_LC),sprintf('  \\omega = %.3f rad/s',w_LC),'FontWeight','bold');

    legend('G(j\omega)N(A,\omega)','Critical Point (-1,0)','Limit Cycle','Location','best');

    %% Figure 3: Nichols chart of G(jw) versus -1/N(A,w)

    figure;

    nicholsplot(G_frd,minus_inv_N_frd);

    grid on;

    title(sprintf('Nichols: G(j\\omega) vs -1/N(A,\\omega), K_p = %.1f, A = %.2f deg',Kp,A_LC));

    legend('G(j\omega)','-1/N(A,\omega)','Location','best');

end

%% Local function: Region III limit-cycle error

function error = region3_error(x, G, R, A_min, A_max, w_min, w_max)

% Frequency and amplitude
w = exp(x(1));
A = exp(x(2));

% Check search limits
if A < A_min || A > A_max || w < w_min || w > w_max

    error = 1e6;
    return;

end

% Onset frequency
w_onset = R/A;

% Critical frequency
w_crit = 1.862*w_onset;

% Region III condition
if w < w_crit

    error = 1e6;
    return;

end

% Region III describing function
varpi = w_onset/w;
M = (4/pi)*varpi;
phi = -acos((pi/2)*varpi);
N = M*exp(1j*phi);

% Frequency response
Gjw = squeeze(freqresp(G,w));

% Limit-cycle condition
error = abs(1 + Gjw*N)^2;

end
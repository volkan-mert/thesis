%% AIRCRAFT LIMIT-CYCLE ANALYSIS USING RLE DESCRIBING FUNCTION
% Three-region actuator rate-limiter describing function.
%
% Linear system:
% G(s) = Kf * sys_ctrl * sys_dyn
%
% RLE regions are defined with
%
%       alpha = w / w_onset = Ai*w/R
%       w_onset = R/Ai
%
% Region I   : alpha < 1       -> no saturation
% Region II  : 1 <= alpha < 1.862 -> transition
% Region III : alpha >= 1.862  -> fully developed saturation
%
% Limit-cycle condition:
%
%       G(jw)*N(jw,Ai) = -1
%
% Ai is the sinusoidal input amplitude of the rate limiter.
% R is the maximum actuator surface rate.

clear; clc; close all;


%% 1. Forward Gain
Kf = 13.68;

%% 2. Control Law
% Controller used in the attached script.
numGc = [5.21,-273.7855,-1425.2,-700.1952];
denGc  = [1,21.3594,545.5538,605.6621,0];

Gc = tf(numGc,denGc);

sys_ctrl = Gc;

% sys_ctrl = tf(5.21*conv([1 -57.36],conv([1 4.26],[1 .55])), conv([1 2*.442*22.85 22.85^2],conv([1 0],[1 1.16])));


%% 3. Longitudinal Dynamics of the Aircraft
% Aircraft dynamics used in the attached script.

numGac = [-10.5240,-16.8384,-0.6247,0];
denGac = [1,2.3473,-5.3061,-0.1836,-0.0418];

Gac = tf(numGac,denGac);

sys_dyn = Gac;

% sys_dyn = tf(-10.524*conv([1 1.562],conv([1 .038],[1 0])), conv([1 2*.212*.088 .088^2],conv([1 3.75],[1 -44])));


%% 4. Complete Linear System
sys = Kf * sys_ctrl * sys_dyn;


%% 5. Rate-Limiter Parameters
% Maximum surface rate of the RLE.
% R = 40 deg/s is the value shown in the attached RLE example.
% Change this value if a different actuator rate limit is required.
R = 60;                         % deg/s

alpha_crit = 1.862;             % Region II / Region III boundary

% alpha = w/w_onset = Ai*w/R
% w_onset = R/Ai


%% 6. Frequency Response G(jw)
w = logspace(-3, 3, 30000);       % rad/s

G = squeeze(freqresp(sys,w));
G = G(:).';

phase_G = angle(G)*180/pi;

% Use the negative phase representation used in a Nichols chart.
phase_G(phase_G > 0) = phase_G(phase_G > 0) - 360;

mag_G_dB = 20*log10(abs(G));


%% 7. Three-Region RLE Describing Function
%
% REGION I: no saturation
%
%       alpha < 1
%       M = 1
%       phi = 0
%       N = 1
%
% REGION II: transition
%
%       1 <= alpha < 1.862
%
%       M = 0.2908*alpha^3 - 1.4396*alpha^2
%           + 1.9232*alpha + 0.223
%
%       phi = 0.528*alpha^3 - 2.6213*alpha^2
%             + 3.5056*alpha - 1.4171
%
%       N = M*exp(j*phi)
%
% The Region-II phase polynomial is in radians.
%
% REGION III: fully developed saturation
%
%       alpha >= 1.862
%
%       M = 4/(pi*alpha)
%
%       phi = -acos(pi/(2*alpha))
%
%       N = M*exp(j*phi)


% ----- Region I -----
N_I = 1;
minus_inv_N_I = -1/N_I;


% ----- Region II -----
alpha_II = linspace(1, alpha_crit, 3000);

M_II = 0.2908*alpha_II.^3 - 1.4396*alpha_II.^2 + 1.9232*alpha_II + 0.223;

phi_II = 0.528*alpha_II.^3 - 2.6213*alpha_II.^2 + 3.5056*alpha_II - 1.4171;

% The rounded polynomial coefficients do not give exactly M=1 and phi=0
% at alpha=1. Force the Region-I boundary to be exact.
M_II(1) = 1;
phi_II(1) = 0;

N_II = M_II .* exp(1i*phi_II);
minus_inv_N_II = -1 ./ N_II;


% ----- Region III -----
alpha_III = logspace(log10(alpha_crit), 3, 6000);

M_III = 4 ./ (pi*alpha_III);
phi_III = -acos(pi ./ (2*alpha_III));

N_III = M_III .* exp(1i*phi_III);
minus_inv_N_III = -1 ./ N_III;


%% 8. Find Limit Cycles in Regions II and III
%
% At a limit cycle:
%
%       G(jw)*N = -1
%
% Therefore
%
%       phase(G) + phase(N) = -180 deg
%
% and
%
%       |G|*|N| = 1
%
% The phase condition is first used to find alpha.
% The magnitude condition is then used to find w.


% Required RLE phase at every frequency.
phi_required_deg = -180 - phase_G;

% Residual = 0 dB at the limit cycle.
res_II  = nan(size(w));
res_III = nan(size(w));

alpha_req_II  = nan(size(w));
alpha_req_III = nan(size(w));


% ----- Region II candidates -----
phi_II_deg = phi_II*180/pi;

% Region II phase runs approximately from 0 deg to -32.6 deg.
valid_II = phi_required_deg <= 0 & phi_required_deg >= phi_II_deg(end);

if any(valid_II)

    % interp1 requires the phase vector in increasing order.
    alpha_req_II(valid_II) = interp1( phi_II_deg(end:-1:1), alpha_II(end:-1:1), phi_required_deg(valid_II), 'linear');

    a = alpha_req_II(valid_II);

    M = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.223;

    % Magnitude condition |G|*M = 1.
    res_II(valid_II) = mag_G_dB(valid_II) + 20*log10(M);
end


% ----- Region III candidates -----
% At the Region-III boundary:
phi_III_start_deg = phi_III(1)*180/pi;

% Region III phase approaches -90 deg as alpha increases.
valid_III = phi_required_deg <= phi_III_start_deg & phi_required_deg > -90;

if any(valid_III)

    phi_req_rad = phi_required_deg(valid_III)*pi/180;

    % From
    % phi = -acos(pi/(2*alpha))
    %
    % alpha = pi/(2*cos(phi))
    a = pi ./ (2*cos(phi_req_rad));

    alpha_req_III(valid_III) = a;

    M = 4 ./ (pi*a);

    % Magnitude condition |G|*M = 1.
    res_III(valid_III) = mag_G_dB(valid_III) + 20*log10(M);
end


% Result columns:
% 1  = Ai
% 2  = w
% 3  = alpha = w/w_onset
% 4  = M
% 5  = phi (deg)
% 6  = region number
% 7  = real(G*N)
% 8  = imag(G*N)
% 9  = real(G)
% 10 = imag(G)
% 11 = w_onset

result = [];


% Search Region II and Region III separately.
for region_now = 2:3

    if region_now == 2
        residual = res_II;
    else
        residual = res_III;
    end

    % Find sign changes of the magnitude residual.
    k_cross = find( ...
        isfinite(residual(1:end-1)) & ...
        isfinite(residual(2:end)) & ...
        residual(1:end-1).*residual(2:end) <= 0);

    for k = 1:length(k_cross)

        i1 = k_cross(k);
        i2 = i1 + 1;

        % Interpolate in log-frequency for a more accurate crossing.
        r1 = residual(i1);
        r2 = residual(i2);

        lw1 = log10(w(i1));
        lw2 = log10(w(i2));

        if abs(r2-r1) > eps
            lwc = lw1 + (0-r1)*(lw2-lw1)/(r2-r1);
        else
            lwc = 0.5*(lw1+lw2);
        end

        wc = 10^lwc;

        % Frequency response at the crossing.
        Gc = squeeze(freqresp(sys,wc));

        phase_Gc = angle(Gc)*180/pi;

        if phase_Gc > 0
            phase_Gc = phase_Gc - 360;
        end

        phi_req_deg = -180 - phase_Gc;


        if region_now == 2

            alpha_limit = interp1( phi_II_deg(end:-1:1), alpha_II(end:-1:1), phi_req_deg, 'linear');

            M_limit = 0.2908*alpha_limit^3 - 1.4396*alpha_limit^2 + 1.9232*alpha_limit + 0.223;

            phi_limit = 0.528*alpha_limit^3 - 2.6213*alpha_limit^2 + 3.5056*alpha_limit;         
        
        else

            phi_limit = phi_req_deg*pi/180;

            alpha_limit = pi/(2*cos(phi_limit));

            M_limit = 4/(pi*alpha_limit);

        end


        % Input amplitude of the sinusoid entering the rate limiter.
        %
        % alpha = Ai*w/R
        %
        % therefore
        %
        % Ai = alpha*R/w
        Ai_limit = alpha_limit*R/wc;

        w_onset = R/Ai_limit;

        N_limit = M_limit*exp(1i*phi_limit);
        L_limit = Gc*N_limit;


        % Do not save the same numerical crossing twice.
        duplicate = false;

        if ~isempty(result)

            same_w = abs(result(:,2)-wc) < 1e-4*max(1,wc);
            same_A = abs(result(:,1)-Ai_limit) < 1e-4*max(1,Ai_limit);

            duplicate = any(same_w & same_A);
        end

        if ~duplicate

            result = [result; Ai_limit, wc, alpha_limit, M_limit, phi_limit*180/pi, region_now, real(L_limit), imag(L_limit), real(Gc), imag(Gc), w_onset];
        end
    end
end


%% 9. Display Results
fprintf('\nRLE LIMIT-CYCLE RESULTS\n');
fprintf('=======================\n');
fprintf('Rate limit R = %.3f deg/s\n\n',R);

fprintf('Region definitions:\n');
fprintf('Region I   : alpha < 1      (no saturation, N = 1)\n');
fprintf('Region II  : 1 <= alpha < 1.862\n');
fprintf('Region III : alpha >= 1.862\n\n');

fprintf('Region I is the unsaturated linear case.\n');
fprintf('A rate-limit-induced nonlinear limit cycle is searched in Regions II and III.\n\n');


if isempty(result)

    fprintf('NO REGION-II OR REGION-III LIMIT CYCLE FOUND!\n\n');

else

    for k = 1:size(result,1)

        if result(k,6) == 2
            region_name = 'II - transition';
        else
            region_name = 'III - fully developed saturation';
        end

        fprintf('Candidate %d\n',k);
        fprintf('-----------\n');
        fprintf('Region       = %s\n',region_name);
        fprintf('Ai           = %.6f deg\n',result(k,1));
        fprintf('w            = %.6f rad/s\n',result(k,2));
        fprintf('f            = %.6f Hz\n',result(k,2)/(2*pi));
        fprintf('alpha        = %.6f\n',result(k,3));
        fprintf('w_onset      = %.6f rad/s\n',result(k,11));
        fprintf('M            = %.6f\n',result(k,4));
        fprintf('phi          = %.6f deg\n',result(k,5));
        fprintf('G(jw)N       = %.6f %+.6fj\n\n', ...
            result(k,7),result(k,8));
    end
end


%% 10. Nyquist Plot: G(jw) and -1/N
figure('Name','Nyquist Plot of G(jw) and RLE -1/N', ...
       'NumberTitle','off');

plot(real(G),imag(G),'LineWidth',1.5);
hold on;

% Region I is a single point because N = 1.
plot(real(minus_inv_N_I),imag(minus_inv_N_I), 'ks','MarkerFaceColor','k','MarkerSize',7);

% Region II and Region III negative inverse describing functions.
plot(real(minus_inv_N_II),imag(minus_inv_N_II), '--','LineWidth',1.3);

plot(real(minus_inv_N_III),imag(minus_inv_N_III), ':','LineWidth',1.5);


% Mark calculated limit-cycle intersections.
if ~isempty(result)

    for k = 1:size(result,1)

        Gc = result(k,9) + 1i*result(k,10);

        plot(real(Gc),imag(Gc), ...
             'ko','MarkerFaceColor','k','MarkerSize',7);

        text(real(Gc),imag(Gc), ...
            sprintf('  LC: A_i = %.3f, \\omega = %.3f', ...
            result(k,1),result(k,2)));
    end
end

grid on;
box on;

xlabel('Real');
ylabel('Imaginary');

title('RLE Describing-Function Limit-Cycle Condition');
subtitle('G(j\omega) = -1/N(j\omega,A_i)');

legend('G(j\omega)', 'Region I: -1/N', 'Region II: -1/N', 'Region III: -1/N', 'Location','best');


%% 11. Plot G(jw)N(jw,Ai) Around the Critical Point (-1,0)
figure('Name','Nyquist Plot of G(jw)N(jw,Ai)', ...
       'NumberTitle','off');

hold on;


if isempty(result)

    plot(-1,0,'ko','MarkerFaceColor','k');
    text(-1,0,'  (-1,0)');

    text(0,0, ...
        'NO REGION-II OR REGION-III LIMIT CYCLE FOUND!', ...
        'HorizontalAlignment','center');

else

    for k = 1:size(result,1)

        Ai = result(k,1);

        % alpha changes with frequency for a fixed input amplitude Ai.
        alpha_w = Ai*w/R;

        N_w = ones(size(w));


        % Region II
        id2 = alpha_w >= 1 & alpha_w < alpha_crit;

        a = alpha_w(id2);

        M = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.223;

        phi = 0.528*a.^3 - 2.6213*a.^2 + 3.5056*a - 1.4171;

        N_w(id2) = M .* exp(1i*phi);


        % Region III
        id3 = alpha_w >= alpha_crit;

        a = alpha_w(id3);

        M = 4 ./ (pi*a);
        phi = -acos(pi ./ (2*a));

        N_w(id3) = M .* exp(1i*phi);


        % Region I automatically remains N = 1.
        L = G .* N_w;

        plot(real(L),imag(L),'LineWidth',1.5);

        plot(result(k,7),result(k,8), 'ko','MarkerFaceColor','k','MarkerSize',7);

        text(result(k,7),result(k,8), sprintf('  A_i = %.3f, \\omega = %.3f', result(k,1),result(k,2)));
    end


    % Critical point.
    plot(-1,0,'ks','MarkerFaceColor','k');
    text(-1,0,'  (-1,0)');
end

grid on;
box on;

xlabel('Real');
ylabel('Imaginary');

title('Nyquist Plot of G(j\omega)N(j\omega,A_i)');
subtitle('Limit cycle occurs when G(j\omega)N(j\omega,A_i) = -1');


%% 12. Nichols Chart: G(jw) and the Three-Region RLE -1/N
%
% This is the Nichols-chart form of the negative inverse describing
% function technique.

figure('Name','Nichols Chart - Three-Region RLE', ...
       'NumberTitle','off');

hold on;


% Linear system Nichols coordinates.
plot(phase_G,mag_G_dB,'LineWidth',1.5);


% Region I: -1/N = -1.
plot(-180,0,'ks','MarkerFaceColor','k','MarkerSize',7);


% Region II Nichols coordinates.
phase_inv_II = angle(minus_inv_N_II)*180/pi;
phase_inv_II(phase_inv_II > 0) = phase_inv_II(phase_inv_II > 0) - 360;

mag_inv_II_dB = 20*log10(abs(minus_inv_N_II));

plot(phase_inv_II,mag_inv_II_dB,'--','LineWidth',1.5);


% Region III Nichols coordinates.
phase_inv_III = angle(minus_inv_N_III)*180/pi;
phase_inv_III(phase_inv_III > 0) = phase_inv_III(phase_inv_III > 0) - 360;

mag_inv_III_dB = 20*log10(abs(minus_inv_N_III));

plot(phase_inv_III,mag_inv_III_dB,':','LineWidth',1.5);


% Mark Region II / Region III boundary.
plot(phase_inv_III(1),mag_inv_III_dB(1), ...
     'ko','MarkerFaceColor','w','MarkerSize',7);

text(phase_inv_III(1),mag_inv_III_dB(1), ...
     '  \alpha = 1.862');


% Mark limit cycles.
if ~isempty(result)

    for k = 1:size(result,1)

        Gc = result(k,9) + 1i*result(k,10);

        phase_limit = angle(Gc)*180/pi;

        if phase_limit > 0
            phase_limit = phase_limit - 360;
        end

        mag_limit_dB = 20*log10(abs(Gc));

        plot(phase_limit,mag_limit_dB, 'ko','MarkerFaceColor','k','MarkerSize',7);

        text(phase_limit+2,mag_limit_dB, sprintf('limit cycle (%.3f rad/s)',result(k,2)), 'VerticalAlignment','bottom');
    end
end


grid on;
box on;

xlabel('phase, deg');
ylabel('amplitude, dB');

title('Nichols Chart: G(j\omega) and RLE Negative Inverse Describing Function');

legend('G(j\omega)', 'Region I', 'Region II: -1/N', 'Region III: -1/N', 'Region II/III boundary', 'Location','best');


%% 13. MATLAB Built-In Nichols Plot
figure('Name','MATLAB Nichols Plot of G(jw)', 'NumberTitle','off');

h = nicholsplot(sys);
setoptions(h,'PhaseMatching','on','Grid','on');

title('Nichols Plot of G(j\omega)');


%% 14. Open Control System Designer in Nichols View
controlSystemDesigner('nichols',sys);

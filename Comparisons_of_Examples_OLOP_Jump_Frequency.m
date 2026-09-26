%% Comparison of the OLOP Jump Frequency in PIO for the examples of X-15 Soft Glide in 1959 vs the Fatal Crash in 1967 X-15 Flight 3-65-97
clear; clc; close all

%% 1. Inputs
% COMMON PARAMETERS
n   = 10000;                % number of frequency points
w   = logspace(-1, 2, n);   % rad/s

% PARAMETERS OF THE SOFT GLIDE LANDING PIO INCIDENT
% ---o--- % : CTRL+R / CTRL+T to comment and uncomment
Kp_Yp   = 13.68;            % controller gain (here pilot gain has been inherited as an internal parameter of the X-15 Soft Glide Landing PIO Incident in 1959)
qco = 1.1;                  % pilot command amplitude of the X-15 Soft Glide Landing PIO Incident in 1959
R   = 15;                   % actuator rate limit of the X-15 Soft Glide Landing PIO Incident in 1959
Kp = 1;  % the pilot gain representing the gradient of the sidestick of the pilot model for the X-15 Soft Glide Landing PIO Incident in 1959

% Gc: Soft Glide PIO Incident's Transfer function of Y_p(s)
num_Gc  = Kp_Yp*tf(1, 1); 
den_Gc = 1;
Gc = tf(num_Gc, den_Gc);

% Gac: X-15 Soft Glide PIO Incident's Transfer function of theta/delta_e(s)
num_Gac = 0.537 * [1 0.82];
den_Gac = [1 1.42 5.29 0];
Gac = tf(num_Gac, den_Gac);

% PARAMETERS OF THE X-15 Flight 3-65-97, 1967 FATAL CRASH PIO INCIDENT
% ---o--- % : CTRL+R / CTRL+T to comment and uncomment
% R   = 60; % actuator rate limit of X-15 Flight 3-65-97, 1967 Fatal Crash PIO Incident's Transfer function of the Aircraft
% Kp = 13.68;                 % the pilot gain representing the gradient of the sidestick of the pilot model for X-15 Flight 3-65-97, 1967 Fatal Crash PIO Incident's Transfer function of the Aircraft
% qco = 1.1; of X-15 Flight 3-65-97, 1967 Fatal Crash PIO Incident's Transfer function of the Aircraft
% X-15 Flight 3-65-97, 1967 Fatal Crash PIO Incident's Transfer function of the Controller
% num_Gc = 5.21 * conv([1 -57.36], conv([1 4.26], [1 0.55])); 
% den_Gc = conv([1 2*0.442*22.85 22.85^2], conv([1 0], [1 1.16]));
% Gc  = tf(num_Gc, den_Gc);
% X-15 Flight 3-65-97, 1967 Fatal Crash PIO Incident's Transfer function of the Aircraft
% num_Gac = -10.524 * conv([1 1.562], conv([1 0.038], [1 0]));
% den_Gac = conv([1 2*0.212*0.088 0.088^2], conv([1 3.75], [1 -1.44]));
% Gac = tf(num_Gac, den_Gac);
% ---o--- %

% Frequency responses on the whole grid, computed once
Gcw  = squeeze(freqresp(Gc , w));
Gacw = squeeze(freqresp(Gac, w));

%% 2. Initial guess from the linear closed loop
pcl = Kp*Gc / (1 + Gc*Gac);
p1  = squeeze(freqresp(pcl, w(1)));
dGuess = qco * abs(p1);     % amplitude guess
fGuess = angle(p1);         % phase guess [rad]

%% 3. Harmonic-balance sweep
NoMag = zeros(1, n);        % open-loop DF magnitude [dB]
NoPh  = zeros(1, n);        % open-loop DF phase [rad]
res   = zeros(1, n);        % solver residual
kOn   = NaN;                % index of the onset point
kOns  = [];                 % array to store onset point indices

for k = 1:n
    mGc = abs(Gcw(k));   pGc = angle(Gcw(k));
    mGa = abs(Gacw(k));  pGa = angle(Gacw(k));
    
    % Solve the harmonic balance at this frequency
    cost = @(x) hbCost(x, w(k), R, Kp, qco, mGc, pGc, mGa, pGa);
    [x, res(k)] = fminsearch(cost, [dGuess fGuess], optimset('Display','off'));
    d = x(1);               % surface command amplitude
    f = x(2);               % surface command phase [rad]
    [mN, pN] = rateLimDF(w(k), R, d);
    
    % Closed loop -> open loop
    A   = d * mGa * mN / (Kp*qco);
    phi = f + pGa + pN;
    NoMag(k) = 20*log10( A / sqrt(1 - 2*A*cos(phi) + A^2));
    NoPh(k)  = phi - atan2(-A*sin(phi), 1 - A*cos(phi));
    
    % Next initial guess: continue, or restart if the solver did not converge
    if res(k) > 1e-4
        if isnan(kOn)
            kOn = k;
            kOns = [kOns, k];
            fprintf('Onset frequency: %.4g rad/s\n', w(k));
        end
        
        % Guesses of the Fatal Crash, X-15 Flight 3-65-97
        dGuess = 24;                % restart values for this aircraft
        fGuess = -237*pi/180;
        
        % Guesses of the Soft-Glide Landing PIO Incident
        % dGuess = 15;                
        % fGuess = -180*pi/180;
        
        % Default Guesses
        % dGuess = d;
        % fGuess = f;
    else
        dGuess = d;
        fGuess = f;
    end
end

%% 4. Nichols-type plot (plain plot, phase wrapped into [-360, 0] deg)
phLin  = mod(angle(Gcw.*Gacw)*180/pi, -360);
magLin = 20*log10(abs(Gcw.*Gacw));
phDF   = mod(NoPh*180/pi, -360);

% OLOP stability boundary
vPh  = [-60 -90 -100 -120 -140 -160 -180];
vMag = [13.5 7.5  5.5   2.5   1.1    0    0];

figure(1); hold on; grid on
plot(phLin, magLin, 'b-')
plot(phDF , NoMag , 'r*')
plot(phDF , NoMag , 'r--') % used only for showing the continuity of the line

% Plot all onset points and label them
if ~isempty(kOns)
    plot(phDF(kOns), NoMag(kOns), 'mp', 'MarkerSize', 14, 'MarkerFaceColor', 'm')
    for i = 1:length(kOns)
        idx = kOns(i);
        text(phDF(idx), NoMag(idx) + 1.5, sprintf('(%d.) \\omega_{onset} = %.2f rad/s', i, w(idx)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontWeight', 'bold', 'Color', 'm', 'FontSize', 10);
    end
end

% Plot stability boundary with updated LineWidth to match Figure 3
plot(vPh, vMag, 'k-', 'LineWidth', 4)

% Add stability boundary text
text(vPh(1)+20, vMag(1), 'Stability Boundary of OLOP', 'Color', 'k', ...
     'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);

% Critical point marker
plot(-180, 0, 'r+', 'MarkerSize', 8, 'LineWidth', 1.5)
xline(-180, 'k--'); yline(0, 'k--')

% Add Stable / Unstable text labels around the boundary 
text(-120, 15, 'Unstable', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);
text(-120, -15, 'Stable', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);

xlabel('Open-loop phase [deg]')
ylabel('Open-loop magnitude [dB]')
title('OLOP: linear loop vs. describing function')
legend('Linear G_c G_{ac}', 'OLOP DF', '', 'Onset point(s)', 'Stability boundary', ...
       'Location', 'best')

% Adjusted limits so the new text labels comfortably fit
xlim([-300 -20]); 
ylim([-25 25]);
hold off

%% 5. Convergence check
figure(2)
semilogx(w, res); grid on
xlabel('\omega [rad/s]'); ylabel('residual')
title('fminsearch residual')

%% 6. Nichols Plot with the Built-In nicholsplot()
figure(3)
% Linear open-loop transfer function
sys_open = Kp * Gc * Gac;
% Plot the linear system's Nichols chart
nicholsplot(sys_open, 'b');
hold on;
% Add classic Nichols grid (M and N circles) to the background
ngrid;
% Convert the non-linear OLOP magnitude (dB) and phase (rad) into a complex frequency response
resp_nl = 10.^(NoMag / 20) .* exp(1i * NoPh);
% Create the Frequency Response Data (FRD) object
sys_nl = frd(resp_nl, w);
% Plot the non-linear FRD model using nicholsplot with red asterisks
nicholsplot(sys_nl, 'r*');
nicholsplot(sys_nl, 'r--');  % this plot is only to show continuity by using a dashed line.

% Re-declare phase in degrees for standard plot markers (onset point)
phDF_deg = NoPh * 180 / pi; 

% Plot all onset points and number them on the Nichols Chart
if ~isempty(kOns)
    plot(phDF_deg(kOns), NoMag(kOns), 'mp', 'MarkerSize', 14, 'MarkerFaceColor', 'm')
    for i = 1:length(kOns)
        idx = kOns(i);
        text(phDF_deg(idx), NoMag(idx) + 1.5, sprintf('(%d.) \\omega_{onset} = %.2f rad/s', i, w(idx)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontWeight', 'bold', 'Color', 'm', 'FontSize', 10);
    end
end

% Add stability boundary and critical point
plot(vPh, vMag, 'k-', 'LineWidth', 4)
text(vPh(1)+20, vMag(1), 'Stability Boundary of OLOP', 'Color', 'k', ...
     'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);
plot(-180, 0, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5) 

% Add Stable / Unstable text labels around the boundary 
text(-120, 20, 'Unstable', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);
text(-120, -20, 'Stable', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);

title('Nichols Chart','The OLOP Jump Phenomena Analysis of the X-15''s Soft Glide Landing PIO Incident in 1959')
legend('Linear K_p*G_c*G_{ac}', 'OLOP DF', '', 'Onset point(s)', 'Stability boundary','', 'Location', 'best')
xlim auto
ylim auto
hold off

%% Local functions
function J = hbCost(x, freq, R, Kp, qco, mGc, pGc, mGa, pGa)
% Harmonic-balance residual (real and imaginary parts squared)
    d = x(1);
    f = x(2);
    [mN, pN] = rateLimDF(freq, R, d);
    
    re = d/mGc*cos(f - pGc) + d*mGa*mN*cos(f + pGa + pN) - Kp*qco;
    im = d/mGc*sin(f - pGc) + d*mGa*mN*sin(f + pGa + pN);
    J = re^2 + im^2;
end

function [mN, pN] = rateLimDF(freq, rate, amp)
% Duda's three-region describing function of a rate limiter.
%   mN : gain (1 below onset, 4/(pi*x) far above)
%   pN : phase in radians
    x = freq*amp/rate;
    if x < 1                            % Region I: no rate limiting
        mN = 1;
        pN = 0;
    elseif x < 1.862                    % Region II: polynomial fit
        mN = polyval([ 0.2908 -1.4396 1.9232  0.2230], x);
        pN = polyval([ 0.5280 -2.6213 3.5056 -1.4171], x);
    else                                % Region III: fully rate limited
        mN = 4/(pi*x);
        pN = -acos(pi/(2*x));
    end
end
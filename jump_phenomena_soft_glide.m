%% OLOP describing function for a rate-limited actuator (Scott Crossfield's X-15 Soft Glide Landing Incident in 1959)
%  Updated for Ashkenas (1964) PIO jump phenomena analysis with Kp correction
%  High-Resolution Version with multiple onset points
clear; clc; close all

%% 1. Inputs
n   = 10000;                % number of frequency points 
w   = logspace(-1, 2, n);   % rad/s
qco = 1.1;                  % pilot command amplitude [deg]

% Parameters derived from the Ashkenas X-15 analysis
R   = 15;                   % actuator rate limit [deg/s]
Kp_Yp  = 13.68;             % pilot gain given inside the transfer function as Kp of Y_{p} (s)
Kp = 1;                     % the pilot gain representing the gradient of the sidestick of the pilot model

% Gc = Transfer function of Y_p(s)
Gc  = Kp_Yp*tf(1, 1);             
          
% Gac = Transfer function of theta/delta_e(s)
num_Gac = 0.537 * [1 0.82];
den_Gac = [1 1.42 5.29 0];
Gac = tf(num_Gac, den_Gac);
          
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

% Variables for multiple onset detection
kOns   = [];                % Array to store all onset point indices
inJump = false;             % Flag to track contiguous jump regions

% Solver precision (resolution) increased
opt = optimset('Display','off', 'TolX', 1e-6, 'TolFun', 1e-6);

for k = 1:n
    mGc = abs(Gcw(k));   pGc = angle(Gcw(k));
    mGa = abs(Gacw(k));  pGa = angle(Gacw(k));
    
    % Solve the harmonic balance at this frequency
    cost = @(x) hbCost(x, w(k), R, Kp, qco, mGc, pGc, mGa, pGa);
    [x, res(k)] = fminsearch(cost, [dGuess fGuess], opt);
    
    d = x(1);               % surface command amplitude
    f = x(2);               % surface command phase [rad]
    [mN, pN] = rateLimDF(w(k), R, d);
    
    % Closed loop -> open loop 
    A   = d * mGa * mN / qco;
    phi = f + pGa + pN;
    NoMag(k) = 20*log10( A / sqrt(1 - 2*A*cos(phi) + A^2) );
    NoPh(k)  = phi - atan2(-A*sin(phi), 1 - A*cos(phi));
    
    % Next initial guess: continue, or restart if the solver did not converge
    if res(k) > 1e-4
        if ~inJump
            kOns(end+1) = k;
            inJump = true;
        end
        % dGuess = 15;                
        % fGuess = -180*pi/180;
        dGuess = d;
        fGuess = f;
    else
        inJump = false;
        dGuess = d;
        fGuess = f;
    end
end

% Print Onset Points and Stability Information to Command Window
if ~isempty(kOns)
    fprintf('\nDetected Onset Point(s):\n');
    for i = 1:length(kOns)
        fprintf('(%d.) Onset frequency: %.4g rad/s\n', i, w(kOns(i)));
    end
end
fprintf('\nNotice about the Stability of the OLOP:\n');
fprintf('If the OLOP is above 0 dB, rate saturation induces phase delay which amplifies closed-loop amplitude, creating a destabilizing snowball effect.\n');
fprintf('If the OLOP is below 0 dB, this amplitude growth is prevented, demonstrating that lower feedback gains mitigate Category II PIO risk.\n\n');

%% 4. Nichols Chart Plot (plain plot, phase wrapped into [-360, 0] deg)
phLin  = mod(angle(Gcw.*Gacw)*180/pi, -360);
magLin = 20*log10(abs(Kp*Gcw.*Gacw)); 
phDF   = mod(NoPh*180/pi, -360);

% OLOP stability boundary
vPh  = [-60 -90 -100 -120 -140 -160 -180];
vMag = [13.5 7.5  5.5   2.5   1.1    0    0];

figure(1); 
hold on; 
grid on
plot(phLin, magLin, 'b-', 'LineWidth', 1.5)
plot(phDF , NoMag , 'r*')
plot(phDF , NoMag , 'r--')  % this plot is only to show continuity by using a dashed line.

% Plot all onset points and number them
if ~isempty(kOns)
    plot(phDF(kOns), NoMag(kOns), 'mp', 'MarkerSize', 14, 'MarkerFaceColor', 'm')
    for i = 1:length(kOns)
        idx = kOns(i);
        text(phDF(idx), NoMag(idx) + 1.5, sprintf('(%d.) \\omega_{onset} = %.2f rad/s', i, w(idx)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'Color', 'm', 'FontSize', 10);  
    end
end

% Plot stability boundary and critical point
plot(vPh, vMag, 'k-', 'LineWidth', 4)
text(vPh(1)+20, vMag(1), 'Stability Boundary of OLOP', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);
plot(-180, 0, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5) 

% Add Stable / Unstable text labels around the boundary (-120 deg is a good anchor)
text(-120, 20, 'Unstable', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);
text(-120, -20, 'Stable', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);
xline(-180, 'k--'); yline(0, 'k--')
xlabel('Open-loop phase [deg]')
ylabel('Open-loop magnitude [dB]')
title('The Nichols Chart','The OLOP Jump Phenomena Analysis of the X-15''s Soft Glide Landing PIO Incident in 1959')
legend('Linear K_p*G_c*G_{ac}', 'OLOP DF', 'Onset point(s)', 'Stability boundary', 'Location', 'best')

axis auto 
xlim auto
ylim auto

hold off

%% 5. Nichols Plot with the Built-In nicholsplot()
figure(2)
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
        text(phDF_deg(idx), NoMag(idx) + 1.5, sprintf('(%d.) \\omega_{onset} = %.2f rad/s', i, w(idx)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'Color', 'm', 'FontSize', 10);
    end
end

% Add stability boundary and critical point
plot(vPh, vMag, 'k-', 'LineWidth', 4)
text(vPh(1)+20, vMag(1), 'Stability Boundary of OLOP', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);
plot(-180, 0, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5) 

% Add Stable / Unstable text labels around the boundary 
text(-120, 20, 'Unstable', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);
text(-120, -20, 'Stable', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontSize', 11);
title('Nichols Chart','The OLOP Jump Phenomena Analysis of the X-15''s Soft Glide Landing PIO Incident in 1959')
legend('Linear K_p*G_c*G_{ac}', 'OLOP DF', '', 'Onset point(s)', 'Stability boundary','', 'Location', 'best')

xlim auto
ylim auto

hold off

%% 6. Convergence check
figure(3)
semilogx(w, res, 'LineWidth', 1.2); grid on
xlabel('\omega [rad/s]'); ylabel('residual')
title('fminsearch residual')
axis tight; 


%% Figure Order Priority 
figure(2)


%% Local functions
function J = hbCost(x, freq, R, Kp, qco, mGc, pGc, mGa, pGa)
% Harmonic-balance residual (real and imaginary parts squared)
    d = x(1);
    f = x(2);
    [mN, pN] = rateLimDF(freq, R, d);
    re = d/mGc*cos(f - pGc) + d*Kp*mGa*mN*cos(f + pGa + pN) - Kp*qco;
    im = d/mGc*sin(f - pGc) + d*Kp*mGa*mN*sin(f + pGa + pN);
    J = re^2 + im^2;
end

function [mN, pN] = rateLimDF(freq, rate, amp)
% Duda's three-region describing function of a rate limiter.
    x = freq*amp/rate;
    if x < 1                            
        mN = 1;
        pN = 0;
    elseif x < 1.862                    
        mN = polyval([ 0.2908 -1.4396 1.9232  0.2230], x);
        pN = polyval([ 0.5280 -2.6213 3.5056 -1.4171], x);
    else                                
        mN = 4/(pi*x);
        pN = -acos(pi/(2*x));
    end
end
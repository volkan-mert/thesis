%% OLOP describing function for a rate-limited actuator (X-15 Flight 3-65-97, Michael J. Adams' Fatal Crash in 1967)
%  Based on Gilbreath (AFIT/GAE/ENY/01M-02) and Duda's rate-limiter
%  describing function. Everything is kept in radians to avoid deg/rad mix.

clear; clc; close all

%% 1. Inputs
n   = 2000;                 % number of frequency points
w   = logspace(-1, 2, n);   % rad/s
qco = 1.1;                  % pilot command amplitude
R   = 60;                   % actuator rate limit
K   = 13.68;                % pilot gain

Gc  = tf( 5.21 * conv([1 -57.36], conv([1 4.26], [1 0.55])), conv([1 2*0.442*22.85 22.85^2], conv([1 0], [1 1.16])) );

Gac = tf(-10.524 * conv([1 1.562], conv([1 0.038], [1 0])), conv([1 2*0.212*0.088 0.088^2], conv([1 3.75], [1 -1.44])) );

% Frequency responses on the whole grid, computed once
Gcw  = squeeze(freqresp(Gc , w));
Gacw = squeeze(freqresp(Gac, w));

%% 2. Initial guess from the linear closed loop
pcl = K*Gc / (1 + Gc*Gac);
p1  = squeeze(freqresp(pcl, w(1)));
dGuess = qco * abs(p1);     % amplitude guess
fGuess = angle(p1);         % phase guess [rad]

%% 3. Harmonic-balance sweep
NoMag = zeros(1, n);        % open-loop DF magnitude [dB]
NoPh  = zeros(1, n);        % open-loop DF phase [rad]
res   = zeros(1, n);        % solver residual
kOn   = NaN;                % index of the onset point

for k = 1:n

    mGc = abs(Gcw(k));   pGc = angle(Gcw(k));
    mGa = abs(Gacw(k));  pGa = angle(Gacw(k));

    % Solve the harmonic balance at this frequency
    cost = @(x) hbCost(x, w(k), R, K, qco, mGc, pGc, mGa, pGa);
    [x, res(k)] = fminsearch(cost, [dGuess fGuess], optimset('Display','off'));

    d = x(1);               % surface command amplitude
    f = x(2);               % surface command phase [rad]

    [mN, pN] = rateLimDF(w(k), R, d);

    % Closed loop -> open loop
    A   = d * mGa * mN / (K*qco);
    phi = f + pGa + pN;

    NoMag(k) = 20*log10( A / sqrt(1 - 2*A*cos(phi) + A^2) );
    NoPh(k)  = phi - atan2(-A*sin(phi), 1 - A*cos(phi));

    % Next initial guess: continue, or restart if the solver did not converge
    if res(k) > 1e-4
        if isnan(kOn)
            kOn = k;
            fprintf('Onset frequency: %.4g rad/s\n', w(k));
        end
        dGuess = 24;                % restart values for this aircraft
        fGuess = -237*pi/180;
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
plot(phDF , NoMag , 'r-')
if ~isnan(kOn)
    plot(phDF(kOn), NoMag(kOn), 'mp', 'MarkerSize', 14, 'MarkerFaceColor', 'm')
end
plot(vPh, vMag, 'k-', 'LineWidth', 2)
plot(-180, 0, 'r+')
xline(-180, 'k--'); yline(0, 'k--')

xlabel('Open-loop phase [deg]')
ylabel('Open-loop magnitude [dB]')
title('OLOP: linear loop vs. describing function')
legend('Linear G_c G_{ac}', 'OLOP DF', 'Onset point', 'Stability boundary', ...
       'Location', 'best')
xlim([-300 -50]); ylim([-20 20])
hold off

%% 5. Convergence check
figure(2)
semilogx(w, res); grid on
xlabel('\omega [rad/s]'); ylabel('residual')
title('fminsearch residual')

%% Local functions
function J = hbCost(x, freq, R, K, qco, mGc, pGc, mGa, pGa)
% Harmonic-balance residual (real and imaginary parts squared)
    d = x(1);
    f = x(2);
    [mN, pN] = rateLimDF(freq, R, d);

    re = d/mGc*cos(f - pGc) + d*mGa*mN*cos(f + pGa + pN) - K*qco;
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

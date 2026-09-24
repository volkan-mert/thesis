clear; clc; close all; clear functions

%% 1. Linear system G(s)

Kp = 13.68;        % Pilot gain
M_del_e = 0.537;     % Input signal entering the rate limiter
omega_n = 2.3;     % Natural frequency

zeta_sp = 1.42 / omega_n / 2;   % Damping ratio

num = Kp*M_del_e*[1 0.82];
den = [1 1.42 omega_n^2 0];

Gs = tf(num,den);

%% 2. Frequency range for G(jw)

w_G = linspace(0.3,3.5,1000);

%% 3. Rate Limiter Element parameters

R  = 15;           % Rate limit, deg/s
Ai = 13.68;        % Input amplitude, deg

w_N = linspace(0.01, 10, 1001);   % Frequencies, rad/s

w_onset = R / Ai;

%% 4. Describing function N(Ai,w)

N = zeros(size(w_N));

for k = 1:length(w_N)

    alpha = w_N(k) / w_onset;

    % Region I
    if alpha <= 1

        M = 1;
        phi = 0;

    % Region II
    elseif alpha < 1.862

        M = 0.2908*alpha^3 - 1.4396*alpha^2 + 1.9232*alpha + 0.223;

        phi = 0.528*alpha^3 - 2.6213*alpha^2 + 3.5056*alpha - 1.4171;

    % Region III
    else

        varpi = w_onset / w_N(k);

        M = (4/pi)*varpi;

        phi = -acos((pi/2)*varpi);

    end

    N(k) = M*exp(1j*phi);

end

%% 5. Calculate -1/N(Ai,w)

minus_inv_N = -1 ./ N;

%% 6. Convert -1/N to an FRD model

response = reshape(minus_inv_N,1,1,[]);

sys_N = frd(response,w_N);

%% 7. Nichols Chart

figure(Name='Nichols Chart',NumberTitle='off');

p1 = nicholsplot(Gs,w_G);

hold on

p2 = nicholsplot(sys_N);

%% 8. Shift -1/N from +180 deg to -180 deg

p2.PhaseMatchingEnabled = 'on';
p2.PhaseMatchingFrequency = w_N(1);

phase_first = rad2deg(angle(minus_inv_N(1)));

p2.PhaseMatchingValue = phase_first - 360;

%% 9. Figure settings

grid on

xlim([-180 -50])
ylim([-5 20])

yline(0,'--','DisplayName','0 dB');

title('Nichols Chart of G(j\omega) and -1/N(A_i,\omega)')

legend('G(j\omega)','-1/N(A_i,\omega)','0 dB')
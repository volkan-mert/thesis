clear;
clc;
close all;
clear functions


%% 1. Figure 8 - Ashkenas 1964 PIO Paper
% Linear system G(s)

Kp = 13.68;          % Pilot gain
M_del_e = 0.537;     % Input signal entering the rate limiter
omega_n = 2.3;       % Natural frequency

zeta_sp = 1.42 / omega_n / 2;   % Damping ratio

num = Kp*M_del_e*[1 0.82];
den = [1 1.42 omega_n^2 0];

Gs = tf(num,den);

Gs.Name = 'G(j\omega)';


%% 2. Frequency ranges

% Frequency range for the describing function
w_N = linspace(0.01,10,1001);

% Frequency range for G(jw)
w_G = linspace(0.3,3.5,1001);


%% 3. Rate Limiter Element parameters

R = 15;              % Rate limit, deg/s
Ai = 13.68;          % Input amplitude, deg

w_onset = R/Ai;


%% 4. Describing function N(Ai,w)

N = zeros(size(w_N));

for k = 1:length(w_N)

    alpha = w_N(k)/w_onset;


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

        varpi = w_onset/w_N(k);

        M = (4/pi)*varpi;

        phi = -acos((pi/2)*varpi);

    end


    N(k) = M*exp(1j*phi);

end


%% 5. Calculate -1/N(Ai,w)

minus_inv_N = -1./N;


%% 6. Convert -1/N(Ai,w) to an FRD model

response_N = reshape(minus_inv_N,1,1,length(w_N));

sys_N = frd(response_N,w_N);

sys_N.Name = '-1/N(A_i,\omega)';


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


%% 9. Nichols Chart settings

grid on

xlim([-180 -50])
ylim([-5 20])

yline(0,'--');

title('Nichols Chart of G(j\omega) and -1/N(A_i,\omega)')

legend('G(j\omega)','-1/N(A_i,\omega)','0 dB','Location','best')


%% 10. Nyquist Plot - Figure 5.24 Method
%
% G(jw) versus -1/N(A,w)
%
% G(jw):
% w changes from 0.3 to 3.5 rad/s
%
% -1/N(A,w):
% Each curve has a fixed frequency.
% Amplitude A changes along the curve.
%
% Limit-cycle condition:
%
% G(jw) = -1/N(A,w)
%
% The frequency of G(jw) at the intersection must be equal
% to the fixed frequency of the corresponding -1/N curve.


%% 10.1 Frequency response of G(jw)

G_response = freqresp(Gs,w_G);

G_response = reshape(G_response,1,1,length(w_G));

G_frd = frd(G_response,w_G);

G_frd.Name = 'G(j\omega)';


%% 10.2 Fixed frequencies for -1/N(A,w)

w_fixed = [0.3 0.6 1 2 3 3.5];


%% 10.3 Amplitude range

A_range = linspace(0.1,50,length(w_N));


%% 10.4 Store systems

systems_524 = cell(1,length(w_fixed)+1);

systems_524{1} = G_frd;


%% 10.5 Calculate -1/N(A,w) curves

for i = 1:length(w_fixed)

    w_fixed_i = w_fixed(i);

    N_curve = zeros(1,length(A_range));


    for k = 1:length(A_range)

        A = A_range(k);

        alpha = A*w_fixed_i/R;


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

            varpi = 1/alpha;

            M = (4/pi)*varpi;

            phi = -acos((pi/2)*varpi);

        end


        N_curve(k) = M*exp(1j*phi);

    end


    %% Calculate -1/N

    minus_inv_N_curve = -1./N_curve;


    %% Convert to FRD
    %
    % nyquistplot() requires an FRD frequency vector.
    %
    % w_N is used only to parameterize this curve for plotting.
    %
    % The physical frequency of the curve is w_fixed_i.
    % A_range is the quantity changing along the curve.

    response_curve = reshape(minus_inv_N_curve,1,1,length(w_N));

    sys_curve = frd(response_curve,w_N);

    sys_curve.Name = ...
        ['-1/N, \omega = ' num2str(w_fixed_i) ' rad/s'];


    systems_524{i+1} = sys_curve;

end


%% 10.6 Draw Figure 5.24 using nyquistplot()

figure(Name='Nyquist Plot - Figure 5.24 Method',NumberTitle='off');

nyquistplot(systems_524{:});

grid on

title('Nyquist Plot of G(j\omega) and -1/N(A,\omega)')

legend('show','Location','best')


%% 11. Nyquist Plot - Figure 5.25 Method
%
% G(jw)N(A,w)
%
% Each curve has a fixed amplitude A.
%
% Frequency changes from:
%
% 0.3 <= w <= 3.5 rad/s
%
% Limit-cycle condition:
%
% G(jw)N(A,w) = -1


%% 11.1 Frequency range

w = w_G;


%% 11.2 Calculate G(jw)
%
% reshape() is used so that Gjw is a 1 x 1001 row vector.
%
% This avoids the 1001 x 1001 matrix problem caused by
% multiplying a column vector by a row vector.

G_response = freqresp(Gs,w);

Gjw = reshape(G_response,1,length(w));


%% 11.3 Fixed amplitudes

A_values = [5 7.5 10 12 13 13.68 14 15 17.5 20 25];


%% 11.4 Store G(jw)N(A,w) systems

systems_525 = cell(1,length(A_values));


%% 11.5 Calculate G(jw)N(A,w)

for i = 1:length(A_values)

    A = A_values(i);

    N_curve = zeros(1,length(w));


    for k = 1:length(w)

        alpha = A*w(k)/R;


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

            varpi = 1/alpha;

            M = (4/pi)*varpi;

            phi = -acos((pi/2)*varpi);

        end


        N_curve(k) = M*exp(1j*phi);

    end


    %% Calculate G(jw)N(A,w)
    %
    % Gjw       : 1 x 1001
    % N_curve   : 1 x 1001
    % GN        : 1 x 1001

    GN = Gjw .* N_curve;


    %% Convert to FRD
    %
    % Required response dimensions:
    %
    % 1 x 1 x 1001

    response_GN = reshape(GN,1,1,length(w));

    sys_GN = frd(response_GN,w);

    sys_GN.Name = ['A = ' num2str(A) ' deg'];


    systems_525{i} = sys_GN;

end


%% 11.6 Draw Figure 5.25 using nyquistplot()

figure(Name='Nyquist Plot - Figure 5.25 Method',NumberTitle='off');

nyquistplot(systems_525{:});

grid on

title('Nyquist Plot of G(j\omega)N(A,\omega)')

legend('show','Location','best')


%% 12. Display basic information

fprintf('\n')
fprintf('ASHKENAS 1964 PIO EXAMPLE\n')
fprintf('--------------------------\n')

fprintf('Kp                    = %.3f\n',Kp)
fprintf('M_del_e               = %.3f\n',M_del_e)
fprintf('R                     = %.3f deg/s\n',R)
fprintf('Ai                    = %.3f deg\n',Ai)
fprintf('omega_onset           = %.4f rad/s\n',w_onset)

fprintf('\n')

fprintf('Frequency range for N(A,w):\n')
fprintf('w_N = %.2f to %.2f rad/s\n',w_N(1),w_N(end))

fprintf('\n')

fprintf('Frequency range for G(jw):\n')
fprintf('w_G = %.2f to %.2f rad/s\n',w_G(1),w_G(end))

fprintf('\n')

fprintf('Figure 5.24 limit-cycle condition:\n')
fprintf('G(jw) = -1/N(A,w)\n')

fprintf('\n')

fprintf('Figure 5.25 limit-cycle condition:\n')
fprintf('G(jw)N(A,w) = -1\n')

fprintf('\n')

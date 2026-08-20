clear; clc; close all

% Define the Laplace variable for continuous-time transfer functions
s = tf('s');

%% 1. Control Laws Transfer Function
% Shorthand: 5.21 (-57.36) (4.26) (.55) / (.442, 22.85) (1.16) (0.)

% Second-order parameters
zeta_ctrl = 0.442;
wn_ctrl = 22.85;

% Numerator and Denominator
num_ctrl = 5.21 * (s - 57.36) * (s + 4.26) * (s + 0.55);
den_ctrl = s * (s + 1.16) * (s^2 + 2*zeta_ctrl*wn_ctrl*s + wn_ctrl^2);

% Transfer Function Model
G_control = num_ctrl / den_ctrl;


%% 2. Aircraft Transfer Function
% Shorthand: -10.524 (1.562) (.038) (0.) / (.212, .088) (3.75) (-1.44)

% Second-order parameters
zeta_ac = 0.212;
wn_ac = 0.088;

% Numerator and Denominator
num_ac = -10.524 * s * (s + 1.562) * (s + 0.038);
den_ac = (s - 1.44) * (s + 3.75) * (s^2 + 2*zeta_ac*wn_ac*s + wn_ac^2);

% Transfer Function Model
G_aircraft = num_ac / den_ac;

%% 3. Open-loop Transfer Function

Kp = 1; % Pilot gain

G_oltf = Kp * G_control * G_aircraft;

%% 4. Stability Check via Poles (Direct Stability Check)
poles = pole(G_oltf);
disp('Poles of the Open-Loop System:');
disp(poles);

% Check if any pole is in the Right-Half Plane (RHP)
if any(real(poles) > 0)
    fprintf('The Open-Loop System is UNSTABLE (contains poles in the RHP).\n\n');
else
    fprintf('The Open-Loop System is STABLE.\n\n');
end

% Plot the Pole-Zero Map
figure(Name='Direct Stability Check: Pole-Zero Map of OLTF', NumberTitle='off');
pzmap(G_oltf);
grid on;
title('Pole-Zero Map of Open-Loop System');

%% 5. Frequency Domain Stability (Bode & Margins)
figure(Name='Stability Margins: Bode Plot', NumberTitle='off');
margin(G_oltf);
grid on;

% Extract margin values programmatically
[Gm, Pm, Wcg, Wcp] = margin(G_oltf);
Gm_dB = 20*log10(Gm);

fprintf('--- Stability Margins ---\n');
fprintf('Gain Margin: %.2f dB (at %.2f rad/s)\n', Gm_dB, Wcg);
fprintf('Phase Margin: %.2f deg (at %.2f rad/s)\n', Pm, Wcp);

%% 6. Nyquist Stability Plot
figure(Name='Nyquist Criterion', NumberTitle='off');
nyquist(G_oltf);
grid on;
title('Nyquist Plot');

%% 7. Root Locus Analysis for Pilot Gain (K_p)
% Define the baseline plant (G_control * G_aircraft) without pilot gain
G_plant = G_control * G_aircraft;

% Plot the Root Locus to see how Kp changes the pole locations
figure(Name='Root Locus Analysis', NumberTitle='off');
rlocus(G_plant);
grid on;
title('Root Locus for Pilot Gain (K_p) Design');

%% 8. Closed-Loop Stability Analysis
% Define a test gain based on what you observe in the Root Locus
Kp_test = 1; 

% Calculate Closed-Loop Transfer Function: T = (Kp*G) / (1 + Kp*G)
G_cltf = feedback(Kp_test * G_plant, 1);

cl_poles = pole(G_cltf);

fprintf('--- Closed-Loop Analysis (with Kp = %g) ---\n', Kp_test);
disp('Closed-Loop Poles:');
disp(cl_poles);

% Check Closed-Loop Stability
if any(real(cl_poles) > 0)
    fprintf('Result: The Closed-Loop system is still UNSTABLE with Kp = %g.\n', Kp_test);
    fprintf('Look at the Root Locus plot to find a gain value that pulls all branches into the left-half plane.\n');
else
    fprintf('Result: SUCCESS! The Closed-Loop system is STABLE with Kp = %g.\n', Kp_test);
    
    % If stable, plot the closed-loop step response
    figure(Name='Closed-Loop Stability Analysis', NumberTitle='off');
    step(G_cltf);
    grid on;
    title(['Closed-Loop Step Response (K_p = ', num2str(Kp_test), ')']);
    G_cltf
end
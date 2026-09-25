% 1. Define two arbitrary systems
s = tf('s');
sys1 = 10 / (s * (s + 1) * (s + 2));      % Curve 1 (Plant to be shifted)
sys2 = 1 / (s^2 + 0.5*s + 1);             % Curve 2 (Target boundary/plant)

% 2. Generate Frequency Response Data (FRD) objects
w = logspace(-1, 1, 1000); 
sys1_frd = frd(sys1, w);
sys2_frd = frd(sys2, w);

% Extract Magnitude (absolute) and Phase (degrees)
[mag1, phase1] = bode(sys1_frd, w);
[mag2, phase2] = bode(sys2_frd, w);

% Squeeze arrays and convert magnitude to dB
mag1_dB = 20*log10(squeeze(mag1));
mag2_dB = 20*log10(squeeze(mag2));

% Unwrap phase to prevent 360-degree jump artifacts during differentiation
p1_unwrapped = unwrap(deg2rad(squeeze(phase1))) * 180/pi;
p2_unwrapped = unwrap(deg2rad(squeeze(phase2))) * 180/pi;

% 3. Calculate Local Slopes (dGain / dPhase)
slope1 = diff(mag1_dB) ./ diff(p1_unwrapped);
slope2 = diff(mag2_dB) ./ diff(p2_unwrapped);

% 4. Find points of tangency (Matching Slopes)
% Let's arbitrarily pick a point on Curve 1 to touch Curve 2 (e.g., index 400)
idx1 = 400; 
target_slope = slope1(idx1);

% Find the point on Curve 2 that has the closest slope to the target
[~, idx2] = min(abs(slope2 - target_slope));

% 5. Calculate required shifts
Gain_Shift_dB = mag2_dB(idx2) - mag1_dB(idx1);

% Phase shift with wrapping to find the shortest path [-180, 180]
Phase_Shift_deg = p2_unwrapped(idx2) - p1_unwrapped(idx1);
Phase_Shift_deg = wrapTo180(Phase_Shift_deg); 

% 6. Apply shifts to the FRD object
% Convert Gain shift back to linear scalar, and Phase shift to complex exponential
K_linear = 10^(Gain_Shift_dB / 20);
Phase_complex = exp(1i * deg2rad(Phase_Shift_deg));

% Multiply the FRD object directly to apply constant gain and phase shift
sys1_shifted_frd = sys1_frd * (K_linear * Phase_complex);

% 7. Plot using nicholsplot()
figure;
nicholsplot(sys2_frd, 'b', sys1_shifted_frd, 'r--');
grid on;
title('Tangential Touching on Nichols Plot via Gain/Phase Shifting');
legend('Target Curve', 'Shifted Curve (Tangent)', 'Location', 'best');

fprintf('Applied Gain Shift: %.2f dB\n', Gain_Shift_dB);
fprintf('Applied Phase Shift: %.2f degrees\n', Phase_Shift_deg);
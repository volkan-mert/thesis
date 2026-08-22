clear; clc; close all

numGac = [-10.5240,-16.8384,-0.6247,0];
denGac = [1,2.3473,-5.3061,-0.1836,-0.0418];

Gac = tf(numGac,denGac);

numGc = [5.21,-273.7855,-1425.2,-700.1952];
denGc  = [1,21.3594,545.5538,605.6621,0];

Gc = tf(numGc,denGc);

Kp = 1;

sysCL = feedback(Kp*Gc*Gac, 1)

%% Stability Check

disp('Poles:')
p_ol = pole(sysCL)      % closed-loop poles

disp('Open-Loop Zeros:')
z_ol = zero(Gc*Gac)     % open-loop zeros of cascade

disp('Closed-Loop Zeros:')
z_cl = roots(numGc)     % zeros from numerator coeffs

disp('Closed-Loop Poles:')
p_cl = roots(denGc)     % poles from denominator coeffs

% pos_roots = 0;
% 
% for k = 1:length(p_cl)
%     if(p_cl(k) > 0 == 'true')
%         pos_roots = pos_roots + 1;
%     end
% end
% 
% if(pos_roots > 0)
%     disp('The closed-loop system is UNSTABLE!');
% else
%     disp('The closed-loop system is STABLE.');
% end

%  to represent only the observable and controllable dynamics, wrap your closed-loop system in the minreal() function to force the pole-zero cancellation:

sysCL_minimal = minreal(sysCL);

disp('Minimal Closed-Loop Poles:')
pole(sysCL_minimal)

pos_roots_m = 0;

for l = 1:length(p_cl)
    if(p_cl(l) > 0 == 'true')
        pos_roots_m = pos_roots_m + 1;
    end
end

if(pos_roots_m > 0)
    disp('The observable and controllable dynamics have been checked. The closed-loop system is UNSTABLE!');
else
    disp('The observable and controllable dynamics have been checked. The closed-loop system is STABLE.');
end


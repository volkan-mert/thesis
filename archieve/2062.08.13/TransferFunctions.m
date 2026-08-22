clear; clc; close all; clear function

%% TRANSFER FUNCTION OF GILBREATH 2001

K_glb = 13.68;
Gc_glb = tf(5.21*conv([1 -57.36],conv([1 4.26],[1 .55])),conv([1 2*.442*22.85 22.85^2],conv([1 0],[1 1.16])));
Gac_glb = tf(-10.524*conv([1 1.562],conv([1 .038],[1 0])),conv([1 2*.212*.088 .088^2],conv([1 3.75], [1 -1.44])));

num_glb = K_glb*Gc_glb*Gac_glb;
den_glb = zpk([], [], 1) + series(Gc_glb, Gac_glb);

pcl1 = tf(num_glb / den_glb)
pcl2 = K_glb*feedback(Gc_glb*Gac_glb, 1)
pcl3 = K_glb*Gc_glb*Gac_glb / (1 + Gc_glb*Gac_glb)

pcl1_minimal = minreal(pcl1)
pcl2_minimal = minreal(pcl2)
pcl3_minimal = minreal(pcl3)

% Extract and sort poles
p1 = sort(pole(pcl1_minimal))
p2 = sort(pole(pcl2_minimal))
p3 = sort(pole(pcl3_minimal))

% Extract and sort zeros (assuming they have the same number of zeros)
z1 = sort(zero(pcl1_minimal))
z2 = sort(zero(pcl2_minimal))
z3 = sort(zero(pcl3_minimal))
%
figure;

% Step Response Comparison
subplot(2,1,1);
step(pcl1_minimal, 'r-', pcl2_minimal, 'b--', pcl3_minimal, 'k:');
title('Step Response Comparison');
legend('pcl1', 'pcl2', 'pcl3');

% Bode Response Comparison
subplot(2,1,2);
bode(pcl1_minimal, 'r-', pcl2_minimal, 'b--', pcl3_minimal, 'k:');
title('Bode Plot Comparison');
legend('pcl1', 'pcl2', 'pcl3');
%% TRANSFER FUNCTION OF DUDA 1995 PAPER

Kp = 13.68;

numGac = [-10.5240,-16.8384,-0.6247,0];
denGac = [1, 2.3473, -5.3061, -0.1836, -0.0418];

numGc = [5.21, -273.7855, -1425.2, -700.1952];
denGc  = [1, 21.3594, 545.5538, 605.6621, 0];

TF_Gc = tf(numGc, denGc);
TF_Gac = tf(numGac, denGac);
TF_CL = Kp*feedback(TF_Gc*TF_Gac, 1);
TF_CL_minimal = minreal(TF_CL)

p4 = sort(pole(TF_CL_minimal))

z4 = sort(zero(TF_CL_minimal))
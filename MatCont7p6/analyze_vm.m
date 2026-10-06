clearvars -except 'exported'
close all;
clc
%% MatCont Limit-Cycle Data Extraction

X = exported.x;

nStates = 4;
nLC     = size(X,2);

% Number of points used to represent each periodic orbit
nMesh = (size(X,1)-2)/nStates;

fprintf('Number of limit cycles       = %d\n',nLC);
fprintf('Points per periodic orbit    = %d\n',nMesh);

%% Last two rows
T  = X(end-1,:);     % Period [s]
Kp = X(end,:);       % Active continuation parameter

%% Limit-cycle angular frequency
omega = 2*pi./T;     % [rad/s]

% Extract \(x_1,x_2,x_3,x_4\) for every limit cycle

thetaMax = zeros(1,nLC);
thetaMin = zeros(1,nLC);
thetaAmp = zeros(1,nLC);

x1Max = zeros(1,nLC);
x1Min = zeros(1,nLC);

for j = 1:nLC

    % One complete periodic orbit
    orbit = reshape( ...
        X(1:nStates*nMesh,j), ...
        nStates,nMesh);

    x1 = orbit(1,:);
    x2 = orbit(2,:);
    x3 = orbit(3,:);
    x4 = orbit(4,:);

    % Aircraft output
    theta = 6.02372*x2 + 7.346*x3;

    % Maximum/minimum values
    thetaMax(j) = max(theta);
    thetaMin(j) = min(theta);

    % Oscillation amplitude
    thetaAmp(j) = ...
        (thetaMax(j)-thetaMin(j))/2;

    x1Max(j) = max(x1);
    x1Min(j) = min(x1);
end

% Create a complete table

LCresults = table( ...
    (1:nLC)', ...
    Kp', ...
    T', ...
    omega', ...
    thetaMin', ...
    thetaMax', ...
    thetaAmp', ...
    'VariableNames', ...
    {'Point', ...
    'Kp', ...
    'Period_s', ...
    'Omega_rad_s', ...
    'ThetaMin', ...
    'ThetaMax', ...
    'ThetaAmplitude'});

disp(LCresults);

% Find the point closest to 2.8 rad/s

targetOmega = 2.8;

[errorOmega,idx] = min(abs(omega-targetOmega));

fprintf('\n====================================\n');
fprintf('LC closest to omega = %.3f rad/s\n',targetOmega);
fprintf('====================================\n');

fprintf('Continuation point : %d\n',idx);
fprintf('Kp                 : %.8f\n',Kp(idx));
fprintf('Period             : %.8f s\n',T(idx));
fprintf('Omega              : %.8f rad/s\n',omega(idx));
fprintf('Theta minimum      : %.8f\n',thetaMin(idx));
fprintf('Theta maximum      : %.8f\n',thetaMax(idx));
fprintf('Theta amplitude    : %.8f\n',thetaAmp(idx));
fprintf('Frequency error    : %.8f rad/s\n',errorOmega);

% Find the LC closest to Kp = 13.68

targetKp = 13.68;

[errorKp,idxKp] = min(abs(Kp-targetKp));

fprintf('\n====================================\n');
fprintf('LC closest to Kp = %.3f\n',targetKp);
fprintf('====================================\n');

fprintf('Continuation point : %d\n',idxKp);
fprintf('Kp                 : %.8f\n',Kp(idxKp));
fprintf('Period             : %.8f s\n',T(idxKp));
fprintf('Omega              : %.8f rad/s\n',omega(idxKp));
fprintf('Theta minimum      : %.8f\n',thetaMin(idxKp));
fprintf('Theta maximum      : %.8f\n',thetaMax(idxKp));
fprintf('Theta amplitude    : %.8f\n',thetaAmp(idxKp));
fprintf('Kp error           : %.8f\n',errorKp);

% Plot the actual LC bifurcation branch

figure;

plot(Kp,thetaMax,'LineWidth',1.5);
hold on;
plot(Kp,thetaMin,'LineWidth',1.5);

xlabel('K_p');
ylabel('\theta');
title('Limit-Cycle Bifurcation Diagram');
legend('\theta_{max}','\theta_{min}', ...
    'Location','best');
grid on;

% Plot the LC frequency

figure;

plot(Kp,omega,'LineWidth',1.5);

xlabel('K_p');
ylabel('\omega_{LC} [rad/s]');
title('Limit-Cycle Frequency vs. K_p');
grid on;


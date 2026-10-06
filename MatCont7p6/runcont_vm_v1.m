%% X-15 PIO BIFURCATION ANALYSIS USING MATCONT7P6
% Nonlinear model:
%
% x1_dot = S*tanh(K*(Kp*(thetac-theta)-x1)/S)
% x2_dot = x3
% x3_dot = x4
% x4_dot = x1 - 5.29*x3 - 1.42*x4
%
% Aircraft output:
%
% theta     = 6.02372*x2 + 7.346*x3
% theta_dot = 6.02372*x3 + 7.346*x4
%
% MatCont parameters:
%
% p(1) = Kp
% p(2) = K
% p(3) = S
% p(4) = thetac
%
% Continuation parameter: Kp

clear;
clc;
close all;

%% 1. LATEX FORMAT
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 2. MATCONT PATH
matcontPath = 'C:\Users\t0900\Documents\MATLAB\Volkan\MatCont7p6';

if ~isfolder(matcontPath)
    error('MatCont7p6 folder was not found.');
end

addpath(genpath(matcontPath));

fprintf('\n============================================================\n');
fprintf(' X-15 PIO MATCONT7P6 ANALYSIS\n');
fprintf('============================================================\n');
fprintf('MatCont path:\n%s\n',matcontPath);

%% 3. CHECK MATCONT
if exist('cont','file') ~= 2
    error('cont.m could not be found.');
end

if exist('init_EP_EP','file') ~= 2
    error('init_EP_EP.m could not be found.');
end

if exist('init_H_LC','file') ~= 2
    error('init_H_LC.m could not be found.');
end

fprintf('MatCont successfully found.\n');

%% 4. MATCONT SYSTEM FILE
systemName = 'X15PIOSoftGlide';

if exist(systemName,'file') ~= 2
    error('X15PIOSoftGlide.m could not be found.');
end

odefile = @X15PIOSoftGlide;

fprintf('System file:\n%s\n',which(systemName));

%% 5. MODEL PARAMETERS
Kp0 = 1.0;
K = 20.0;
S = 15.0;
thetac = 0.0;

% Parameter vector:
% p(1) = Kp
% p(2) = K
% p(3) = S
% p(4) = thetac

p = [Kp0 K S thetac];

% Active continuation parameter = Kp
ap = 1;

fprintf('\nInitial parameters:\n');
fprintf('Kp     = %.4f\n',Kp0);
fprintf('K      = %.4f\n',K);
fprintf('S      = %.4f deg/s\n',S);
fprintf('thetac = %.4f deg\n',thetac);

%% 6. INITIAL EQUILIBRIUM
xEq0 = [0;0;0;0];

nStates = length(xEq0);

%% 7. AIRCRAFT OUTPUT EQUATIONS
Ctheta = [0 6.02372 7.346 0];

CthetaDot = [0 0 6.02372 7.346];

%% 8. INITIALIZE EQUILIBRIUM CONTINUATION
[x0EP,v0EP] = init_EP_EP(odefile,xEq0,p,ap);

%% 9. EQUILIBRIUM CONTINUATION OPTIONS
optEP = contset;

optEP = contset(optEP,'Singularities',1);
optEP = contset(optEP,'Eigenvalues',1);
optEP = contset(optEP,'MaxNumPoints',600);
optEP = contset(optEP,'InitStepsize',1e-3);
optEP = contset(optEP,'MinStepsize',1e-6);
optEP = contset(optEP,'MaxStepsize',0.02);
optEP = contset(optEP,'FunTolerance',1e-7);
optEP = contset(optEP,'VarTolerance',1e-7);
optEP = contset(optEP,'TestTolerance',1e-6);
optEP = contset(optEP,'Backward',0);

%% 10. EQUILIBRIUM CONTINUATION
fprintf('\n============================================================\n');
fprintf(' EQUILIBRIUM CONTINUATION\n');
fprintf('============================================================\n');

[xEP,vEP,sEP,hEP,fEP] = cont(@equilibrium,x0EP,v0EP,optEP);

%% 11. DISPLAY EQUILIBRIUM SPECIAL POINTS
fprintf('\nEquilibrium special points:\n');

for k = 1:length(sEP)
    label = strtrim(sEP(k).label);

    if isempty(label)
        label = '-';
    end

    fprintf('%2d   Label = %-5s   Index = %d\n',k,label,sEP(k).index);
end

%% 12. FIND HOPF BIFURCATION
HopfNumber = [];

for k = 1:length(sEP)
    if strcmp(strtrim(sEP(k).label),'H')
        HopfNumber = k;
        break;
    end
end

%% 13. TRY BACKWARD EQUILIBRIUM CONTINUATION IF NEEDED
if isempty(HopfNumber)
    fprintf('\nNo Hopf point found in forward direction.\n');
    fprintf('Trying backward equilibrium continuation...\n');

    optEPback = contset(optEP,'Backward',1);

    [xEP,vEP,sEP,hEP,fEP] = cont(@equilibrium,x0EP,v0EP,optEPback);

    for k = 1:length(sEP)
        if strcmp(strtrim(sEP(k).label),'H')
            HopfNumber = k;
            break;
        end
    end
end

if isempty(HopfNumber)
    error('No Hopf bifurcation was found.');
end

%% 14. EXTRACT HOPF POINT
iH = sEP(HopfNumber).index;

xH = xEP(1:nStates,iH);

KpH = xEP(end,iH);

pH = p;

pH(ap) = KpH;

thetaH = Ctheta*xH;

fprintf('\n============================================================\n');
fprintf(' HOPF BIFURCATION\n');
fprintf('============================================================\n');
fprintf('Hopf index = %d\n',iH);
fprintf('Kp_H = %.10f\n',KpH);
fprintf('x1 = %.6e\n',xH(1));
fprintf('x2 = %.6e\n',xH(2));
fprintf('x3 = %.6e\n',xH(3));
fprintf('x4 = %.6e\n',xH(4));

%% 15. LINEARIZED SYSTEM AT HOPF
AH = [-K -6.02372*K*KpH -7.346*K*KpH 0;
0 0 1 0;
0 0 0 1;
1 0 -5.29 -1.42];

lambdaH = eig(AH);

fprintf('\nEigenvalues at Hopf:\n');

disp(lambdaH);

%% 16. HOPF FREQUENCY
[~,index] = max(abs(imag(lambdaH)));

omegaHopf = abs(imag(lambdaH(index)));

THopf = 2*pi/omegaHopf;

fprintf('omega_H = %.6f rad/s\n',omegaHopf);
fprintf('T_H     = %.6f s\n',THopf);

%% 17. INITIALIZE LIMIT CYCLE
ntst = 40;

ncol = 4;

LCamp = 1e-4;

try
    [x0LC,v0LC] = init_H_LC(odefile,xH,pH,ap,LCamp,ntst,ncol);
catch
    LCamp = 1e-3;

    fprintf('Retrying LC initialization with amplitude %.2e\n',LCamp);

    [x0LC,v0LC] = init_H_LC(odefile,xH,pH,ap,LCamp,ntst,ncol);
end

%% 18. LIMIT-CYCLE CONTINUATION OPTIONS
optLC = contset;

optLC = contset(optLC,'Singularities',1);
optLC = contset(optLC,'Multipliers',1);
optLC = contset(optLC,'MaxNumPoints',600);
optLC = contset(optLC,'InitStepsize',1e-4);
optLC = contset(optLC,'MinStepsize',1e-7);
optLC = contset(optLC,'MaxStepsize',0.02);
optLC = contset(optLC,'FunTolerance',1e-7);
optLC = contset(optLC,'VarTolerance',1e-7);
optLC = contset(optLC,'TestTolerance',1e-6);
optLC = contset(optLC,'Adapt',1);

%% 19. FORWARD LIMIT-CYCLE CONTINUATION
fprintf('\n============================================================\n');
fprintf(' FORWARD LIMIT-CYCLE CONTINUATION\n');
fprintf('============================================================\n');

optForward = contset(optLC,'Backward',0);

[xLCf,vLCf,sLCf,hLCf,fLCf] = cont(@limitcycle,x0LC,v0LC,optForward);

%% 20. BACKWARD LIMIT-CYCLE CONTINUATION
fprintf('\n============================================================\n');
fprintf(' BACKWARD LIMIT-CYCLE CONTINUATION\n');
fprintf('============================================================\n');

optBackward = contset(optLC,'Backward',1);

try
    [xLCb,vLCb,sLCb,hLCb,fLCb] = cont(@limitcycle,x0LC,v0LC,optBackward);
catch ME
    fprintf('Backward LC continuation could not be completed.\n');
    fprintf('%s\n',ME.message);

    xLCb = [];
    vLCb = [];
    sLCb = [];
    hLCb = [];
    fLCb = [];
end

%% 21. PROCESS LIMIT-CYCLE DATA
resultF = processLCBranch(xLCf,nStates,Ctheta,CthetaDot,K,S,thetac);

if ~isempty(xLCb)
    resultB = processLCBranch(xLCb,nStates,Ctheta,CthetaDot,K,S,thetac);
else
    resultB = [];
end

%% 22. EQUILIBRIUM BRANCH
KpEP = xEP(end,:);

thetaEP = zeros(size(KpEP));

for k = 1:size(xEP,2)
    xNow = xEP(1:nStates,k);

    thetaEP(k) = Ctheta*xNow;
end

%% 23. FIND LIMIT CYCLE CLOSEST TO omega = 2.8 rad/s
targetOmega = 2.8;

[indexF,errorF] = closestPoint(resultF.omega,targetOmega);

if ~isempty(resultB)
    [indexB,errorB] = closestPoint(resultB.omega,targetOmega);
else
    indexB = [];
    errorB = inf;
end

if isempty(indexF) && isempty(indexB)
    error('No valid limit-cycle frequency was found.');
end

if errorF <= errorB
    selected = resultF;
    selectedIndex = indexF;
    selectedBranch = 'Forward';
else
    selected = resultB;
    selectedIndex = indexB;
    selectedBranch = 'Backward';
end

%% 24. SELECTED LIMIT CYCLE
KpSelected = selected.Kp(selectedIndex);

TSelected = selected.period(selectedIndex);

fSelected = selected.frequency(selectedIndex);

omegaSelected = selected.omega(selectedIndex);

thetaAmpSelected = selected.thetaAmp(selectedIndex);

thetaDotAmpSelected = selected.thetaDotAmp(selectedIndex);

thetaSelected = selected.thetaOrbit{selectedIndex};

thetaDotSelected = selected.thetaDotOrbit{selectedIndex};

x1dotSelected = selected.x1dotOrbit{selectedIndex};

XSelected = selected.stateOrbit{selectedIndex};

x1Selected = XSelected(1,:);

fprintf('\n============================================================\n');
fprintf(' SELECTED LIMIT CYCLE\n');
fprintf('============================================================\n');
fprintf('Branch      = %s\n',selectedBranch);
fprintf('Kp          = %.6f\n',KpSelected);
fprintf('T           = %.6f s\n',TSelected);
fprintf('f           = %.6f Hz\n',fSelected);
fprintf('omega       = %.6f rad/s\n',omegaSelected);
fprintf('A_theta     = %.6f deg\n',thetaAmpSelected);
fprintf('A_theta_dot = %.6f deg/s\n',thetaDotAmpSelected);

%% 25. TIME VECTOR FOR DISPLAY
t = linspace(0,TSelected,length(thetaSelected));

%% 26. FIGURE SCREEN LAYOUT
totalFigures = 14;

figureRows = 4;

figureColumns = 4;

%% 27. FIGURE 1 : BIFURCATION DIAGRAM
[fig1,ax1] = makeFigure('01 - Bifurcation Diagram',1,totalFigures,figureRows,figureColumns);

hold(ax1,'on');

plot(ax1,KpEP,thetaEP,'k-','LineWidth',1.5,'DisplayName','$\mathrm{Equilibrium}$');

plot(ax1,resultF.Kp,resultF.thetaMax,'LineWidth',1.8,'DisplayName','$\theta_{\max}$');

plot(ax1,resultF.Kp,resultF.thetaMin,'LineWidth',1.8,'DisplayName','$\theta_{\min}$');

if ~isempty(resultB)
    plot(ax1,resultB.Kp,resultB.thetaMax,'--','LineWidth',1.5,'DisplayName','$\theta_{\max}$, backward');

    plot(ax1,resultB.Kp,resultB.thetaMin,'--','LineWidth',1.5,'DisplayName','$\theta_{\min}$, backward');
end

plot(ax1,KpH,thetaH,'ko','MarkerFaceColor','k','MarkerSize',6,'DisplayName','$H$');

xlabel(ax1,'$K_p$');

ylabel(ax1,'$\theta\;[\mathrm{deg}]$');

title(ax1,'$\mathrm{X\!-\!15\ PIO\ Bifurcation\ Diagram}$');

legend(ax1,'Location','best');

finishFigure(fig1,ax1);

%% 28. FIGURE 2 : LIMIT-CYCLE AMPLITUDE
[fig2,ax2] = makeFigure('02 - Limit-Cycle Amplitude',2,totalFigures,figureRows,figureColumns);

hold(ax2,'on');

plot(ax2,resultF.Kp,resultF.thetaAmp,'LineWidth',1.8,'DisplayName','$\mathrm{Forward\ LC}$');

if ~isempty(resultB)
    plot(ax2,resultB.Kp,resultB.thetaAmp,'--','LineWidth',1.5,'DisplayName','$\mathrm{Backward\ LC}$');
end

plot(ax2,KpH,0,'ko','MarkerFaceColor','k','MarkerSize',6,'DisplayName','$H$');

plot(ax2,KpSelected,thetaAmpSelected,'ro','MarkerFaceColor','r','MarkerSize',6,'DisplayName','$\mathrm{Selected\ LC}$');

xlabel(ax2,'$K_p$');

ylabel(ax2,'$A_{\theta}\;[\mathrm{deg}]$');

title(ax2,'$\mathrm{Limit\!-\!Cycle\ Amplitude}\quad A_{\theta}(K_p)$');

legend(ax2,'Location','best');

finishFigure(fig2,ax2);

%% 29. FIGURE 3 : LIMIT-CYCLE PERIOD
[fig3,ax3] = makeFigure('03 - Limit-Cycle Period',3,totalFigures,figureRows,figureColumns);

hold(ax3,'on');

plot(ax3,resultF.Kp,resultF.period,'LineWidth',1.8,'DisplayName','$\mathrm{Forward\ LC}$');

if ~isempty(resultB)
    plot(ax3,resultB.Kp,resultB.period,'--','LineWidth',1.5,'DisplayName','$\mathrm{Backward\ LC}$');
end

plot(ax3,KpSelected,TSelected,'ro','MarkerFaceColor','r','MarkerSize',6,'DisplayName','$\mathrm{Selected\ LC}$');

xlabel(ax3,'$K_p$');

ylabel(ax3,'$T\;[\mathrm{s}]$');

title(ax3,'$\mathrm{Limit\!-\!Cycle\ Period}\quad T(K_p)$');

legend(ax3,'Location','best');

finishFigure(fig3,ax3);

%% 30. FIGURE 4 : ANGULAR FREQUENCY
[fig4,ax4] = makeFigure('04 - Angular Frequency',4,totalFigures,figureRows,figureColumns);

hold(ax4,'on');

plot(ax4,resultF.Kp,resultF.omega,'LineWidth',1.8,'DisplayName','$\mathrm{Forward\ LC}$');

if ~isempty(resultB)
    plot(ax4,resultB.Kp,resultB.omega,'--','LineWidth',1.5,'DisplayName','$\mathrm{Backward\ LC}$');
end

plot(ax4,KpH,omegaHopf,'ko','MarkerFaceColor','k','MarkerSize',6,'DisplayName','$\omega_H$');

plot(ax4,KpSelected,omegaSelected,'ro','MarkerFaceColor','r','MarkerSize',6,'DisplayName','$\mathrm{Selected\ LC}$');

yline(ax4,targetOmega,'--','$\omega=2.8\;\mathrm{rad/s}$','Interpreter','latex','HandleVisibility','off');

xlabel(ax4,'$K_p$');

ylabel(ax4,'$\omega\;[\mathrm{rad/s}]$');

title(ax4,'$\mathrm{Limit\!-\!Cycle\ Angular\ Frequency}\quad \omega(K_p)$');

legend(ax4,'Location','best');

finishFigure(fig4,ax4);

%% 31. FIGURE 5 : PITCH-RATE AMPLITUDE
[fig5,ax5] = makeFigure('05 - Pitch-Rate Amplitude',5,totalFigures,figureRows,figureColumns);

hold(ax5,'on');

plot(ax5,resultF.Kp,resultF.thetaDotAmp,'LineWidth',1.8,'DisplayName','$\mathrm{Forward\ LC}$');

if ~isempty(resultB)
    plot(ax5,resultB.Kp,resultB.thetaDotAmp,'--','LineWidth',1.5,'DisplayName','$\mathrm{Backward\ LC}$');
end

plot(ax5,KpSelected,thetaDotAmpSelected,'ro','MarkerFaceColor','r','MarkerSize',6,'DisplayName','$\mathrm{Selected\ LC}$');

xlabel(ax5,'$K_p$');

ylabel(ax5,'$A_{\dot{\theta}}\;[\mathrm{deg/s}]$');

title(ax5,'$\mathrm{Pitch\!-\!Rate\ Limit\!-\!Cycle\ Amplitude}$');

legend(ax5,'Location','best');

finishFigure(fig5,ax5);

%% 32. FIGURE 6 : MAXIMUM ACTUATOR RATE
[fig6,ax6] = makeFigure('06 - Maximum Actuator Rate',6,totalFigures,figureRows,figureColumns);

hold(ax6,'on');

plot(ax6,resultF.Kp,resultF.maxX1dot,'LineWidth',1.8,'DisplayName','$\mathrm{Forward\ LC}$');

if ~isempty(resultB)
    plot(ax6,resultB.Kp,resultB.maxX1dot,'--','LineWidth',1.5,'DisplayName','$\mathrm{Backward\ LC}$');
end

yline(ax6,S,'--','$S=15\;\mathrm{deg/s}$','Interpreter','latex','HandleVisibility','off');

xlabel(ax6,'$K_p$');

ylabel(ax6,'$\max|\dot{x}_1|\;[\mathrm{deg/s}]$');

title(ax6,'$\mathrm{Maximum\ Actuator\ Rate}$');

legend(ax6,'Location','best');

finishFigure(fig6,ax6);

%% 33. FIGURE 7 : PHASE PORTRAIT
[fig7,ax7] = makeFigure('07 - Phase Portrait',7,totalFigures,figureRows,figureColumns);

plot(ax7,thetaSelected,thetaDotSelected,'LineWidth',1.8);

xlabel(ax7,'$\theta\;[\mathrm{deg}]$');

ylabel(ax7,'$\dot{\theta}\;[\mathrm{deg/s}]$');

title(ax7,sprintf('$K_p=%.4f,\\ T=%.4f\\,\\mathrm{s},\\ \\omega=%.4f\\,\\mathrm{rad/s}$',KpSelected,TSelected,omegaSelected));

finishFigure(fig7,ax7);

%% 34. FIGURE 8 : PITCH ANGLE TIME HISTORY
[fig8,ax8] = makeFigure('08 - Pitch Angle Time History',8,totalFigures,figureRows,figureColumns);

plot(ax8,t,thetaSelected,'LineWidth',1.8);

xlabel(ax8,'$t\;[\mathrm{s}]$');

ylabel(ax8,'$\theta(t)\;[\mathrm{deg}]$');

title(ax8,'$\mathrm{Pitch\ Angle\ Time\ History}$');

finishFigure(fig8,ax8);

%% 35. FIGURE 9 : PITCH RATE TIME HISTORY
[fig9,ax9] = makeFigure('09 - Pitch Rate Time History',9,totalFigures,figureRows,figureColumns);

plot(ax9,t,thetaDotSelected,'LineWidth',1.8);

xlabel(ax9,'$t\;[\mathrm{s}]$');

ylabel(ax9,'$\dot{\theta}(t)\;[\mathrm{deg/s}]$');

title(ax9,'$\mathrm{Pitch\ Rate\ Time\ History}$');

finishFigure(fig9,ax9);

%% 36. FIGURE 10 : ACTUATOR RATE
[fig10,ax10] = makeFigure('10 - Actuator Rate',10,totalFigures,figureRows,figureColumns);

hold(ax10,'on');

plot(ax10,t,x1dotSelected,'LineWidth',1.8);

yline(ax10,S,'--','$+S$','Interpreter','latex');

yline(ax10,-S,'--','$-S$','Interpreter','latex');

xlabel(ax10,'$t\;[\mathrm{s}]$');

ylabel(ax10,'$\dot{x}_1\;[\mathrm{deg/s}]$');

title(ax10,'$\mathrm{Rate\!-\!Limited\ Actuator\ Response}$');

finishFigure(fig10,ax10);

%% 37. FIGURE 11 : ACTUATOR STATE VS PITCH ANGLE
[fig11,ax11] = makeFigure('11 - Actuator State',11,totalFigures,figureRows,figureColumns);

plot(ax11,thetaSelected,x1Selected,'LineWidth',1.8);

xlabel(ax11,'$\theta\;[\mathrm{deg}]$');

ylabel(ax11,'$x_1$');

title(ax11,'$\mathrm{Limit\!-\!Cycle\ Projection:}\ x_1\ \mathrm{vs.}\ \theta$');

finishFigure(fig11,ax11);

%% 38. FIGURE 12 : ACTUATOR PHASE PORTRAIT
[fig12,ax12] = makeFigure('12 - Actuator Phase Portrait',12,totalFigures,figureRows,figureColumns);

plot(ax12,x1Selected,x1dotSelected,'LineWidth',1.8);

xlabel(ax12,'$x_1$');

ylabel(ax12,'$\dot{x}_1\;[\mathrm{deg/s}]$');

title(ax12,'$\mathrm{Actuator\ Limit\!-\!Cycle\ Phase\ Portrait}$');

finishFigure(fig12,ax12);

%% 39. FIGURE 13 : HOPF CLOSE-UP
[fig13,ax13] = makeFigure('13 - Hopf Close-Up',13,totalFigures,figureRows,figureColumns);

hold(ax13,'on');

plot(ax13,KpEP,thetaEP,'k-','LineWidth',1.5,'DisplayName','$\mathrm{Equilibrium}$');

plot(ax13,resultF.Kp,resultF.thetaMax,'LineWidth',1.8,'DisplayName','$\theta_{\max}$');

plot(ax13,resultF.Kp,resultF.thetaMin,'LineWidth',1.8,'DisplayName','$\theta_{\min}$');

if ~isempty(resultB)
    plot(ax13,resultB.Kp,resultB.thetaMax,'--','LineWidth',1.5,'HandleVisibility','off');

    plot(ax13,resultB.Kp,resultB.thetaMin,'--','LineWidth',1.5,'HandleVisibility','off');
end

plot(ax13,KpH,thetaH,'ko','MarkerFaceColor','k','MarkerSize',6,'DisplayName','$H$');

xlabel(ax13,'$K_p$');

ylabel(ax13,'$\theta\;[\mathrm{deg}]$');

title(ax13,'$\mathrm{Hopf\ Bifurcation\ Close\!-\!Up}$');

legend(ax13,'Location','best');

finishFigure(fig13,ax13);

xlim(ax13,[KpH-0.1 KpH+0.5]);

autoY(ax13);

%% 40. FIGURE 14 : LIMIT-CYCLE AMPLITUDE NEAR HOPF
[fig14,ax14] = makeFigure('14 - Hopf Amplitude',14,totalFigures,figureRows,figureColumns);

hold(ax14,'on');

plot(ax14,resultF.Kp,resultF.thetaAmp,'LineWidth',1.8,'DisplayName','$\mathrm{Forward\ LC}$');

if ~isempty(resultB)
    plot(ax14,resultB.Kp,resultB.thetaAmp,'--','LineWidth',1.5,'DisplayName','$\mathrm{Backward\ LC}$');
end

plot(ax14,KpH,0,'ko','MarkerFaceColor','k','MarkerSize',6,'DisplayName','$H$');

xlabel(ax14,'$K_p$');

ylabel(ax14,'$A_{\theta}\;[\mathrm{deg}]$');

title(ax14,'$\mathrm{Limit\!-\!Cycle\ Growth\ from\ Hopf\ Bifurcation}$');

legend(ax14,'Location','best');

finishFigure(fig14,ax14);

xlim(ax14,[KpH-0.1 KpH+0.5]);

autoY(ax14);

%% 41. DISPLAY LIMIT-CYCLE SPECIAL POINTS
fprintf('\n============================================================\n');
fprintf(' FORWARD LIMIT-CYCLE SPECIAL POINTS\n');
fprintf('============================================================\n');

for k = 1:length(sLCf)
    label = strtrim(sLCf(k).label);

    if isempty(label)
        label = '-';
    end

    fprintf('%2d   Label = %-5s   Index = %d\n',k,label,sLCf(k).index);
end

if ~isempty(sLCb)
    fprintf('\n============================================================\n');
    fprintf(' BACKWARD LIMIT-CYCLE SPECIAL POINTS\n');
    fprintf('============================================================\n');

    for k = 1:length(sLCb)
        label = strtrim(sLCb(k).label);

        if isempty(label)
            label = '-';
        end

        fprintf('%2d   Label = %-5s   Index = %d\n',k,label,sLCb(k).index);
    end
end

%% 42. FINAL RESULTS
fprintf('\n============================================================\n');
fprintf(' FINAL RESULTS\n');
fprintf('============================================================\n');

fprintf('\nHopf point:\n');
fprintf('Kp_H        = %.6f\n',KpH);
fprintf('omega_H     = %.6f rad/s\n',omegaHopf);
fprintf('T_H         = %.6f s\n',THopf);

fprintf('\nSelected nonlinear limit cycle:\n');
fprintf('Branch      = %s\n',selectedBranch);
fprintf('Kp          = %.6f\n',KpSelected);
fprintf('T           = %.6f s\n',TSelected);
fprintf('f           = %.6f Hz\n',fSelected);
fprintf('omega       = %.6f rad/s\n',omegaSelected);
fprintf('A_theta     = %.6f deg\n',thetaAmpSelected);
fprintf('A_theta_dot = %.6f deg/s\n',thetaDotAmpSelected);
fprintf('2*pi/T      = %.6f rad/s\n',2*pi/TSelected);

%% 43. SAVE RESULTS
save('X15_MatCont_Bifurcation_Results.mat',...
'xEP','vEP','sEP','hEP','fEP',...
'xLCf','vLCf','sLCf','hLCf','fLCf',...
'xLCb','vLCb','sLCb','hLCb','fLCb',...
'resultF','resultB',...
'KpH','omegaHopf','THopf',...
'KpSelected','TSelected','fSelected',...
'omegaSelected','thetaAmpSelected','thetaDotAmpSelected');

fprintf('\nResults saved as X15_MatCont_Bifurcation_Results.mat\n');
fprintf('============================================================\n');
fprintf(' ANALYSIS COMPLETED\n');
fprintf('============================================================\n');

%% LOCAL FUNCTION 1 : PROCESS LIMIT-CYCLE DATA
function result = processLCBranch(xLC,nStates,Ctheta,CthetaDot,K,S,thetac)

if isempty(xLC)
    result = [];
    return;
end

% MatCont LC matrix:
% rows 1:end-2 = periodic orbit states
% row end-1    = period
% row end      = active parameter Kp

nOrbitRows = size(xLC,1)-2;

nPoints = nOrbitRows/nStates;

if abs(nPoints-round(nPoints)) > 1e-10
    error('Unexpected MatCont limit-cycle matrix format.');
end

nPoints = round(nPoints);

nLC = size(xLC,2);

result.period = xLC(end-1,:);

result.Kp = xLC(end,:);

result.frequency = NaN(1,nLC);

result.omega = NaN(1,nLC);

validPeriod = result.period > 0 & isfinite(result.period);

result.frequency(validPeriod) = 1./result.period(validPeriod);

result.omega(validPeriod) = 2*pi./result.period(validPeriod);

result.thetaMax = NaN(1,nLC);

result.thetaMin = NaN(1,nLC);

result.thetaAmp = NaN(1,nLC);

result.thetaDotAmp = NaN(1,nLC);

result.maxX1dot = NaN(1,nLC);

result.thetaOrbit = cell(1,nLC);

result.thetaDotOrbit = cell(1,nLC);

result.x1dotOrbit = cell(1,nLC);

result.stateOrbit = cell(1,nLC);

for k = 1:nLC
    orbitVector = xLC(1:nOrbitRows,k);

    X = reshape(orbitVector,nStates,nPoints);

    x1 = X(1,:);

    theta = Ctheta*X;

    thetaDot = CthetaDot*X;

    Kp = result.Kp(k);

    x1dot = S*tanh(K*(Kp*(thetac-theta)-x1)/S);

    result.thetaMax(k) = max(theta);

    result.thetaMin(k) = min(theta);

    result.thetaAmp(k) = (max(theta)-min(theta))/2;

    result.thetaDotAmp(k) = (max(thetaDot)-min(thetaDot))/2;

    result.maxX1dot(k) = max(abs(x1dot));

    result.thetaOrbit{k} = theta;

    result.thetaDotOrbit{k} = thetaDot;

    result.x1dotOrbit{k} = x1dot;

    result.stateOrbit{k} = X;
end

end

%% LOCAL FUNCTION 2 : FIND POINT CLOSEST TO TARGET
function [index,minError] = closestPoint(data,target)

valid = isfinite(data);

if ~any(valid)
    index = [];
    minError = inf;
    return;
end

validIndex = find(valid);

errorData = abs(data(valid)-target);

[minError,localIndex] = min(errorData);

index = validIndex(localIndex);

end

%% LOCAL FUNCTION 3 : CREATE AND POSITION INTERACTIVE FIGURE
function [fig,ax] = makeFigure(figName,figNumber,totalFigures,nRows,nColumns)

screenSize = get(groot,'ScreenSize');

screenX = screenSize(1);

screenY = screenSize(2);

screenWidth = screenSize(3);

screenHeight = screenSize(4);

gap = 6;

leftMargin = 5;

rightMargin = 5;

bottomMargin = 45;

topMargin = 30;

usableWidth = screenWidth-leftMargin-rightMargin;

usableHeight = screenHeight-bottomMargin-topMargin;

figWidth = floor((usableWidth-(nColumns-1)*gap)/nColumns);

figHeight = floor((usableHeight-(nRows-1)*gap)/nRows);

row = ceil(figNumber/nColumns);

columnInRow = mod(figNumber-1,nColumns)+1;

figuresBeforeRow = (row-1)*nColumns;

figuresInRow = min(nColumns,totalFigures-figuresBeforeRow);

rowWidth = figuresInRow*figWidth+(figuresInRow-1)*gap;

startX = screenX+leftMargin+(usableWidth-rowWidth)/2;

xPosition = startX+(columnInRow-1)*(figWidth+gap);

yPosition = screenY+bottomMargin+(nRows-row)*(figHeight+gap);

fig = figure('Name',figName,...
'NumberTitle','off',...
'Color','w',...
'MenuBar','none',...
'ToolBar','none',...
'Position',[xPosition yPosition figWidth figHeight]);

enableLegacyExplorationModes(fig);

ax = axes('Parent',fig);

hold(ax,'on');

grid(ax,'on');

box(ax,'on');

ax.FontSize = 8;

ax.LineWidth = 0.8;

ax.TickLabelInterpreter = 'latex';

try
    tb = axtoolbar(ax,{'zoomin','zoomout','pan','restoreview','datacursor'});

    tb.Visible = 'on';
catch
    set(fig,'ToolBar','figure');
end

end

%% LOCAL FUNCTION 4 : FINAL FIGURE SETTINGS
function finishFigure(fig,ax)

grid(ax,'on');

box(ax,'on');

axis(ax,'tight');

xLim = ax.XLim;

xRange = diff(xLim);

if xRange > 0
    ax.XLim = xLim + [-0.025 0.025]*xRange;
end

yLim = ax.YLim;

yRange = diff(yLim);

if yRange > 0
    ax.YLim = yLim + [-0.05 0.05]*yRange;
else
    yCenter = mean(yLim);

    margin = 0.05*max(1,abs(yCenter));

    ax.YLim = [yCenter-margin yCenter+margin];
end

z = zoom(fig);

z.ActionPostCallback = @(~,eventData) autoY(eventData.Axes);

p = pan(fig);

p.ActionPostCallback = @(~,eventData) autoY(eventData.Axes);

zoom(fig,'off');

pan(fig,'off');

ax.FontSize = 8;

ax.XLabel.FontSize = 9;

ax.YLabel.FontSize = 9;

ax.Title.FontSize = 10;

if ~isempty(ax.Legend)
    ax.Legend.FontSize = 7;
end

end

%% LOCAL FUNCTION 5 : AUTOMATIC Y-AXIS FOCUS
function autoY(ax)

if isempty(ax) || ~isvalid(ax)
    return;
end

currentX = ax.XLim;

lines = findobj(ax,'Type','line');

visibleY = [];

for k = 1:length(lines)
    X = lines(k).XData;

    Y = lines(k).YData;

    if isempty(X) || isempty(Y)
        continue;
    end

    X = X(:);

    Y = Y(:);

    if length(X) ~= length(Y)
        continue;
    end

    inside = X >= currentX(1) & X <= currentX(2) & isfinite(X) & isfinite(Y);

    if any(inside)
        visibleY = [visibleY;Y(inside)];
    end
end

if isempty(visibleY)
    return;
end

yMin = min(visibleY);

yMax = max(visibleY);

if abs(yMax-yMin) < 1e-12
    margin = 0.05*max(1,abs(yMax));
else
    margin = 0.07*(yMax-yMin);
end

newYLim = [yMin-margin yMax+margin];

if all(isfinite(newYLim)) && newYLim(1) < newYLim(2)
    ax.YLim = newYLim;
end

end
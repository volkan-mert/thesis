
% Initialize MatCont with "init" and add the paths.

%%
clear;
clearvars -global cds
clearvars -global hds

% Discretization parameters
M=40; % degree of collocation polynomial
handles=str2func(['@MackeyGlass',num2str(M)]); % available for M=5,10,20,30,40
opt=contset;

%%
fprintf('Starting continuation of the nontrivial equilibrium\n');
beta=2;
gamma=1;
n=6;
tau=0.5;
par=[beta,gamma,n,tau]';

x=ones(M+1,1); % initial equilibrium estimate
ap=4; % active parameter
TOL=1e-8; 
opt=contset(opt,'FunTolerance',TOL);
opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',0);
[x0,v0]=init_EP_EP(handles,x,par,ap);
[x,v,s,h,f]=cont(@equilibrium,x0,v0,opt);
% check if the continuation has detected Hopf, otherwise compute backward
% branch
if strcmp(s(2).label,'H ')==0
    opt=contset(opt,'Backward',1);
    [x,v,s,h,f]=cont(@equilibrium,x0,v0,opt);
end
figure(1); clf
cpl(x,v,s,[size(x,1) 1]);
xlabel('$\tau$','Interpreter','latex')

xH=x; vH=v; sH=s; hH=h; fH=f; 

if strcmp(s(2).label,'H ')==1
    fprintf('Hopf bifurcation at tau = %.15f\n', x(end,s(2).index));
    if strcmp(s(3).label,'H ')==1
        fprintf('Second Hopf bifurcation at tau = %.15f\n', x(end,s(3).index));
    end
else
    return;
end

%% 
fprintf('Starting continuation of the limit cycle branch from Hopf\n');
par(ap)=xH(end,sH(2).index);
ntst=40; % number of intervals
ncol=4; % degree of polynomial
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'MaxStepsize',0.5);
opt=contset(opt,'Backward',0);
opt=contset(opt,'IgnoreSingularity',[1,3,4]); % ignore what is not PD: [BPC,PD,LPC,NS]
opt=contset(opt,'TSearchOrder',1); % 1: increasing; other: decreasing
[x0,v0]=init_H_LC(handles,xH(1:M+1,sH(2).index),par,ap,1e-6,ntst,ncol);
[x,v,s,h,f]= cont(@limitcycle,x0,v0,opt);
iPD=findPD(s); % auxiliary function to find the index corresponding to PD

if isempty(iPD)
    opt=contset(opt,'Backward',1);
    [x,v,s,h,f]= cont(@limitcycle,x0,v0,opt);
    iPD=findPD(s);
end

if isempty(iPD)
    return
else
    fprintf('PD bifurcation at tau = %.15f\n', x(end,s(iPD).index));
end

cpl(x,v,s,[size(x,1) 1]);
xPD=x; vPD=v; sPD=s; hPD=h; fPD=f; 

%% 
fprintf('Starting continuation of the limit cycle branch from PD\n');
[x0,v0]=init_PD_LC(handles,xPD,sPD(iPD),ntst,ncol,1e-6);
[x,v,s,h,f]= cont(@limitcycle,x0,v0,opt);
iPD2=findPD(s);
if ~isempty(iPD2)
    fprintf('PD bifurcation at tau = %.15f\n', x(end,s(iPD2).index));
end
cpl(x,v,s,[size(x,1) 1]);

% %% 
% fprintf('Starting continuation of Hopf in two parameters\n');
% ap2=3;
% par(ap)=xH(end,sH(2).index);
% [x0,v0]=init_H_H(handles,xH(1:M+1,sH(2).index),par,[ap ap2]);
% opt=contset(opt,'Backward',0);
% [x,v,s,h,f]= cont(@hopf,x0,v0,opt);
% figure(2); clf
% cpl(x,v,s,[size(x,1)-2 size(x,1)-1]); hold on
% xlabel('$\tau$','Interpreter','latex')
% ylabel('$n$','Interpreter','latex')
% axis([0 10 0 10]) 
% 
% opt=contset(opt,'Backward',1);
% [x,v,s,h,f]= cont(@hopf,x0,v0,opt);
% cpl(x,v,s,[size(x,1)-2 size(x,1)-1]);
% 
% %% 
% fprintf('Starting continuation of PD in two parameters\n');
% ap2=3;
% opt=contset(opt,'MaxStepsize',20);
% opt=contset(opt,'IgnoreSingularity',1);
% opt=contset(opt,'Backward',0);
% [x0,v0]=init_PD_PD(handles,xPD,sPD(iPD),[ap ap2],ntst,ncol);
% [x,v,s,h,f]= cont(@perioddoubling,x0,v0,opt);
% cpl(x,v,s,[size(x,1)-1 size(x,1)]);
% 
% opt=contset(opt,'Backward',1);
% [x,v,s,h,f]= cont(@perioddoubling,x0,v0,opt);
% cpl(x,v,s,[size(x,1)-1 size(x,1)]);

%% Auxiliary functions

function i=findPD(s)
% Finds the index corresponding to the label 'PD ' in the vector s
    i=[];
    for is=2:length(s)
        if strcmp(s(is).label,'PD ')==1
            i=is;
            break
        end
    end
end
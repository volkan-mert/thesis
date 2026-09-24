function out = X15PIORelay4states
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_R,par_A,par_omega,par_tau,par_Kp,par_Tr,par_mu)
dydt=[kmrgd(2);
(par_R*tanh((par_Kp/par_R)*(kmrgd(3)*cos(par_omega*par_tau)-kmrgd(4)*sin(par_omega*par_tau)-kmrgd(1)))-kmrgd(2))/par_Tr;
par_omega*kmrgd(4)+par_mu*(1-(kmrgd(3)^2+kmrgd(4)^2)/(par_A^2))*kmrgd(3);
-par_omega*kmrgd(3)+par_mu*(1-(kmrgd(3)^2+kmrgd(4)^2)/(par_A^2))*kmrgd(4);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(X15PIORelay4states);
y0=[0,0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_R,par_A,par_omega,par_tau,par_Kp,par_Tr,par_mu)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_R,par_A,par_omega,par_tau,par_Kp,par_Tr,par_mu)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_R,par_A,par_omega,par_tau,par_Kp,par_Tr,par_mu)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_R,par_A,par_omega,par_tau,par_Kp,par_Tr,par_mu)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_R,par_A,par_omega,par_tau,par_Kp,par_Tr,par_mu)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_R,par_A,par_omega,par_tau,par_Kp,par_Tr,par_mu)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_R,par_A,par_omega,par_tau,par_Kp,par_Tr,par_mu)

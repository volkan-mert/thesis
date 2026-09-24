function out = RLEAshkenas
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
function dydt = fun_eval(t,kmrgd,par_K,par_b,par_A,par_omega,par_lambda)
dydt=[par_b*tanh(par_K*(kmrgd(2)-kmrgd(1))/par_b);
par_omega*kmrgd(3)+par_lambda*(1-(kmrgd(2)^2+kmrgd(3)^2)/(par_A^2))*kmrgd(2);
-par_omega*kmrgd(2)+par_lambda*(1-(kmrgd(2)^2+kmrgd(3)^2)/(par_A^2))*kmrgd(3);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(RLEAshkenas);
y0=[0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_K,par_b,par_A,par_omega,par_lambda)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_K,par_b,par_A,par_omega,par_lambda)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_K,par_b,par_A,par_omega,par_lambda)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_K,par_b,par_A,par_omega,par_lambda)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_K,par_b,par_A,par_omega,par_lambda)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_K,par_b,par_A,par_omega,par_lambda)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_K,par_b,par_A,par_omega,par_lambda)

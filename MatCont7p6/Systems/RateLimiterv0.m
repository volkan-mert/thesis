function out = RateLimiterv0
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
function dydt = fun_eval(t,kmrgd,par_K,par_S,par_yc)
dydt=[par_S*tanh(par_K*(par_yc-kmrgd(1))/par_S);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(RateLimiterv0);
y0=[0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_K,par_S,par_yc)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_K,par_S,par_yc)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_K,par_S,par_yc)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_K,par_S,par_yc)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_K,par_S,par_yc)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_K,par_S,par_yc)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_K,par_S,par_yc)

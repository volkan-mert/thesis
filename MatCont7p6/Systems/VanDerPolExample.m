function out = VanDerPolExample
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
function dydt = fun_eval(t,kmrgd,par_mu)
dydt=[kmrgd(2);
par_mu*(1-kmrgd(1)*kmrgd(1))*kmrgd(2)-kmrgd(1);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(VanDerPolExample);
y0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_mu)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_mu)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_mu)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_mu)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_mu)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_mu)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_mu)

function out = Lab3ProbCPredatorDoublePreySys
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
function dydt = fun_eval(t,kmrgd,par_beta)
dydt=[kmrgd(1)*(2.4-kmrgd(1)-6*kmrgd(2)-4*kmrgd(3));
kmrgd(2)*(par_beta-kmrgd(1)-kmrgd(2)-10*kmrgd(3));
-kmrgd(3)*(1-0.25*kmrgd(1)-4*kmrgd(2)+kmrgd(3));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Lab3ProbCPredatorDoublePreySys);
y0=[0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_beta)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_beta)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_beta)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_beta)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_beta)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_beta)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_beta)

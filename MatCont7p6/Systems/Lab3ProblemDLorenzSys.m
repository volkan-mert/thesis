function out = Lab3ProblemDLorenzSys
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
function dydt = fun_eval(t,kmrgd,par_sigma,par_b,par_r)
dydt=[-kmrgd(1)*par_sigma+kmrgd(2)*par_sigma;
par_r*kmrgd(1)-kmrgd(2)-kmrgd(1)*kmrgd(3);
-par_b*kmrgd(3)+kmrgd(1)*kmrgd(2);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Lab3ProblemDLorenzSys);
y0=[0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_sigma,par_b,par_r)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_sigma,par_b,par_r)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_sigma,par_b,par_r)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_sigma,par_b,par_r)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_sigma,par_b,par_r)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_sigma,par_b,par_r)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_sigma,par_b,par_r)

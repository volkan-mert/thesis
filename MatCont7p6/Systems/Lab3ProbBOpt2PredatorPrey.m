function out = Lab3ProbBOpt2PredatorPrey
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
function dydt = fun_eval(t,kmrgd,par_alpha,par_delta)
dydt=[1-(exp(kmrgd(2)))/(1+par_alpha*exp(kmrgd(1)));
-1+(exp(kmrgd(1)))/(1+par_alpha*exp(kmrgd(1)))-par_delta*exp(kmrgd(2));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Lab3ProbBOpt2PredatorPrey);
y0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_alpha,par_delta)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_alpha,par_delta)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_alpha,par_delta)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_alpha,par_delta)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_alpha,par_delta)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_alpha,par_delta)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_alpha,par_delta)

function out = RoesslerSolExLab3
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
function dydt = fun_eval(t,kmrgd,par_A,par_B,par_C)
dydt=[-kmrgd(2)-kmrgd(3);
kmrgd(1)+par_A*kmrgd(2);
par_B*kmrgd(1)-par_C*kmrgd(3)+kmrgd(1)*kmrgd(3);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(RoesslerSolExLab3);
y0=[0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_A,par_B,par_C)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_A,par_B,par_C)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_A,par_B,par_C)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_A,par_B,par_C)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_A,par_B,par_C)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_A,par_B,par_C)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_A,par_B,par_C)

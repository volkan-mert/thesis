function out = X15PIOSoftGlide
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
function dydt = fun_eval(t,kmrgd,par_Kp,par_K,par_S,par_R,par_thetac)
dydt=[max(par_R,min(par_S,par_K*(par_Kp*(par_thetac-(6.02372*kmrgd(2)+7.346*kmrgd(3)))-kmrgd(1))));
kmrgd(3);
kmrgd(4);
kmrgd(1)-5.29*kmrgd(3)-1.42*kmrgd(4);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(X15PIOSoftGlide);
y0=[0,0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_Kp,par_K,par_S,par_R,par_thetac)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_Kp,par_K,par_S,par_R,par_thetac)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_Kp,par_K,par_S,par_R,par_thetac)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_Kp,par_K,par_S,par_R,par_thetac)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_Kp,par_K,par_S,par_R,par_thetac)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_Kp,par_K,par_S,par_R,par_thetac)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_Kp,par_K,par_S,par_R,par_thetac)

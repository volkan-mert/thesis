function varargout = auto_matlab(mode, varargin)
%AUTO_MATLAB  Single-file MATLAB version of AUTO's core capability.
%
%   Pure MATLAB / Octave. No toolboxes, no external files. A pseudo-arclength
%   predictor-corrector continuation of equilibria of du/dt = f(u,p), with
%   automatic detection of Limit Points (folds, "LP") and Hopf bifurcations
%   ("HB") -- the equilibrium-continuation core of AUTO (Doedel et al.).
%
%   HOW TO USE
%   ----------
%   Run the built-in demos (produces two bifurcation diagrams):
%       >> auto_matlab                 % runs both demos
%       >> auto_matlab fold            % S-curve, detects LP at lam = +/-0.3849
%       >> auto_matlab hopf            % Brusselator, detects HB at B = 2
%
%   Use the engine on YOUR OWN model:
%       f    = @(u,p) [ ... ];         % RHS, n x 1 column vector
%       jac  = @(u,p) [ ... ];         % state Jacobian df/du (optional)
%       opts = struct('jac',jac,'ds',0.02,'p_min',0,'p_max',10);
%       out  = auto_matlab('continue', f, u0, p0, opts);
%       %   out.U (n x N), out.p (1 x N), out.folds, out.hopfs, out.tangent
%
%   METHOD (standard, textbook)
%     Extended system  H(x)=[ f(u,p) ; t_prev'*(x-x_pred) ]=0 , x=[u;p].
%     Tangent = unit null vector of the (n x n+1) extended Jacobian.
%     Fold test  : sign change of t(n+1)  (the dp/ds component).
%     Hopf test  : sign change of det( bialternate-product of df/du ),
%                  confirmed by a genuine purely-imaginary eigenvalue pair.
%
%   Refs: Keller (1977); Seydel, Practical Bifurcation & Stability Analysis;
%   Kuznetsov, Elements of Applied Bifurcation Theory; Doedel et al., AUTO.

if nargin == 0, mode = 'all'; end

switch lower(mode)
    case 'all'
        demo_fold(); demo_hopf();
    case 'fold'
        out = demo_fold();  if nargout, varargout{1} = out; end
    case 'hopf'
        out = demo_hopf();  if nargout, varargout{1} = out; end
    case 'continue'
        % auto_matlab('continue', f, u0, p0, opts)
        f = varargin{1}; u0 = varargin{2}; p0 = varargin{3};
        if numel(varargin) >= 4, opts = varargin{4}; else, opts = struct(); end
        varargout{1} = pac_continue(f, u0, p0, opts);
    otherwise
        error('auto_matlab:mode','Unknown mode "%s". Use fold|hopf|continue.',mode);
end
end

% ========================================================================
% ============================  DEMOS  ===================================
% ========================================================================
function out = demo_fold()
% Limit-point (fold) demo: lam = u1 - u1^3 gives an S-curve.
% Analytic folds: u1 = +/-1/sqrt(3), lam = +/-2/(3*sqrt(3)) = +/-0.3849.
f   = @(u,p) [ -u(1)^3 + u(1) - p ;  -u(2) ];
jac = @(u,p) [ -3*u(1)^2 + 1 , 0 ; 0 , -1 ];
opts = struct('jac',jac,'ds',0.01,'ds_max',0.05, ...
              'n_max',4000,'p_min',-0.5,'p_max',0.5,'direction',+1);
% start on the middle branch (u1=0, lam=0) so both folds are traversed
out = pac_continue(f, [0;0], 0, opts);

figure('Color','w'); hold on; box on;
plot(out.p, out.U(1,:), 'b-', 'LineWidth', 1.6);
for i = 1:numel(out.folds)
    plot(out.folds(i).p, out.folds(i).u(1), 'rs','MarkerSize',9,'MarkerFaceColor','r');
    text(out.folds(i).p, out.folds(i).u(1), sprintf('  LP \\lambda=%.4f',out.folds(i).p),'Color','r');
end
xlabel('\lambda (parameter)'); ylabel('u_1');
title('Fold / limit-point continuation (AUTO-style S-curve)'); grid on;
fprintf('Analytic folds: lam = +/-%.4f at u1 = +/-%.4f\n', 2/(3*sqrt(3)), 1/sqrt(3));
end

function out = demo_hopf()
% Hopf demo: Brusselator continued in B; Hopf at B = 1+A^2 (A=1 -> B=2).
A = 1.0;
f   = @(u,B) [ A - (B+1)*u(1) + u(1)^2*u(2) ;  B*u(1) - u(1)^2*u(2) ];
jac = @(u,B) [ -(B+1) + 2*u(1)*u(2) ,  u(1)^2 ; B - 2*u(1)*u(2) , -u(1)^2 ];
B0 = 0.5;
opts = struct('jac',jac,'ds',0.02,'ds_max',0.05, ...
              'n_max',600,'p_min',0.4,'p_max',4.0,'direction',+1);
out = pac_continue(f, [A; B0/A], B0, opts);

realmax = zeros(1,numel(out.p));
for k = 1:numel(out.p)
    realmax(k) = max(real(eig(jac(out.U(:,k), out.p(k)))));
end
figure('Color','w'); hold on; box on;
plot(out.p, realmax, 'b-', 'LineWidth',1.6);
plot(xlim,[0 0],'k--');
for i = 1:numel(out.hopfs)
    plot(out.hopfs(i).p, 0, 'mo','MarkerSize',10,'MarkerFaceColor','m');
    text(out.hopfs(i).p, 0, sprintf('  HB B=%.4f, \\omega=%.3f',out.hopfs(i).p,out.hopfs(i).omega),'Color','m');
end
xlabel('B (parameter)'); ylabel('max Re(\lambda) of df/du');
title('Hopf detection on the Brusselator equilibrium branch'); grid on;
fprintf('Analytic Hopf: B = 1+A^2 = %.4f, omega = A = %.4f rad/s\n', 1+A^2, A);
end

% ========================================================================
% ==========================  ENGINE  ====================================
% ========================================================================
function out = pac_continue(f, u0, p0, opts)
if nargin < 4, opts = struct(); end
gp = @(s,d) getfielddef(opts,s,d);
ds        = gp('ds',        0.02);
ds_min    = gp('ds_min',    1e-5);
ds_max    = gp('ds_max',    0.20);
n_max     = gp('n_max',     1000);
p_min     = gp('p_min',     -Inf);
p_max     = gp('p_max',     +Inf);
direction = gp('direction', +1);
tol       = gp('tol',       1e-9);
newt_max  = gp('newt_max',  25);
det_fold  = gp('detect_fold', true);
det_hopf  = gp('detect_hopf', true);
verbose   = gp('verbose',   true);
has_jac   = isfield(opts,'jac');

u0 = u0(:);  n = numel(u0);
if has_jac, dfdu = @(u,p) opts.jac(u,p);
else,       dfdu = @(u,p) fd_state_jac(f,u,p,n); end
gext = @(u,p) [dfdu(u,p), fd_dfdp(f,u,p)];

% converge starting equilibrium (Newton in u at fixed p0)
u = u0; p = p0;
for it = 1:50
    r = f(u,p);
    if norm(r) < tol, break; end
    u = u - dfdu(u,p) \ r;
end
if norm(f(u,p)) > 1e-6
    error('pac_continue:startpoint','Could not converge an equilibrium at p0=%g.',p0);
end

% initial tangent
t = null_tangent(gext(u,p), [], n);
if sign(t(n+1)) ~= sign(direction), t = -t; end

U = u;  P = p;  T = t;
folds = struct('u',{},'p',{});
hopfs = struct('u',{},'p',{},'omega',{});
tauLP_prev = t(n+1);
tauHB_prev = hopf_test(dfdu(u,p));
x = [u;p];  stop = 'reached n_max';
if verbose, fprintf('  step      p            |u|        type\n'); end

for k = 1:n_max
    xpred = x + ds*t;
    [xc, converged, iters] = corrector(f,gext,xpred,t,n,tol,newt_max);
    if ~converged
        ds = ds/2;
        if ds < ds_min, stop = 'corrector failed (ds<ds_min)'; break; end
        continue;
    end
    uc = xc(1:n); pc = xc(n+1);
    tnew = null_tangent(gext(uc,pc), t, n);

    tauLP = tnew(n+1);
    if det_fold && tauLP_prev*tauLP < 0
        [uf,pf] = refine_fold(f,gext,x,xc,t,tnew,n,tol,newt_max);
        folds(end+1) = struct('u',uf,'p',pf); 
        if verbose, fprintf('  %4d  %11.5f  %9.4f   LP (fold)\n',k,pf,norm(uf)); end
    end
    tauLP_prev = tauLP;

    tauHB = hopf_test(dfdu(uc,pc));
    if det_hopf && tauHB_prev*tauHB < 0
        [uh,ph,om,isHopf] = refine_hopf(f,dfdu,gext,x,xc,t,tnew,n,tol,newt_max);
        if isHopf
            hopfs(end+1) = struct('u',uh,'p',ph,'omega',om); 
            if verbose, fprintf('  %4d  %11.5f  %9.4f   HB (Hopf, w=%.4f)\n',k,ph,norm(uh),om); end
        end
    end
    tauHB_prev = tauHB;

    x = xc; t = tnew;
    U = [U, x(1:n)]; P = [P, x(n+1)]; T = [T, t]; 

    if iters <= 3, ds = min(ds*1.3, ds_max);
    elseif iters >= 8, ds = max(ds/1.3, ds_min); end

    if x(n+1) < p_min, stop = 'p < p_min'; break; end
    if x(n+1) > p_max, stop = 'p > p_max'; break; end
end

out = struct('U',U,'p',P,'tangent',T,'folds',folds,'hopfs',hopfs,'stop',stop);
if verbose
    fprintf('  done: %d points, %d fold(s), %d Hopf(s). Stop: %s\n', ...
            size(U,2), numel(folds), numel(hopfs), stop);
end
end

% ------------------------- helper functions ------------------------------
function v = getfielddef(s,name,def)
if isfield(s,name), v = s.(name); else, v = def; end
end

function J = fd_state_jac(f,u,p,n)
J = zeros(n,n); f0 = f(u,p);
for j = 1:n
    h = 1e-7*max(1,abs(u(j)));
    up = u; up(j) = up(j)+h;
    J(:,j) = (f(up,p)-f0)/h;
end
end

function d = fd_dfdp(f,u,p)
h = 1e-7*max(1,abs(p));
d = (f(u,p+h)-f(u,p-h))/(2*h);
end

function t = null_tangent(G, tprev, n)
if isempty(tprev), row = [zeros(1,n) 1]; else, row = tprev(:)'; end
A = [G; row];  b = [zeros(n,1); 1];
t = A\b;  t = t/norm(t);
if ~isempty(tprev) && dot(t,tprev) < 0, t = -t; end
end

function [xc,ok,it] = corrector(f,gext,xpred,t,n,tol,newt_max)
xc = xpred; ok = false; it = 0;
for it = 1:newt_max
    u = xc(1:n); p = xc(n+1);
    H = [f(u,p); dot(t, xc-xpred)];
    if norm(H) < tol, ok = true; return; end
    J = [gext(u,p); t(:)'];
    xc = xc - J\H;
end
if norm([f(xc(1:n),xc(n+1)); dot(t,xc-xpred)]) < 1e-7, ok = true; end
end

function [uf,pf] = refine_fold(f,gext,xa,xb,ta,tb,n,tol,newt_max)
fa = ta(n+1); fb = tb(n+1);
xf = xa + (xb-xa)*(fa/(fa-fb));
tmid = null_tangent(gext(xf(1:n),xf(n+1)), tb, n);
[xf,~] = corrector(f,gext,xf,tmid,n,tol,newt_max);
uf = xf(1:n); pf = xf(n+1);
end

function [uh,ph,omega,isHopf] = refine_hopf(f,dfdu,gext,xa,xb,ta,tb,n,tol,newt_max)
fa = hopf_test(dfdu(xa(1:n),xa(n+1)));
fb = hopf_test(dfdu(xb(1:n),xb(n+1)));
xh = xa + (xb-xa)*(fa/(fa-fb));
tmid = null_tangent(gext(xh(1:n),xh(n+1)), tb, n);
[xh,~] = corrector(f,gext,xh,tmid,n,tol,newt_max);
uh = xh(1:n); ph = xh(n+1);
ev = eig(dfdu(uh,ph));
[~,idx] = sort(abs(real(ev)));
lead = ev(idx(1));  omega = abs(imag(lead));
isHopf = (omega > 1e-4) && (abs(real(lead)) < 1e-2);
end

function tau = hopf_test(A)
n = size(A,1);
if n < 2, tau = 1; return; end
tau = det(bialt(A));
end

function B = bialt(A)
n = size(A,1);
idx = [];
for p = 2:n
    for q = 1:p-1, idx = [idx; p q]; end 
end
m = size(idx,1);  B = zeros(m,m);
for i = 1:m
    p = idx(i,1); q = idx(i,2);
    for j = 1:m
        r = idx(j,1); s = idx(j,2);
        if r == q,                    val = -A(p,s);
        elseif (r ~= p) && (s == q),  val =  A(p,r);
        elseif (r == p) && (s == q),  val =  A(p,p) + A(q,q);
        elseif (r == p) && (s ~= q),  val =  A(q,s);
        elseif s == p,                val = -A(q,r);
        else,                         val = 0;
        end
        B(i,j) = val;
    end
end
end

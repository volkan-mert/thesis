function [out, delay] = MackeyGlass5
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
delay{1} = @make_kmrgd_from_func;

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_beta,par_gamma,par_n,par_tau)
M = 5; % degree of collocation polynomials

dim = 1; % total dimension of the delay equation
dim_d = 1; % dimension of the differential part
dim_r = 0; % dimension of the renewal part
ind_d = [1]; % indices of differential equations
%ind_r = []; % indices of renewal equations
% indices of all components of differential part in kmrgd:
ind_d_all = [1,2,3,4,5,6];
% indices of all components of renewal part in kmrgd:
%ind_r_all = [];

delays = -[-par_tau];
if any(delays < 0)
    % the equation is advanced instead of delayed: we need to stop
    dydt = NaN(size(kmrgd));
    return
else
    tau_max = max(delays);
    if tau_max == 0
        % the equation is not delayed: we need to stop
        dydt = NaN(size(kmrgd));
        return
    end
end

% Chebyshev type II nodes (extrema) in [-1, 1]:
unit_nodes = [1;0.809016994374947;0.309016994374947;-0.309016994374947;-0.809016994374947;-1];
nodes = (unit_nodes-1)*tau_max/2; % Chebyshev type II nodes (extrema) in the delay interval [-tau_max, 0]
% weights for barycentric interpolation on Chebyshev type II nodes (extrema) (scaling-independent):
weights = [0.5;-1;1;-1;1;-0.5];
% differentiation matrix for Chebyshev type II nodes (extrema) in [-1, 1]:
unit_DM = [
8.5,-10.4721359549996,2.89442719099992,-1.52786404500042,1.10557280900008,-0.5;
2.6180339887499,-1.17082039324994,-2,0.894427190999916,-0.618033988749895,0.276393202250021;
-0.723606797749979,2,-0.170820393249937,-1.61803398874989,0.894427190999916,-0.381966011250105;
0.381966011250105,-0.894427190999916,1.61803398874989,0.170820393249937,-2,0.723606797749979;
-0.276393202250021,0.618033988749895,-0.894427190999916,2,1.17082039324994,-2.61803398874989;
0.5,-1.10557280900008,1.52786404500042,-2.89442719099992,10.4721359549996,-8.49999999999999;
];
DM = unit_DM*2/tau_max; % differentiation matrix for Chebyshev type II nodes (extrema) in the delay interval [-tau_max, 0]

% kmrgd has length dim_d+M*dim
% kmrgd(i), i in 1:dim_d = value of the i-th differential component of the state
%     (i.e. the ind_d(i)-th component of the full state) at the current time
% kmrgd(dim_d+(m-1)*dim+i), m in 1:M, i in 1:dim_d = value of the i-th component
%     of the (full) state at the m-th node in the past
% kmrgd([i,(dim_d+ind_d(i)):dim:(dim_d+M*dim))]), i in 1:dim_d = values of the
%     i-th differential component of the state (i.e. the ind_d(i)-th component of
%     the full state) at all nodes (length M+1); equivalent to kmrgd(ind_d_all(i:dim_d:end))
% kmrgd(dim_d+(m-1)*dim+(1:dim)), m in 1:M = values of the (full) state at the m-th node in the past
% kmrgd(ind_d_all) = values of the differential components of the state at all the nodes
%     (M+1 blocks of length dim_d)
% kmrgd(ind_r_all) = values of the (integrated) renewal components of the state at
%     the nodes in the past (value at current time is always 0) (M blocks of length dim_r)

% values of the derivative of the renewal part of the (integrated) state at all the nodes
%     (M+1 blocks of length dim_r):
%kmrgd_r_der = reshape(transpose(DM(:,2:end)*transpose(reshape(kmrgd,dim_r,M))),(M+1)*dim_r,1); % only renewal equations
%kmrgd_r_der = reshape(transpose(DM(:,2:end)*transpose(reshape(kmrgd(ind_r_all),dim_r,M))),(M+1)*dim_r,1); % coupled differential and renewal equations
% equivalent to kron(DM(:,2:end),eye(dim_r))*kmrgd(ind_r_all)

% kmrgd_r_der(i), i in 1:dim_r = value of the derivative of the i-th (integrated)
%     renewal component of the state (i.e. the ind_r(i)-th component of the full state)
%     at the current time (i.e. true value of the i-th component)
% kmrgd_r_der(m*dim_r+i), m in 1:M, i in 1:dim_r = value of the same at m-th node in the past
% kmrgd_r_der(i:dim_r:(M+1)*dim_r), i in 1:dim_r = values of the same at all nodes (length M+1)

% renewal part of the right-hand side:
%F_r = [
%];
% F_r extended with 0 for the differential part:
% (ind_r_all-dim_d removes the offset in the indices due to the first block in kmrgd)
%F_r0 = repmat(F_r,M,1); % only renewal equations
%F_r0 = zeros(M*dim,1); F_r0(ind_r_all-dim_d) = repmat(F_r,M,1); % coupled differential and renewal equations

% kmrgd with added zeros for the renewal components at current time (M+1 blocks of length dim):
kmrgd0 = kmrgd; % only differential equations
%kmrgd0 = [zeros(dim_r,1); kmrgd]; % only renewal equations
%kmrgd0 = zeros((M+1)*dim,1); kmrgd0([ind_d,dim+1:end]) = kmrgd; % coupled differential and renewal equations

dydt=[
par_beta*kmrgd(dim_d+(M-1)*dim+1)/(1+kmrgd(dim_d+(M-1)*dim+1)^par_n)-par_gamma*kmrgd(1);
reshape(transpose(DM(2:end,:)*transpose(reshape(kmrgd0,dim,M+1))),M*dim,1); % only differential equations
%reshape(transpose(DM(2:end,:)*transpose(reshape(kmrgd0,dim,M+1))),M*dim,1)-F_r0; % only/also renewal equations
];
% in the last line, reshape(...) is equivalent to kron(DM(2:end,:),eye(dim))*kmrgd0

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(MackeyGlass5);
y0=[0,0,0,0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_beta,par_gamma,par_n,par_tau)
M = 5; % degree of collocation polynomials

dim = 1; % total dimension of the delay equation
dim_d = 1; % dimension of the differential part
dim_r = 0; % dimension of the renewal part
ind_d = [1]; % indices of differential equations
%ind_r = []; % indices of renewal equations
% indices of all components of differential part in kmrgd:
ind_d_all = [1,2,3,4,5,6];
% indices of all components of renewal part in kmrgd:
%ind_r_all = [];

delays = -[-par_tau];
if any(delays < 0)
    % the equation is advanced instead of delayed: we need to stop
    jac = NaN(length(kmrgd));
    return
else
    tau_max = max(delays);
    if tau_max == 0
        % the equation is not delayed: we need to stop
        jac = NaN(length(kmrgd));
        return
    end
end

% Chebyshev type II nodes (extrema) in [-1, 1]:
unit_nodes = [1;0.809016994374947;0.309016994374947;-0.309016994374947;-0.809016994374947;-1];
nodes = (unit_nodes-1)*tau_max/2; % Chebyshev type II nodes (extrema) in the delay interval [-tau_max, 0]
% weights for barycentric interpolation on Chebyshev type II nodes (extrema) (scaling-independent):
weights = [0.5;-1;1;-1;1;-0.5];
% differentiation matrix for Chebyshev type II nodes (extrema) in [-1, 1]:
unit_DM = [
8.5,-10.4721359549996,2.89442719099992,-1.52786404500042,1.10557280900008,-0.5;
2.6180339887499,-1.17082039324994,-2,0.894427190999916,-0.618033988749895,0.276393202250021;
-0.723606797749979,2,-0.170820393249937,-1.61803398874989,0.894427190999916,-0.381966011250105;
0.381966011250105,-0.894427190999916,1.61803398874989,0.170820393249937,-2,0.723606797749979;
-0.276393202250021,0.618033988749895,-0.894427190999916,2,1.17082039324994,-2.61803398874989;
0.5,-1.10557280900008,1.52786404500042,-2.89442719099992,10.4721359549996,-8.49999999999999;
];
DM = unit_DM*2/tau_max; % differentiation matrix for Chebyshev type II nodes (extrema) in the delay interval [-tau_max, 0]

% for more details on kmrgd*, see fun_eval

% compute the Jacobian matrix using finite differences for the
% nonlinear part of the approximating ordinary equation (which
% is defined by the right-hand side of the delay equation) and
% the exact Jacobian matrix for the linear part (which is defined
% by the differentiation matrix)

global cds;
if isfield(cds,'options') && isfield(cds.options,'Increment') && ~isempty(cds.options.Increment)
    increment = cds.options.Increment;
else
    increment = 1e-5;
end

jac_F_d = NaN(dim_d,length(kmrgd));
%jac_F_r = NaN(dim_r,length(kmrgd));

for ii = 1:length(kmrgd)
    orig_kmrgd_ii = kmrgd(ii);
    kmrgd(ii) = orig_kmrgd_ii-increment;

% values of the derivative of the renewal part of the (integrated) state at all the nodes
%     (M+1 blocks of length dim_r):
%kmrgd_r_der = reshape(transpose(DM(:,2:end)*transpose(reshape(kmrgd,dim_r,M))),(M+1)*dim_r,1); % only renewal equations
%kmrgd_r_der = reshape(transpose(DM(:,2:end)*transpose(reshape(kmrgd(ind_r_all),dim_r,M))),(M+1)*dim_r,1); % coupled differential and renewal equations
% equivalent to kron(DM(:,2:end),eye(dim_r))*kmrgd(ind_r_all)

    % differential part of the right-hand side:
    F_d1 = [
par_beta*kmrgd(dim_d+(M-1)*dim+1)/(1+kmrgd(dim_d+(M-1)*dim+1)^par_n)-par_gamma*kmrgd(1);
    ];
    % renewal part of the right-hand side:
    %F_r1 = [
    %];

    kmrgd(ii) = orig_kmrgd_ii+increment;

% values of the derivative of the renewal part of the (integrated) state at all the nodes
%     (M+1 blocks of length dim_r):
%kmrgd_r_der = reshape(transpose(DM(:,2:end)*transpose(reshape(kmrgd,dim_r,M))),(M+1)*dim_r,1); % only renewal equations
%kmrgd_r_der = reshape(transpose(DM(:,2:end)*transpose(reshape(kmrgd(ind_r_all),dim_r,M))),(M+1)*dim_r,1); % coupled differential and renewal equations
% equivalent to kron(DM(:,2:end),eye(dim_r))*kmrgd(ind_r_all)

    % differential part of the right-hand side:
    F_d2 = [
par_beta*kmrgd(dim_d+(M-1)*dim+1)/(1+kmrgd(dim_d+(M-1)*dim+1)^par_n)-par_gamma*kmrgd(1);
    ];
    % renewal part of the right-hand side:
    %F_r2 = [
    %];

    jac_F_d(:,ii) = F_d2-F_d1;
    %jac_F_r(:,ii) = F_r2-F_r1;

    kmrgd(ii) = orig_kmrgd_ii;
end

jac_F_d = jac_F_d/(2*increment);
%jac_F_r = jac_F_r/(2*increment);

% jac_F_r extended vertically with 0 for the differential part:
% (ind_r_all-dim_d removes the offset in the indices due to the first block in kmrgd)
%jac_F_r0 = repmat(jac_F_r,M,1); % only renewal equations
%jac_F_r0 = zeros(M*dim,length(kmrgd)); jac_F_r0(ind_r_all-dim_d,:) = repmat(jac_F_r,M,1); % coupled differential and renewal equations

% use sparse matrices with kron(...,eye) for performance
Dkron = sparse([],[],[],M*dim,dim_d+M*dim,M*(dim_d+M*dim));
speyedim = speye(dim);
Dkron(:,(dim_d+1):(dim_d+M*dim)) = kron(DM(2:end,2:end), speyedim);
% the first block-column is related only to the differential part
speyedim0 = speyedim(:,ind_d);
Dkron(:,1:dim_d) = kron(DM(2:end,1), speyedim0);

jac = [
jac_F_d;
full(Dkron); % only differential equations
%full(Dkron)-jac_F_r0; % only/also renewal equations
];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_beta,par_gamma,par_n,par_tau)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_beta,par_gamma,par_n,par_tau)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_beta,par_gamma,par_n,par_tau)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_beta,par_gamma,par_n,par_tau)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_beta,par_gamma,par_n,par_tau)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_beta,par_gamma,par_n,par_tau)
% --------------------------------------------------------------------------
function kmrgd = make_kmrgd_from_func(state_funcs, parameters)
% Construct the initial value (vector) from the specified initial function.
%
% state_funcs: cell array of initial functions, in the order of the
%              coordinates; each is a function of one variable in the
%              delay interval
% parameters: vector of the values of the parameters (possibly needed to
%             compute the delays)

M = 5; % degree of collocation polynomials

dim = 1; % total dimension of the delay equation
dim_d = 1; % dimension of the differential part
dim_r = 0; % dimension of the renewal part
ind_d = [1]; % indices of differential equations
%ind_r = []; % indices of renewal equations
% indices of all components of differential part in kmrgd:
ind_d_all = [1,2,3,4,5,6];
% indices of all components of renewal part in kmrgd:
%ind_r_all = [];

dim_kmrgd = dim_d*(M+1)+dim_r*M;
kmrgd = NaN(dim_kmrgd, 1);

par_beta = parameters(1);
par_gamma = parameters(2);
par_n = parameters(3);
par_tau = parameters(4);

delays = -[-par_tau];
if any(delays < 0)
    % the equation is advanced instead of delayed: we need to stop
    kmrgd(:) = NaN;
    return
else
    tau_max = max(delays);
    if tau_max == 0
        % the equation is not delayed: we need to stop
        kmrgd(:) = NaN;
        return
    end
end

% Chebyshev type II nodes (extrema) in [-1, 1]:
unit_nodes = [1;0.809016994374947;0.309016994374947;-0.309016994374947;-0.809016994374947;-1];
nodes = (unit_nodes-1)*tau_max/2; % Chebyshev type II nodes (extrema) in the delay interval [-tau_max, 0]
% weights for barycentric interpolation on Chebyshev type II nodes (extrema) (scaling-independent):
weights = [0.5;-1;1;-1;1;-0.5];
% differentiation matrix for Chebyshev type II nodes (extrema) in [-1, 1]:
unit_DM = [
8.5,-10.4721359549996,2.89442719099992,-1.52786404500042,1.10557280900008,-0.5;
2.6180339887499,-1.17082039324994,-2,0.894427190999916,-0.618033988749895,0.276393202250021;
-0.723606797749979,2,-0.170820393249937,-1.61803398874989,0.894427190999916,-0.381966011250105;
0.381966011250105,-0.894427190999916,1.61803398874989,0.170820393249937,-2,0.723606797749979;
-0.276393202250021,0.618033988749895,-0.894427190999916,2,1.17082039324994,-2.61803398874989;
0.5,-1.10557280900008,1.52786404500042,-2.89442719099992,10.4721359549996,-8.49999999999999;
];
DM = unit_DM*2/tau_max; % differentiation matrix for Chebyshev type II nodes (extrema) in the delay interval [-tau_max, 0]

for i = 1:dim_d
    for j = 1:M+1
        kmrgd(ind_d_all(i+(j-1)*dim_d)) = state_funcs{ind_d(i)}(nodes(j));
    end
end
for i = 1:dim_r
    for j = 1:M
        kmrgd(ind_r_all(i+(j-1)*dim_r)) = state_funcs{ind_r(i)}(nodes(j+1));
    end
    kmrgd(ind_r_all(i:dim_r:end)) = DM(2:end,2:end) \ kmrgd(ind_r_all(i:dim_r:end));
end
% --------------------------------------------------------------------------
function ff = barint(x, w, f, xx)
%BARINT Barycentric interpolation
%  ff = BARINT(x, w, f, xx) computes the values ff of a function at xx
%  using the barycentric interpolation formula with x interpolation
%  nodes, w barycentric weights and f values of the function at x.
%
%  Reference:
%    J.-P. Berrut and L. N. Trefethen,
%    Barycentric Lagrange interpolation,
%    SIAM Review, 46(3):501-517, 2004,
%    DOI: 10.1137/S0036144502417715
if isscalar(xx)
    % most frequent cases
    if xx == x(end)
        ff = f(end);
    elseif xx == x(1)
        ff = f(1);
    else
        temp = w ./ (xx - x);
        numer = dot(temp, f);
        denom = sum(temp);
        ff = numer / denom;
        if isnan(ff) || isinf(ff)
            j = find(xx == x);
            if ~isempty(j)
                ff = f(j);
            end
        end
    end
else
    n = length(x);
    numer = zeros(size(xx));
    denom = zeros(size(xx));
    exact = zeros(size(xx));
    for j = 1:n
        xdiff = xx - x(j);
        temp = w(j) ./ xdiff;
        numer = numer + temp * f(j);
        denom = denom + temp;
        exact(xdiff == 0) = j;
    end
    jj = find(exact);
    ff = numer ./ denom;
    ff(jj) = f(exact(jj));
end
function value = DE_int(f, a, b)
%DE_INT Integral
%  value = DE_INT(f, a, b) computes the integral of f in [a, b]
%  using the Clenshaw-Curtis quadrature formula of degree M = 10
%
%  Reference:
%    L. N. Trefethen,
%    Spectral Methods in MATLAB,
%    SIAM, 2000,
%    DOI: 10.1137/1.9780898719598

% Clenshaw-Curtis quadrature nodes and weights in [-1, 1]:
nodes = [1;0.951056516295154;0.809016994374947;0.587785252292473;0.309016994374947;6.12323399573677e-17;-0.309016994374947;-0.587785252292473;-0.809016994374947;-0.951056516295154;-1];
weights = [0.0101010101010101,0.0945790548837016,0.185635214424248,0.253588333283687,0.299213270424237,0.313766233766234,0.299213270424237,0.253588333283687,0.185635214424248,0.0945790548837016,0.0101010101010101];
% nodes rescaled to [a, b]:
nodes = (b-a)/2*nodes+(a+b)/2;
n = length(nodes);
ff = NaN(n, 1);
for i = 1:length(nodes)
    ff(i) = f(nodes(i));
end
% scaling constant for weights: (b-a)/2
value = (b-a)/2*(weights*ff);

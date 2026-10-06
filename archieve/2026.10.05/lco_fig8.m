%% pio_rle_nidf_gsnm.m
%  Limit Cycle Prediction for Category II PIO Using the Negative Inverse
%  Describing Function of a Rate Limiter: Golden-Section Search for the
%  Common Frequency and Amplitude on the Nyquist Plot
%
%  pio  : Pilot-Induced Oscillation (Category II)
%  rle  : Rate Limiter Element
%  nidf : Negative Inverse Describing Function, -1/N(A,w)
%  gsnm : Golden-Section (Search) Numerical Method
%
%  Cases: X-15 Soft Glide Landing PIO (1959), X-15 Flight 3-65-97 (1967)

clear; clc; close all; clear functions
%% 1. Linear system G(s)
% CTRL+R / CTRL+T to comment and uncomment to switch between two incidents:
% 1)
% -------------------------------------------------------------------------
% THE PARAMETERS OF THE X-15 SOFT GLIDE LANDING PIO INCIDENT in 1959
modelName = 'X-15, The Soft Glide PIO Case, 1959:';
Kp = 1;   % pilot gain (sidestick gradient)
R  = 15;  % rate limit, deg/s
qco     = 1;               % pilot command amplitude
Kp_Yp   = 13.68;           % controller gain (pilot gain inherited)
M_del_e = 0.537;           % input signal entering the rate limiter
omega_n = 2.3;             % natural frequency
zeta_sp = 1.42 / omega_n / 2;   % damping ratio
num = Kp_Yp*M_del_e*[1 0.82];
den = [1 1.42 omega_n^2 0];
Gac = qco*Kp*tf(num,den);  % aircraft transfer function
% -------------------------------------------------------------------------
% 2)
% -------------------------------------------------------------------------
% THE PARAMETERS OF THE X-15 Flight 3-65-97, 1967 FATAL CRASH PIO INCIDENT
% modelName = 'X-15, The Flight\_3-65-97 (The PIO Case caused to the Fatal Crash):';
% Kp = 13.68;   % pilot gain (sidestick gradient)
% R  = 60;      % actuator rate limit, deg/s
% qco = 1;
% num_Gc = 5.21 * conv([1 -57.36], conv([1 4.26], [1 0.55]));
% den_Gc = conv([1 2*0.442*22.85 22.85^2], conv([1 0], [1 1.16]));
% Gc  = tf(num_Gc, den_Gc);
% num_Gac = -10.524 * conv([1 1.562], conv([1 0.038], [1 0]));
% den_Gac = conv([1 2*0.212*0.088 0.088^2], conv([1 3.75], [1 -1.44]));
% Gac = qco*Kp*tf(num_Gac, den_Gac);
% -------------------------------------------------------------------------
% The Linearized Transfer Function
Gs = qco*Kp*Gac;

%% 2. Resolution and the Frequency Range
n = 1000;                  % resolution
w = logspace(-2, 2, n);    % common frequency range, rad/s

%% 3. Rate Limiter Element parameters (OLOP, Category II PIO)
Ai      = qco*Kp;          % input amplitude, deg
w_onset = R / Ai;

%% 4. Describing function N(Ai,w)
N = DF(w / w_onset);       % alpha = w / w_onset = w*Ai/R

%% 5. Calculate -1/N(Ai,w)
minus_inv_N = -1 ./ N;

%% 6. Find Intersections on the Nyquist plane (Polyline Crossing Method, used as starting guess)
% G(jw) and -1/N(Ai,w) are treated as polylines in the complex plane.
% Every segment pair that crosses gives one starting frequency on G.
G_curve = squeeze(freqresp(Gs, w)).';     % G(jw) as a row vector
w_cross = polyCross(G_curve, minus_inv_N, w);

%% 7. Golden-Section Search for the Common Frequency
tol = 1e-8;                % search tolerance
r   = (sqrt(5) - 1) / 2;   % golden ratio factor (about 0.618)

w_cross = sort(w_cross);   % starting guesses, low to high
nX      = length(w_cross);

w_star_all = zeros(1, nX);
A_star_all = zeros(1, nX);
G_star_all = zeros(1, nX);
dist_all   = zeros(1, nX);
w_N_Ai_all = zeros(1, nX);

for i = 1:nX
    % Bracket: each side ends halfway to the neighbouring crossing
    if i == 1
        a = 0.7 * w_cross(1);
    else
        a = (w_cross(i-1) + w_cross(i)) / 2;
    end
    if i == nX
        b = 1.4 * w_cross(nX);
    else
        b = (w_cross(i) + w_cross(i+1)) / 2;
    end

    % Outer search over w; inner search (bestAlpha) over alpha = w*A/R
    while (b - a) > tol
        w1 = b - r*(b - a);
        w2 = a + r*(b - a);
        d1 = curveDistance(Gs, w1, tol);
        d2 = curveDistance(Gs, w2, tol);
        if d1 < d2
            b = w2;
        else
            a = w1;
        end
    end

    w_star     = (a + b) / 2;
    G_star     = squeeze(freqresp(Gs, w_star));
    alpha_star = bestAlpha(G_star, tol);

    w_star_all(i) = w_star;
    A_star_all(i) = alpha_star * R / w_star;    % amplitude giving the common frequency
    G_star_all(i) = G_star;                     % intersection point in the complex plane
    dist_all(i)   = curveDistance(Gs, w_star, tol);
    w_N_Ai_all(i) = alpha_star * R / Ai;        % where the fixed-Ai curve hits this point
end

%% 8. Data tip values and Command Window output
% The same segment data is used for every data tip, so the printed values
% are exactly the values shown in the data tips on the figures.
seg = cell(1, nX);
for i = 1:nX
    seg{i} = matchedSegments(Gs, w_star_all(i), A_star_all(i), R, 0.02);
end

disp(modelName)
fprintf('\n--- Limit Cycle Intersections (Golden-Section Search) ---\n');
if nX == 0
    fprintf('No intersections found.\n');
end
for i = 1:nX
    s  = seg{i};
    k  = s.i;
    nm = sprintf('-1/N(A*=%.4f)', s.A);
    fprintf('\nIntersection %d   (A* = %.4f deg)\n', i, s.A);
    fprintf('  %-20s %14s %16s %20s\n', 'Response', 'Real', 'Imaginary', 'Frequency (rad/s)');
    fprintf('  %-20s %14.4f %16.4f %20.4f\n', 'G(jw)', real(s.G(k)), imag(s.G(k)), s.w(k));
    fprintf('  %-20s %14.4f %16.4f %20.4f\n', nm, real(s.N(k)), imag(s.N(k)), s.w(k));
    fprintf('  Distance at solution : %.2e\n', dist_all(i));
    fprintf('  With fixed Ai = %.4f, -1/N(Ai,w) reaches this point at %.4f rad/s\n', Ai, w_N_Ai_all(i));
end
fprintf('\n');

%% 9. Screen layout: two figures side by side
scr = get(groot, 'ScreenSize');      % [1 1 width height] of the main screen
W   = scr(3);
yb  = 40;                            % space left for the taskbar
H   = scr(4) - yb;

pos1 = [1,          yb, 0.45*W, H];  % left  : Intersections
pos2 = [0.45*W + 1, yb, 0.55*W, H];  % right : Zoom In / Zoom Out

% Nyquist options: positive frequencies only (no mirrored branch)
opts = nyquistoptions;
opts.ShowFullContour = 'off';

if nX > 0
    %% 10. Figure 1: Intersection graphs, one row per intersection
    figure(Name='Intersections',NumberTitle='off',OuterPosition=pos1);
    t1 = tiledlayout(nX, 1, 'TileSpacing', 'tight', 'Padding', 'compact');
    title(t1, {'The Nyquist Plot of G(j\omega) and -1/N(A_i,\omega)', '(The Negative Inverse Describing Function Technique)'})
    subtitle(t1, modelName, 'Interpreter', 'latex')

    for i = 1:nX
        w_star = w_star_all(i);
        A_star = A_star_all(i);
        Gp     = G_star_all(i);

        w_plot     = unique([w, w_star]);
        sys_G      = frd(squeeze(freqresp(Gs, w_plot)), w_plot);
        sys_N_star = frd(reshape(-1 ./ DF(w_plot * A_star / R), 1, 1, []), w_plot);

        nexttile(t1)
        nyquistplot(sys_G, w_plot, opts);
        hold on
        nyquistplot(sys_N_star, w_plot, opts);
        axI = gca;
        grid on
        xlim(real(Gp) + [-3 3])
        ylim(imag(Gp) + [-3 3])
        addTips(axI, seg{i});
        plot(axI, real(Gp), imag(Gp), 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Intersection')
        title(sprintf('Intersection %d:  \\omega^* = %.4f rad/s,  A^* = %.4f', i, w_star, A_star))
        legend('G(j\omega)','-1/N(A^*,\omega)','Intersection','Location','best','FontSize',7)
    end

    %% 11. Figure 2: Zoom In (left) and Zoom Out (right), one row per intersection
    figure(Name='Zoom In / Zoom Out',NumberTitle='off',OuterPosition=pos2);
    t2 = tiledlayout(nX, 2, 'TileSpacing', 'tight', 'Padding', 'compact');
    title(t2, {'The Nyquist Plot of G(j\omega) and -1/N(A_i,\omega)', '(The Negative Inverse Describing Function Technique)'})
    subtitle(t2, modelName, 'Interpreter', 'latex')

    for i = 1:nX
        w_star = w_star_all(i);
        A_star = A_star_all(i);
        Gp     = G_star_all(i);

        % Zoom In i (left tile)
        w_zoom     = unique([w_star*(1 + (-0.01:0.0005:0.01)), w_star]);
        sys_G_zoom = frd(squeeze(freqresp(Gs, w_zoom)), w_zoom);
        sys_N_zoom = frd(reshape(-1 ./ DF(w_zoom * A_star / R), 1, 1, []), w_zoom);

        nexttile(t2)
        nyquistplot(sys_G_zoom, w_zoom, opts);
        hold on
        nyquistplot(sys_N_zoom, w_zoom, opts);
        axZ = gca;
        addTips(axZ, seg{i});
        plot(axZ, real(Gp), imag(Gp), 'ko', 'MarkerSize', 12, 'LineWidth', 1.5, 'DisplayName', 'Intersection')
        grid on
        xlim(real(Gp) + [-0.01 0.01])
        ylim(imag(Gp) + [-0.01 0.01])
        title(sprintf('Zoom In %d', i))
        legend('G(j\omega)','-1/N(A^*,\omega)','Intersection','Location','best','FontSize',7)

        % Zoom Out i (right tile)
        w_plot     = unique([w, w_star]);
        sys_G      = frd(squeeze(freqresp(Gs, w_plot)), w_plot);
        sys_N_star = frd(reshape(-1 ./ DF(w_plot * A_star / R), 1, 1, []), w_plot);

        nexttile(t2)
        nyquistplot(sys_G, w_plot, opts);
        hold on
        nyquistplot(sys_N_star, w_plot, opts);
        axO = gca;
        grid on
        xlim([-10 10])
        ylim([-10 10])
        addTips(axO, seg{i});
        plot(axO, real(Gp), imag(Gp), 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Intersection')
        title(sprintf('Zoom Out %d', i))
        legend('G(j\omega)','-1/N(A^*,\omega)','Intersection','Location','best','FontSize',7)
    end
end

%% Functions

% Rate limiter describing function N(alpha), alpha = w*A/R
function N = DF(alpha)
    N = zeros(size(alpha));
    for k = 1:length(alpha)
        al = alpha(k);
        if al <= 1                      % Region I
            M = 1;
            phi = 0;
        elseif al < 1.862               % Region II
            M   = 0.2908*al^3 - 1.4396*al^2 + 1.9232*al + 0.223;
            phi = 0.528*al^3 - 2.6213*al^2 + 3.5056*al - 1.4171;
        else                            % Region III
            varpi = 1 / al;
            M   = (4/pi)*varpi;
            phi = -acos((pi/2)*varpi);
        end
        N(k) = M*exp(1j*phi);
    end
end

% Crossings of two polylines P and Q in the complex plane.
% Returns the frequency on P (interpolated along the segment) at each crossing.
function w_cross = polyCross(P, Q, wP)
    % Segment end points of P (columns) and Q (rows)
    x1 = real(P(1:end-1)).';  y1 = imag(P(1:end-1)).';
    x2 = real(P(2:end)).';    y2 = imag(P(2:end)).';
    x3 = real(Q(1:end-1));    y3 = imag(Q(1:end-1));
    x4 = real(Q(2:end));      y4 = imag(Q(2:end));

    % Line-segment intersection parameters t (along P) and u (along Q)
    den = (x1 - x2).*(y3 - y4) - (y1 - y2).*(x3 - x4);
    t   = ((x1 - x3).*(y3 - y4) - (y1 - y3).*(x3 - x4)) ./ den;
    u   = -((x1 - x2).*(y1 - y3) - (y1 - y2).*(x1 - x3)) ./ den;

    hit = (den ~= 0) & (t >= 0) & (t < 1) & (u >= 0) & (u < 1);
    [iP, ~] = find(hit);
    tt = t(hit);

    w_cross = wP(iP) + tt.' .* (wP(iP + 1) - wP(iP));
    w_cross = unique(w_cross);
end

% Golden-section search over alpha: point on -1/N closest to G(jw)
function alpha = bestAlpha(Gpt, tol)
    r = (sqrt(5) - 1) / 2;
    a = 1;                              % start of rate limiting
    b = 20;                             % upper bound on alpha
    while (b - a) > tol
        al1 = b - r*(b - a);
        al2 = a + r*(b - a);
        d1 = abs(Gpt - (-1 / DF(al1)));
        d2 = abs(Gpt - (-1 / DF(al2)));
        if d1 < d2
            b = al2;
        else
            a = al1;
        end
    end
    alpha = (a + b) / 2;
end

% Smallest distance between G(jw) and the -1/N curve at frequency w
function d = curveDistance(Gs, w, tol)
    Gpt   = squeeze(freqresp(Gs, w));
    alpha = bestAlpha(Gpt, tol);
    d     = abs(Gpt - (-1 / DF(alpha)));
end

% Short segments of G(jw) and -1/N(A*,w) around w* in the complex plane
function s = matchedSegments(Gs, w_star, A_star, R, frac)
    s.w = unique([w_star*(1 + linspace(-frac, frac, 81)), w_star]);
    s.i = find(s.w == w_star);
    s.A = A_star;
    s.G = squeeze(freqresp(Gs, s.w)).';
    s.N = -1 ./ DF(s.w * A_star / R);
end

% Draw the segments and put data tips at w* on both curves
function addTips(ax, s)
    blue   = [0 0.447 0.741];
    orange = [0.850 0.325 0.098];

    hG = plot(ax, real(s.G), imag(s.G), '-', 'Color', blue, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    hN = plot(ax, real(s.N), imag(s.N), '-', 'Color', orange, 'LineWidth', 2.5, 'HandleVisibility', 'off');

    setTemplate(hG, 'G(jw)', s.w);
    setTemplate(hN, sprintf('-1/N(A*=%.4f)', s.A), s.w);

    dtG = datatip(hG, 'DataIndex', s.i);
    dtN = datatip(hN, 'DataIndex', s.i);
    dtG.Location = 'northwest';
    dtN.Location = 'southeast';
    dtG.FontSize = 7;                   % smaller tips so they fit inside the tiles
    dtN.FontSize = 7;
end

% Data tip rows: Response, Real, Imaginary, Frequency (4 decimals)
function setTemplate(h, name, wv)
    rows = h.DataTipTemplate.DataTipRows;
    rows(1).Label  = 'Real';
    rows(1).Format = '%.4f';
    rows(2).Label  = 'Imaginary';
    rows(2).Format = '%.4f';
    respRow = dataTipTextRow('Response', repmat({name}, 1, numel(wv)));
    freqRow = dataTipTextRow('Frequency (rad/s)', wv, '%.4f');
    h.DataTipTemplate.DataTipRows = [respRow; rows(1); rows(2); freqRow];
end
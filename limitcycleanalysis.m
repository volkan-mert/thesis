%% limitcycleanalysis.m
%  Limit Cycle Prediction for Category II PIO Using the Negative Inverse
%  Describing Function of a Rate Limiter: Golden-Section Search for the
%  Common Frequency and Amplitude on the Nyquist Plot
%
%  pio  : Pilot-Induced Oscillation (Category II)
%  rle  : Rate Limiter Element
%  nidf : Negative Inverse Describing Function, -1/N(A,w)
%  gsnm : Golden-Section (Search) Numerical Method
%
%  Cases: 
%   1.  X-15 Soft Glide Landing PIO (1959) 
%   2.  X-15 Flight 3-65-97 (1967)

clear; clc; close all; clear functions
tic                                  % total run time
%% 1. Linear system G(s)
% CTRL+R / CTRL+T to comment and uncomment to switch between two incidents:
% 1)
% -------------------------------------------------------------------------
% THE PARAMETERS OF THE X-15 SOFT GLIDE LANDING PIO INCIDENT in 1959
modelName = 'X-15, The Soft Glide PIO Case, 1959';
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
% modelName = 'X-15, The Flight 3-65-97 (The PIO Case caused to the Fatal Crash)';
% Kp = 13.68;   % pilot gain (sidestick gradient)
% R  = 60;      % actuator rate limit, deg/s
% qco = 1;
% num_Gc = 5.21 * conv([1 -57.36], conv([1 4.26], [1 0.55]));
% den_Gc = conv([1 2*0.442*22.85 22.85^2], conv([1 0], [1 1.16]));
% Gc  = tf(num_Gc, den_Gc);
% num_Gac = -10.524 * conv([1 1.562], conv([1 0.038], [1 0]));
% den_Gac = conv([1 2*0.212*0.088 0.088^2], conv([1 3.75], [1 -1.44]));
% Gac = qco*Kp*tf(num_Gac, den_Gac);
% % -------------------------------------------------------------------------
% The Linearized Transfer Function
Gs = qco*Kp*Gac;

% Fast evaluation of G(jw) from its polynomials (much cheaper than freqresp per call)
[numG, denG] = tfdata(Gs, 'v');
evalG = @(wv) polyval(numG, 1j*wv) ./ polyval(denG, 1j*wv);

%% 2. Resolution and the Frequency Range (for the plots)
n = 1000;                  % resolution
w = logspace(-2, 2, n);    % plotting frequency range, rad/s

%% 3. Rate Limiter Element parameters (OLOP, Category II PIO)
Ai = qco*Kp;               % input amplitude, deg (used for comparison only)

%% 4. Scan of Frequencies and Amplitudes (vectorized over the whole (w, A) grid)
w_scan = 0.01:0.1:100;     % scanned frequencies, rad/s
A_scan = 0.01:0.1:100;     % scanned amplitudes, deg
nW = numel(w_scan);
nA = numel(A_scan);

G_scan = evalG(w_scan);                      % G(jw) at the scanned frequencies (1 x nW)

alpha_grid = (w_scan.' * A_scan) / R;        % alpha = w*A/R  (rows: w, columns: A)
N_scan     = -1 ./ DF(alpha_grid);           % -1/N(A,w) for every (w, A)
D          = abs(G_scan.' - N_scan);         % distance |G(jw) - (-1/N(A,w))|

[d_min, kmin] = min(D, [], 2);               % smallest distance and its A index at each w
d_min  = d_min.';
A_best = A_scan(kmin);

side   = NaN(1, nW);       % side of the -1/N curve that G lies on (+1 or -1)
seglen = zeros(1, nW);     % local length of the -1/N curve near the best A
for iw = 1:nW
    k = kmin(iw);
    if k > 1 && k < nA
        tang = N_scan(iw,k+1) - N_scan(iw,k-1);
        if abs(tang) > 0
            side(iw)   = 1 - 2*(imag(conj(tang) * (G_scan(iw) - N_scan(iw,k))) < 0);
            seglen(iw) = abs(tang);
        end
    end
end

% A side change between two neighbouring frequencies brackets one intersection
brackets = zeros(0, 2);
A_scanX  = zeros(1, 0);
for iw = 1:nW-1
    if ~isnan(side(iw)) && ~isnan(side(iw+1)) && side(iw) ~= side(iw+1)
        lim = abs(G_scan(iw+1) - G_scan(iw)) + max(seglen(iw), seglen(iw+1));
        if d_min(iw) <= lim && d_min(iw+1) <= lim
            brackets(end+1, :) = [w_scan(iw), w_scan(iw+1)];
            if d_min(iw) <= d_min(iw+1)
                A_scanX(end+1) = A_best(iw);
            else
                A_scanX(end+1) = A_best(iw+1);
            end
        end
    end
end
nX = size(brackets, 1);

%% 5. Golden-Section Search for the Common Frequency
tol = 1e-8;                % search tolerance
r   = (sqrt(5) - 1) / 2;   % golden ratio factor (about 0.618)

w_star_all = zeros(1, nX);
A_star_all = zeros(1, nX);
G_star_all = zeros(1, nX);
dist_all   = zeros(1, nX);
w_N_Ai_all = zeros(1, nX);

for i = 1:nX
    % Bracket from the scan
    a = brackets(i, 1);
    b = brackets(i, 2);

    % Outer search over w (one new evaluation per iteration);
    % inner search (bestAlpha) over alpha = w*A/R
    w1 = b - r*(b - a);   d1 = curveDistance(numG, denG, w1, tol);
    w2 = a + r*(b - a);   d2 = curveDistance(numG, denG, w2, tol);
    while (b - a) > tol
        if d1 < d2                      % minimum is in [a, w2]
            b  = w2;
            w2 = w1;  d2 = d1;
            w1 = b - r*(b - a);  d1 = curveDistance(numG, denG, w1, tol);
        else                            % minimum is in [w1, b]
            a  = w1;
            w1 = w2;  d1 = d2;
            w2 = a + r*(b - a);  d2 = curveDistance(numG, denG, w2, tol);
        end
    end

    w_star     = (a + b) / 2;
    G_star     = evalG(w_star);
    alpha_star = bestAlpha(G_star, tol);

    w_star_all(i) = w_star;
    A_star_all(i) = alpha_star * R / w_star;    % amplitude giving the common frequency
    G_star_all(i) = G_star;                     % intersection point in the complex plane
    dist_all(i)   = abs(G_star - (-1 / DF(alpha_star)));
    w_N_Ai_all(i) = alpha_star * R / Ai;        % where the fixed-Ai curve hits this point
end

%% 6. Data tip values and Command Window output
% The same segment data is used for every data tip, so the printed values
% are exactly the values shown in the data tips on the figure.
seg = cell(1, nX);
for i = 1:nX
    seg{i} = matchedSegments(evalG, w_star_all(i), A_star_all(i), R, 0.02);
end

disp(modelName)
fprintf('\n----------------------------------------------------------------------------------------------------------------\n');
fprintf('\n    *************************************************   \n');
fprintf('\n    * Limit Cycle Occurences: (Intersection Points) *   \n');
fprintf('\n    *************************************************   \n');
fprintf('\n----------------------------------------------------------------------------------------------------------------\n');
if nX == 0
    fprintf('No intersections found.\n');
end
for i = 1:nX
    s  = seg{i};
    k  = s.i;
    nm = sprintf('-1/N(A*=%.4f)', s.A);
    fprintf('\nLimit Cycle Occurence  %d \n', i);
    fprintf('\nIntersection Point %d   (A* = %.4f deg)\n', i, s.A);
    fprintf('  Scan bracket         : [%.2f, %.2f] rad/s,  A from scan = %.2f deg\n', brackets(i,1), brackets(i,2), A_scanX(i));
    fprintf('  %-20s %14s %16s %20s\n', 'Response', 'Real', 'Imaginary', 'Frequency (rad/s)');
    fprintf('  %-20s %14.4f %16.4f %20.4f\n', 'G(jw)', real(s.G(k)), imag(s.G(k)), s.w(k));
    fprintf('  %-20s %14.4f %16.4f %20.4f\n', nm, real(s.N(k)), imag(s.N(k)), s.w(k));
    fprintf('  Distance at solution : %.2e\n', dist_all(i));
    fprintf('  With fixed Ai = %.4f, -1/N(Ai,w) reaches this point at %.4f rad/s\n', Ai, w_N_Ai_all(i));
    fprintf('\n----------------------------------------------------------------------------------------------------------------\n');
end
fprintf('\n');

%% 7. Screen layout, Nyquist options, scan colors and title style
scr = get(groot, 'ScreenSize');      % [1 1 width height] of the main screen
W   = scr(3);
yb  = 40;                            % space left for the taskbar
H   = scr(4) - yb;

posZ = [1, yb, W, H];                % whole screen : Zoom In / Zoom Out

% Nyquist options: positive frequencies only (no mirrored branch)
opts = nyquistoptions;
opts.ShowFullContour = 'off';

nBins  = 64;                         % frequency color bands (one line object per band)
cmapB  = parula(nBins);
famLbl = '-1/N(A,\omega),  A = 0.01:0.1:100';

% Title style (LaTeX interpreter)
fsMain = 13;                         % main figure title
fsSub  = 11;                         % model name subtitle
fsTile = 11;                         % intersection and zoom titles
mainTitle = {'The Nyquist Plot of $G(j\omega)$ and $-1/N(A_i,\omega)$', '(The Negative Inverse Describing Function Technique)'};

if nX > 0
    %% 8. Figure: Zoom In (left) and Zoom Out (right), one titled row per intersection
    % Scanned curves for the Zoom Out window are the same for every row: build them once
    S_out = buildScan(N_scan, w_scan, A_scan, 0, 10, nBins);

    fig2 = figure(Name='Zoom In / Zoom Out',NumberTitle='off',OuterPosition=posZ);
    t2 = tiledlayout(fig2, nX, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    t2.OuterPosition = [0 0 0.92 1];
    title(t2, mainTitle, 'Interpreter', 'latex', 'FontSize', fsMain)
    subtitle(t2, modelName, 'Interpreter', 'latex', 'FontSize', fsSub)

    for i = 1:nX
        w_star = w_star_all(i);
        A_star = A_star_all(i);
        Gp     = G_star_all(i);

        % Row i: inner 1x2 layout carrying the intersection title
        tRow = tiledlayout(t2, 1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
        tRow.Layout.Tile = i;
        title(tRow, sprintf('Intersection %d:\\quad $\\omega^* = %.4f$ rad/s,\\quad $A^* = %.4f$', i, w_star, A_star), 'Interpreter', 'latex', 'FontSize', fsTile)

        % Zoom In i (left tile)
        w_zoom     = unique([w_star*(1 + (-0.01:0.0005:0.01)), w_star]);
        sys_G_zoom = frd(evalG(w_zoom), w_zoom);
        S_in       = buildScan(N_scan, w_scan, A_scan, Gp, 0.01, nBins);

        titledCell(tRow, 1, sprintf('Zoom In %d', i), fsTile);
        nyquistplot(sys_G_zoom, w_zoom, opts);
        hold on
        axZ = gca;
        drawBins(axZ, S_in, cmapB, R);
        plot(axZ, NaN, NaN, '-', 'Color', cmapB(round(nBins/2),:), 'LineWidth', 2)   % legend sample
        addTips(axZ, seg{i});
        plot(axZ, real(Gp), imag(Gp), 'ko', 'MarkerSize', 12, 'LineWidth', 1.5)
        grid on
        xlim(real(Gp) + [-0.01 0.01])
        ylim(imag(Gp) + [-0.01 0.01])
        title('')
        legend(axZ, 'G(j\omega)', famLbl, 'Intersection', 'Location', 'best', 'FontSize', 7)

        % Zoom Out i (right tile)
        w_plot = unique([w, w_star]);
        sys_G  = frd(evalG(w_plot), w_plot);

        titledCell(tRow, 2, sprintf('Zoom Out %d', i), fsTile);
        nyquistplot(sys_G, w_plot, opts);
        hold on
        axO = gca;
        drawBins(axO, S_out, cmapB, R);
        plot(axO, NaN, NaN, '-', 'Color', cmapB(round(nBins/2),:), 'LineWidth', 2)   % legend sample
        addTips(axO, seg{i});
        plot(axO, real(Gp), imag(Gp), 'ko', 'MarkerSize', 10, 'LineWidth', 1.5)
        grid on
        xlim([-10 10])
        ylim([-10 10])
        title('')
        legend(axO, 'G(j\omega)', famLbl, 'Intersection', 'Location', 'best', 'FontSize', 7)
    end
    addOmegaColorbar(fig2, [w_scan(1) w_scan(end)], cmapB);
    drawnow
end

fprintf('Total run time: %.2f s\n\n', toc);

%% Functions

% Rate limiter describing function N(alpha), alpha = w*A/R (vectorized)
function N = DF(alpha)
    M   = ones(size(alpha));            % Region I: alpha <= 1
    phi = zeros(size(alpha));

    r2 = alpha > 1 & alpha < 1.862;     % Region II
    al = alpha(r2);
    M(r2)   = 0.2908*al.^3 - 1.4396*al.^2 + 1.9232*al + 0.223;
    phi(r2) = 0.528*al.^3 - 2.6213*al.^2 + 3.5056*al - 1.4171;

    r3 = alpha >= 1.862;                % Region III
    varpi   = 1 ./ alpha(r3);
    M(r3)   = (4/pi)*varpi;
    phi(r3) = -acos((pi/2)*varpi);

    N = M .* exp(1j*phi);
end

% Golden-section search over alpha: point on -1/N closest to G(jw)
% (one new evaluation per iteration)
function alpha = bestAlpha(Gpt, tol)
    r = (sqrt(5) - 1) / 2;
    a = 1;                              % start of rate limiting
    b = 20;                             % upper bound on alpha
    al1 = b - r*(b - a);   d1 = abs(Gpt + 1 / DF(al1));
    al2 = a + r*(b - a);   d2 = abs(Gpt + 1 / DF(al2));
    while (b - a) > tol
        if d1 < d2
            b   = al2;
            al2 = al1;  d2 = d1;
            al1 = b - r*(b - a);  d1 = abs(Gpt + 1 / DF(al1));
        else
            a   = al1;
            al1 = al2;  d1 = d2;
            al2 = a + r*(b - a);  d2 = abs(Gpt + 1 / DF(al2));
        end
    end
    alpha = (a + b) / 2;
end

% Smallest distance between G(jw) and the -1/N curve at frequency w
function d = curveDistance(numG, denG, w, tol)
    Gpt   = polyval(numG, 1j*w) / polyval(denG, 1j*w);
    alpha = bestAlpha(Gpt, tol);
    d     = abs(Gpt + 1 / DF(alpha));
end

% Short segments of G(jw) and -1/N(A*,w) around w* in the complex plane
function s = matchedSegments(evalG, w_star, A_star, R, frac)
    s.w = unique([w_star*(1 + linspace(-frac, frac, 81)), w_star]);
    s.i = find(s.w == w_star);
    s.A = A_star;
    s.G = evalG(s.w);
    s.N = -1 ./ DF(s.w * A_star / R);
end

% One tile holding a 1x1 layout with a LaTeX title; the next chart goes inside it
function titledCell(parent, tileIdx, txt, fs)
    tc = tiledlayout(parent, 1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
    tc.Layout.Tile = tileIdx;
    title(tc, txt, 'Interpreter', 'latex', 'FontSize', fs)
    nexttile(tc)
end

% Collect the scanned -1/N(A,w) curves near a view window (center ctr, half-width hw)
% into nBins frequency bands. Each band is one line with NaN gaps between curves.
function S = buildScan(N_scan, w_scan, A_scan, ctr, hw, nBins)
    [nW, nA] = size(N_scan);
    sz   = [nW, nA];
    rows = (1:nW).';

    % Closest scanned point to the window center and the local segment length (all curves at once)
    [dk, k] = min(abs(N_scan - ctr), [], 2);
    kL = max(k - 1, 1);
    kR = min(k + 1, nA);
    Nk = N_scan(sub2ind(sz, rows, k));
    sg = max(abs(N_scan(sub2ind(sz, rows, kR)) - Nk), abs(Nk - N_scan(sub2ind(sz, rows, kL))));
    keep = find(dk <= sqrt(2)*hw + sg);         % curves that come near the window

    inBox = abs(real(N_scan) - real(ctr)) <= 2*hw & abs(imag(N_scan) - imag(ctr)) <= 2*hw;
    band  = min(floor((w_scan - w_scan(1)) / (w_scan(end) - w_scan(1)) * nBins) + 1, nBins);

    parts = cell(1, nBins);
    for iw = keep.'
        inb = find(inBox(iw, :));
        i1  = max(min([inb, kL(iw)]) - 1, 1);
        i2  = min(max([inb, kR(iw)]) + 1, nA);
        idx = i1:i2;
        nP  = numel(idx);
        Nc  = N_scan(iw, idx);
        parts{band(iw)}{end+1} = [real(Nc), NaN; imag(Nc), NaN; repmat(w_scan(iw), 1, nP), NaN; A_scan(idx), NaN];
    end

    S = struct('X', cell(1, nBins), 'Y', cell(1, nBins), 'W', cell(1, nBins), 'A', cell(1, nBins));
    for b = 1:nBins
        if ~isempty(parts{b})
            M = [parts{b}{:}];
            S(b).X = M(1, :);
            S(b).Y = M(2, :);
            S(b).W = M(3, :);
            S(b).A = M(4, :);
        end
    end
end

% Draw the frequency bands (high w first, so low-w curves stay visible on top)
function drawBins(ax, S, cmapB, R)
    for b = numel(S):-1:1
        if isempty(S(b).X)
            continue
        end
        h = plot(ax, S(b).X, S(b).Y, '-', 'Color', cmapB(b,:), 'LineWidth', 1, 'HandleVisibility', 'off');
        setScanTemplate(h, S(b).W, S(b).A, R);
    end
end

% Colorbar for the scanned frequency, on the right side of a figure
function addOmegaColorbar(fig, wlim, cmapB)
    axCB = axes(fig, 'Position', [0.93 0.12 0.001 0.76], 'Visible', 'off', 'HitTest', 'off');
    colormap(axCB, cmapB)
    clim(axCB, wlim)
    cb = colorbar(axCB, 'Position', [0.945 0.12 0.015 0.76]);
    cb.Label.String = '\omega (rad/s)';
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

% Data tip rows for scanned -1/N(A,w) points (per-point w and A, NaN at curve gaps)
function setScanTemplate(h, wv, Av, R)
    nP   = numel(wv);
    rows = h.DataTipTemplate.DataTipRows;
    rows(1).Label  = 'Real';
    rows(1).Format = '%.4f';
    rows(2).Label  = 'Imaginary';
    rows(2).Format = '%.4f';
    respRow  = dataTipTextRow('Response', repmat({'-1/N(A,w) scan'}, 1, nP));
    freqRow  = dataTipTextRow('Frequency (rad/s)', wv, '%.4f');
    ampRow   = dataTipTextRow('Amplitude A (deg)', Av, '%.4f');
    alphaRow = dataTipTextRow('alpha = wA/R', wv .* Av / R, '%.4f');
    h.DataTipTemplate.DataTipRows = [respRow; rows(1); rows(2); freqRow; ampRow; alphaRow];
end
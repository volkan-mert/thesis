%  Rate-limiting element on the NYQUIST plane:  G(jw)  vs.  -1/N
%
%  Key property:  -1/N = (pi/(4*varpi)) * exp(1i*(pi + psi)),  psi = acos((pi/2)*varpi)
%                 Re{-1/N} = -(pi/(4*varpi)) * cos(psi) = -(pi/(4*varpi))*(pi/2)*varpi
%                          = -pi^2/8   ==  CONSTANT
%  => the -1/N locus is a VERTICAL LINE at Re = -pi^2/8 = -1.2337, going down
%     from the real axis to -j*inf.  It does not depend on A, on w, or on R.
%
%  Limit cycle  <=>  G(jw) = -1/N  <=>  Re{G(jw)} = -pi^2/8  AND  Im{G(jw)} < 0
clear; clc; close all

%% ---------------------------------------------------------------- 1. Plant
num_claw = 5.21 * [1, -52.55, -273.6, -134.4];
den_claw = [1, 21.36, 545.6, 605.7, 0];
G_claw   = tf(num_claw, den_claw);

num_ac = -10.524 * [1, 1.6, 0.059, 0];
den_ac = [1, 2.35, -5.31, 0.184, -0.041];
G_ac   = tf(num_ac, den_ac);

Kp = 1;                                   % pilot gain (pilot model inside the loop)
G  = Kp * G_claw * G_ac;

R     = 15;                               % rate (slew) limit [deg/s]
% u_rle = [0.2, 0.3, 0.31, 0.5, 1, 5];      % R.L.E. input amplitudes
u_rle = 0.3;      % R.L.E. input amplitudes
w     = logspace(-1, 2, 2000);

CRIT = -pi^2/8;                           % = -1.2337, the -1/N vertical line

%% ------------------------------------------------- 2. -1/N as FRD objects
% NOTE the min(1,...) clamp.  Without it, acos() of an argument > 1 returns a
% COMPLEX number for varpi in (2/pi, 1]  (i.e. A*w between R and pi*R/2), and
% the "N(A*w <= R) = 1" line does NOT overwrite that band.  That was the bug.
sys_inv_list = cell(1, numel(u_rle));
for k = 1:numel(u_rle)
    A_i   = u_rle(k);
    varpi = R ./ (A_i * w);
    N     = (4/pi) * varpi .* exp(-1i * acos(min(1, (pi/2)*varpi)));   % clamped
    N((A_i * w) <= R) = 1;                                            % linear region
    sys_inv_list{k} = frd(reshape(-1./N, 1, 1, []), w);
end

%% ------------------------------------------------------- 3. Nyquist figure
figure('Color','w');

opt = nyquistoptions;
opt.Grid            = 'on';
opt.ShowFullContour = 'off';              % suppress the negative-frequency mirror

h = nyquistplot(G, sys_inv_list{:}, w);
setoptions(h, opt);

ax = gca; hold(ax,'on');

% analytic -1/N locus drawn on top (all FRD curves collapse onto this line)
plot(ax, [CRIT CRIT], [-6 0], 'r--', 'LineWidth', 2);          % saturated branch
plot(ax, [-pi/4 CRIT], [0 0],  'r--', 'LineWidth', 2);         % onset segment
plot(ax, -1, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 1.5);     % critical point
text(ax, CRIT, -5.6, sprintf('  Re = -\\pi^2/8 = %.4f', CRIT), ...
     'Color','r','FontWeight','bold','FontSize',9);

xlim(ax, [-6 3]); ylim(ax, [-5 3]);
xlabel(ax,'Real Axis'); ylabel(ax,'Imaginary Axis');
title(ax, {'Nyquist Plot of a Rate Limiting Element', ...
           'Open-Loop Plant G(j\omega) vs. the -1/N Locus', ...
           ['\color{red}(Slew Rate, R = ' num2str(R) ' deg/s)']});

%% ------------------------------- 4. Exact intersection (harmonic balance)
% On the Nyquist plane this is a SINGLE scalar equation:  Re{G(jw)} + pi^2/8 = 0
Gv  = squeeze(freqresp(G, w)).';
res = real(Gv) - CRIT;
idx = find(sign(res(1:end-1)).*sign(res(2:end)) < 0);

fprintf('\n=============== LIMIT CYCLE INTERSECTION  (Re{G} = -pi^2/8) ===============\n');
fprintf('%-14s %-12s %-10s %-22s %-10s\n', ...
        'w_lc (rad/s)','f_lc (Hz)','varpi','G(jw_lc)','A_lc');
fprintf('---------------------------------------------------------------------------\n');

n_lc = 0;  w_lc = NaN;
for k = 1:numel(idx)
    fun = @(x) real(squeeze(freqresp(G, x))) - CRIT;
    try
        wk = fzero(fun, [w(idx(k)) w(idx(k)+1)]);
    catch
        continue
    end

    Gk = squeeze(freqresp(G, wk));

    % the locus is only the LOWER half of the line -> reject Im >= 0
    if imag(Gk) >= 0
        fprintf('   (w = %.5f rejected: Im{G} = %+.4f > 0, upper half-plane)\n', wk, imag(Gk));
        continue
    end

    vk = pi/(4*abs(Gk));                  % from |G| = pi/(4*varpi)
    Ak = R/(vk*wk);                       % amplitude that closes the loop

    Nk = (4/pi)*vk*exp(-1i*acos(min(1,(pi/2)*vk)));
    if abs(Gk*Nk + 1) > 1e-6 || ~isfinite(Ak) || Ak <= 0
        continue
    end

    n_lc = n_lc + 1;  w_lc = wk;
    fprintf('%-14.5f %-12.5f %-10.5f %+8.5f %+8.5fj      %-10.4f\n', ...
            wk, wk/(2*pi), vk, real(Gk), imag(Gk), Ak);

    plot(ax, real(Gk), imag(Gk), 'rp', 'MarkerSize', 16, 'MarkerFaceColor', 'r');
    text(ax, real(Gk), imag(Gk), sprintf( ...
        ['  \\leftarrow Limit Cycle\n  \\omega_{lc} = %.4f rad/s (%.4f Hz)\n' ...
         '  A_{lc} = %.4f,  \\varpi = %.4f'], wk, wk/(2*pi), Ak, vk), ...
        'Color','r','FontWeight','bold','FontSize',9);
end

if n_lc == 0
    fprintf('%s\n','   -- G(jw) never reaches the -1/N line in the lower half-plane --');
    text(ax, -1.5, 2, 'NO LIMIT CYCLE!', 'Color','r','FontSize',20,'FontWeight','bold', ...
         'HorizontalAlignment','center','BackgroundColor','w', ...
         'EdgeColor','r','LineWidth',1.5,'Margin',6);
end
fprintf('===========================================================================\n');

legend(ax, {'G(j\omega)','-1/N locus (Re = -\pi^2/8)'}, 'Location','southwest');

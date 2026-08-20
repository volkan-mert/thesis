%% RLE DESCRIBING FUNCTION AND LIMIT CYCLE ANALYSIS
%
% Limit-cycle condition:
%
%       1 + G(jw)*N(A,w) = 0
%
% where
%
%       G(s) = Kp*Gc(s)*Gac(s)
%
% Equivalent forms:
%
%       G(jw) = -1/N(A,w)
%
%       G(jw)*N(A,w) = -1
%

clear;
clc;
close all;


%% ==============================================================
%  1. TRANSFER FUNCTION OF THE CONTROLLER
% ==============================================================

numGc = [5.21, -273.7855, -1425.2, -700.1952];

denGc = [1, 21.3594, 545.5538, 605.6621, 0];

Gc = tf(numGc,denGc)


%% ==============================================================
%  2. TRANSFER FUNCTION OF THE AIRCRAFT LONGITUDINAL DYNAMICS
% ==============================================================

numGac = [-10.5240, -16.8384, -0.6247, 0];

denGac = [1, 2.3473, -5.3061, -0.1836, -0.0418];

Gac = tf(numGac,denGac)


%% ==============================================================
%  3. FORWARD GAIN
% ==============================================================

Kp = 13.68;


%% ==============================================================
%  4. OPEN-LOOP AND CLOSED-LOOP SYSTEM
% ==============================================================

G = Kp*Gc*Gac;

sysCL = feedback(G,1);

sysCL_minimal = minreal(sysCL);


%% ==============================================================
%  5. CLOSED-LOOP STABILITY CHECK
% ==============================================================

poles_CL = pole(sysCL_minimal);


fprintf('\n');
fprintf('============================================================\n');
fprintf('CLOSED-LOOP STABILITY CHECK\n');
fprintf('============================================================\n\n');

fprintf('Poles of sysCL_minimal:\n\n');

disp(poles_CL);


if all(real(poles_CL) < 0)

    fprintf('All poles have negative real parts.\n');
    fprintf('sysCL IS STABLE.\n');

else

    fprintf('At least one pole has zero or positive real part.\n');
    fprintf('sysCL IS NOT STABLE.\n');

end


%% ==============================================================
%  6. RATE LIMITER PARAMETERS
% ==============================================================

R = 60;          % Maximum actuator rate, deg/s

Ai = 13.68;      % Initial input amplitude, deg


w_onset = R/Ai;

w_crit = 1.862*w_onset;


fprintf('\n');
fprintf('============================================================\n');
fprintf('RATE LIMITER PARAMETERS\n');
fprintf('============================================================\n\n');

fprintf('R          = %.2f deg/s\n',R);
fprintf('Ai         = %.2f deg\n',Ai);
fprintf('w_onset    = %.4f rad/s\n',w_onset);
fprintf('w_crit     = %.4f rad/s\n',w_crit);


%% ==============================================================
%  7. CHECK Kp*Gc*Gac*N FOR FIXED Ai = 13.68 deg
% ==============================================================

w_fixed = logspace(-1,2,4000);


G_fixed = squeeze(freqresp(G,w_fixed));

G_fixed = G_fixed(:).';


N_fixed = rle_df(Ai,w_fixed,R);


error_fixed = abs(1 + G_fixed.*N_fixed);


[min_error_fixed,index_fixed] = min(error_fixed);

w_best_fixed = w_fixed(index_fixed);


fprintf('\n');
fprintf('============================================================\n');
fprintf('FIXED AMPLITUDE LIMIT-CYCLE CHECK\n');
fprintf('============================================================\n\n');

fprintf('Kp                 = %.2f\n',Kp);
fprintf('Ai                 = %.2f deg\n',Ai);
fprintf('Closest frequency  = %.4f rad/s\n',w_best_fixed);
fprintf('Minimum error      = %.6e\n',min_error_fixed);


if min_error_fixed < 1e-3

    fprintf('\nLIMIT CYCLE FOUND FOR Ai = %.2f deg\n',Ai);

else

    fprintf('\nNO LIMIT CYCLE FOUND FOR Ai = %.2f deg\n',Ai);

end


%% ==============================================================
%  8. Kp SWEEP SETTINGS
% ==============================================================

Kp_values = 0.1:0.1:100;


A_min = 1;

A_max = 100;


w_min = 0.1;

w_max = 100;


G_base = Gc*Gac;


fprintf('\n');
fprintf('============================================================\n');
fprintf('LIMIT CYCLE SWEEP SETTINGS\n');
fprintf('============================================================\n\n');

fprintf('Kp range = 0.1 to 100 with 0.1 increments\n');
fprintf('A range  = %.1f to %.1f deg\n',A_min,A_max);
fprintf('w range  = %.1f to %.1f rad/s\n',w_min,w_max);


%% ==============================================================
%  9. SEARCH REGION II
%
%  Region II:
%
%       1 < alpha < 1.862
%
%       alpha = w*A/R
%
%  Therefore:
%
%       w_onset < w < w_crit
% ==============================================================

fprintf('\n');
fprintf('Searching Region II limit cycles.\n');


Region2 = find_region_cycles(G_base,R,Kp_values,A_min,A_max,w_min,w_max,2);


%% ==============================================================
%  10. SEARCH REGION III
%
%  Region III:
%
%       alpha >= 1.862
%
%  Therefore:
%
%       w >= w_crit
% ==============================================================

fprintf('Searching Region III limit cycles.\n');


Region3 = find_region_cycles(G_base,R,Kp_values,A_min,A_max,w_min,w_max,3);


%% ==============================================================
%  11. COMBINE REGION II AND REGION III
% ==============================================================

LimitCycleResults = [Region2; Region3];


if ~isempty(LimitCycleResults)

    LimitCycleResults = sortrows(LimitCycleResults,{'Kp','Region'});

end


%% ==============================================================
%  12. CHECK REGION CONDITIONS
% ==============================================================

RegionCheck = strings(height(LimitCycleResults),1);


for i = 1:height(LimitCycleResults)

    w_now = LimitCycleResults.w(i);

    w_onset_now = LimitCycleResults.w_onset(i);

    w_crit_now = LimitCycleResults.w_crit(i);


    if LimitCycleResults.Region(i) == "Region II"

        if w_now > w_onset_now && w_now < w_crit_now

            RegionCheck(i) = "OK";

        else

            RegionCheck(i) = "FAILED";

        end

    end


    if LimitCycleResults.Region(i) == "Region III"

        if w_now >= w_crit_now

            RegionCheck(i) = "OK";

        else

            RegionCheck(i) = "FAILED";

        end

    end

end


%% ==============================================================
%  13. DISPLAY ALL LIMIT CYCLES AS ONE TABLE
% ==============================================================

fprintf('\n');
fprintf('=============================================================================\n');
fprintf('              REGION II AND REGION III LIMIT CYCLE RESULTS\n');
fprintf('=============================================================================\n\n');


if isempty(LimitCycleResults)

    fprintf('NO LIMIT CYCLES FOUND.\n');

else

    LimitCycleTable = LimitCycleResults(:,{'Region','Kp','A','w','w_onset','w_crit'});

    LimitCycleTable.RegionCheck = RegionCheck;

    format short g

    disp(LimitCycleTable);

    fprintf('Total number of limit cycles = %d\n',height(LimitCycleTable));
    fprintf('Number of Region II limit cycles = %d\n',sum(LimitCycleTable.Region == "Region II"));
    fprintf('Number of Region III limit cycles = %d\n',sum(LimitCycleTable.Region == "Region III"));

end


%% ==============================================================
%  14. LIMIT-CYCLE AMPLITUDE AGAINST Kp
% ==============================================================

if ~isempty(LimitCycleResults)

    index_R2 = LimitCycleResults.Region == "Region II";

    index_R3 = LimitCycleResults.Region == "Region III";


    figure;

    hold on;


    if any(index_R2)

        plot(LimitCycleResults.Kp(index_R2),LimitCycleResults.A(index_R2),'o-','LineWidth',1.5,'MarkerSize',7);

    end


    if any(index_R3)

        plot(LimitCycleResults.Kp(index_R3),LimitCycleResults.A(index_R3),'s-','LineWidth',1.5,'MarkerSize',7);

    end


    grid on;

    xlabel('K_p');

    ylabel('Limit Cycle Amplitude A (deg)');

    title('Limit Cycle Amplitude vs K_p');

    legend('Region II','Region III','Location','best');

end


%% ==============================================================
%  15. LIMIT-CYCLE FREQUENCY AGAINST Kp
% ==============================================================

if ~isempty(LimitCycleResults)

    figure;

    hold on;


    if any(index_R2)

        plot(LimitCycleResults.Kp(index_R2),LimitCycleResults.w(index_R2),'o-','LineWidth',1.5,'MarkerSize',7);

    end


    if any(index_R3)

        plot(LimitCycleResults.Kp(index_R3),LimitCycleResults.w(index_R3),'s-','LineWidth',1.5,'MarkerSize',7);

    end


    grid on;

    xlabel('K_p');

    ylabel('\omega (rad/s)');

    title('Limit Cycle Frequency vs K_p');

    legend('Region II','Region III','Location','best');

end


%% ==============================================================
%  16. w, w_onset AND w_crit AGAINST Kp
% ==============================================================

if ~isempty(LimitCycleResults)

    figure;

    hold on;


    plot(LimitCycleResults.Kp,LimitCycleResults.w,'o-','LineWidth',1.5);

    plot(LimitCycleResults.Kp,LimitCycleResults.w_onset,'--','LineWidth',1.5);

    plot(LimitCycleResults.Kp,LimitCycleResults.w_crit,'-.','LineWidth',1.5);


    grid on;

    xlabel('K_p');

    ylabel('Frequency (rad/s)');

    title('\omega, \omega_{onset} and \omega_{crit}');

    legend('\omega','\omega_{onset}','\omega_{crit}','Location','best');

end


%% ==============================================================
%  17. NYQUIST AND NICHOLS PLOTS
%
%  For every limit cycle:
%
%  Figure 1:
%
%       G(jw) versus -1/N(A,w)
%
%  Figure 2:
%
%       G(jw)*N(A,w)
%
%  Figure 3:
%
%       Nichols chart of G(jw) versus -1/N(A,w)
% ==============================================================

if ~isempty(LimitCycleResults)

    w_plot = logspace(log10(w_min),log10(w_max),4000);


    for i = 1:height(LimitCycleResults)


        %% Current limit-cycle values

        Region_now = LimitCycleResults.Region(i);

        Region_char = char(Region_now);

        Kp_now = LimitCycleResults.Kp(i);

        A_now = LimitCycleResults.A(i);

        w_now = LimitCycleResults.w(i);

        w_onset_now = LimitCycleResults.w_onset(i);

        w_crit_now = LimitCycleResults.w_crit(i);


        %% Linear system for current Kp

        G_now = Kp_now*Gc*Gac;


        %% Frequency response of G(jw)

        G_response_plot = squeeze(freqresp(G_now,w_plot));

        G_response_plot = G_response_plot(:).';


        %% Describing function N(A,w)

        N_response_plot = rle_df(A_now,w_plot,R);


        %% -1/N(A,w)

        minus_inv_N_plot = -1./N_response_plot;


        %% G(jw)*N(A,w)

        GN_response_plot = G_response_plot.*N_response_plot;


        %% Convert responses to FRD models

        G_frd = frd(reshape(G_response_plot,1,1,length(w_plot)),w_plot);

        minus_inv_N_frd = frd(reshape(minus_inv_N_plot,1,1,length(w_plot)),w_plot);

        GN_frd = frd(reshape(GN_response_plot,1,1,length(w_plot)),w_plot);


        %% Values exactly at the limit-cycle frequency

        G_LC = squeeze(freqresp(G_now,w_now));

        N_LC = rle_df(A_now,w_now,R);

        minus_inv_N_LC = -1/N_LC;

        GN_LC = G_LC*N_LC;


        %% ========================================================
        %  FIGURE 1
        %
        %  NYQUIST:
        %
        %       G(jw) vs -1/N(A,w)
        %
        %  Limit-cycle condition:
        %
        %       G(jw) = -1/N(A,w)
        % ========================================================

        figure('Name',sprintf('Nyquist G vs -1/N | Kp = %.1f',Kp_now),'NumberTitle','off');


        nyquist(G_frd,minus_inv_N_frd);


        grid on;

        hold on;


        plot(real(G_LC),imag(G_LC),'ko','MarkerFaceColor','k','MarkerSize',7);


        plot(real(minus_inv_N_LC),imag(minus_inv_N_LC),'kx','MarkerSize',10,'LineWidth',1.5);


        title(sprintf('Nyquist: G(j\\omega) vs -1/N(A,\\omega) | %s | K_p = %.1f | A = %.2f deg | \\omega = %.3f rad/s',Region_char,Kp_now,A_now,w_now));


        text(real(G_LC),imag(G_LC),sprintf('  LC: K_p = %.1f, A = %.2f deg, \\omega = %.3f rad/s',Kp_now,A_now,w_now),'FontWeight','bold');


        legend('G(j\omega)','-1/N(A,\omega)','G(j\omega_{LC})','-1/N(A,\omega_{LC})','Location','best');


        %% ========================================================
        %  FIGURE 2
        %
        %  NYQUIST:
        %
        %       G(jw)*N(A,w)
        %
        %  Limit-cycle condition:
        %
        %       G(jw)*N(A,w) = -1
        % ========================================================

        figure('Name',sprintf('Nyquist G*N | Kp = %.1f',Kp_now),'NumberTitle','off');


        nyquist(GN_frd);


        grid on;

        hold on;


        plot(-1,0,'rx','MarkerSize',12,'LineWidth',2);


        plot(real(GN_LC),imag(GN_LC),'ko','MarkerFaceColor','k','MarkerSize',7);


        title(sprintf('Nyquist: G(j\\omega)N(A,\\omega) | %s | K_p = %.1f | A = %.2f deg | \\omega = %.3f rad/s',Region_char,Kp_now,A_now,w_now));


        text(-0.95,0.15,sprintf('LC: K_p = %.1f, A = %.2f deg, \\omega = %.3f rad/s\n\\omega_{onset} = %.3f rad/s\n\\omega_{crit} = %.3f rad/s',Kp_now,A_now,w_now,w_onset_now,w_crit_now),'FontWeight','bold');


        legend('G(j\omega)N(A,\omega)','Critical Point (-1,0)','Limit Cycle','Location','best');


        %% ========================================================
        %  FIGURE 3
        %
        %  NICHOLS:
        %
        %       G(jw) vs -1/N(A,w)
        % ========================================================

        figure('Name',sprintf('Nichols G vs -1/N | Kp = %.1f',Kp_now),'NumberTitle','off');


        nicholsplot(G_frd,minus_inv_N_frd);


        grid on;


        title(sprintf('Nichols: G(j\\omega) vs -1/N(A,\\omega) | %s | K_p = %.1f | A = %.2f deg | \\omega = %.3f rad/s',Region_char,Kp_now,A_now,w_now));


        legend('G(j\omega)','-1/N(A,\omega)','Location','best');


    end

end


%% ==============================================================
%  LOCAL FUNCTION
%  FIND REGION II OR REGION III LIMIT CYCLES
% ==============================================================

function ResultTable = find_region_cycles(G_base,R,Kp_values,A_min,A_max,w_min,w_max,region_number)


%% Select alpha range

if region_number == 2

    alpha_values = linspace(1.0001,1.8619,1200);

    RegionName = "Region II";

else

    alpha_max = A_max*w_max/R;

    alpha_part1 = linspace(1.862,min(5,alpha_max),2500);


    if alpha_max > 5

        alpha_part2 = logspace(log10(5.01),log10(alpha_max),800);

        alpha_values = unique([alpha_part1 alpha_part2]);

    else

        alpha_values = alpha_part1;

    end


    RegionName = "Region III";

end


%% Frequency vector

w_values = logspace(log10(w_min),log10(w_max),5000);


%% Frequency response of Gc*Gac

G_response = squeeze(freqresp(G_base,w_values));

G_response = G_response(:).';


%% Storage for solution branches

max_branches = 6;

Kp_branch = NaN(length(alpha_values),max_branches);

w_branch = NaN(length(alpha_values),max_branches);


%% ==============================================================
%  FIND NEGATIVE REAL-AXIS CROSSINGS
% ==============================================================

for ia = 1:length(alpha_values)

    alpha = alpha_values(ia);

    N_alpha = rle_df_alpha(alpha);

    H = G_response*N_alpha;

    imag_H = imag(H);


    crossing_index = find(imag_H(1:end-1).*imag_H(2:end) <= 0);


    branch_number = 0;


    for j = 1:length(crossing_index)

        k = crossing_index(j);

        w1 = w_values(k);

        w2 = w_values(k+1);


        try

            w_zero = fzero(@(ww) imag(squeeze(freqresp(G_base,ww))*N_alpha),[w1 w2]);

        catch

            continue;

        end


        H_zero = squeeze(freqresp(G_base,w_zero))*N_alpha;


        if real(H_zero) < 0

            Kp_required = -1/real(H_zero);

            A_required = alpha*R/w_zero;


            if Kp_required >= min(Kp_values) && Kp_required <= max(Kp_values) && A_required >= A_min && A_required <= A_max && w_zero >= w_min && w_zero <= w_max

                branch_number = branch_number + 1;


                if branch_number <= max_branches

                    Kp_branch(ia,branch_number) = Kp_required;

                    w_branch(ia,branch_number) = w_zero;

                end

            end

        end

    end

end


%% ==============================================================
%  SEARCH THE DISCRETE Kp GRID
% ==============================================================

Kp_result = [];

A_result = [];

w_result = [];

w_onset_result = [];

w_crit_result = [];

error_result = [];


for Kp_test = Kp_values

    for branch = 1:max_branches

        K_curve = Kp_branch(:,branch);

        W_curve = w_branch(:,branch);


        for ia = 1:length(alpha_values)-1

            if isnan(K_curve(ia)) || isnan(K_curve(ia+1))

                continue;

            end


            d1 = K_curve(ia) - Kp_test;

            d2 = K_curve(ia+1) - Kp_test;


            if d1*d2 <= 0

                if abs(K_curve(ia+1) - K_curve(ia)) > 1e-12

                    alpha0 = interp1([K_curve(ia) K_curve(ia+1)],[alpha_values(ia) alpha_values(ia+1)],Kp_test);

                else

                    alpha0 = (alpha_values(ia) + alpha_values(ia+1))/2;

                end


                w0 = interp1([alpha_values(ia) alpha_values(ia+1)],[W_curve(ia) W_curve(ia+1)],alpha0);


                if ~isfinite(w0) || w0 <= 0

                    continue;

                end


                x0 = [log(w0) alpha0];


                options = optimset('Display','off','TolX',1e-10,'TolFun',1e-12);


                x = fminsearch(@(x) region_error(x,Kp_test,G_base,R,region_number,A_min,A_max,w_min,w_max),x0,options);


                w_final = exp(x(1));

                alpha_final = x(2);

                A_final = alpha_final*R/w_final;


                if region_number == 2

                    region_ok = alpha_final > 1 && alpha_final < 1.862;

                else

                    region_ok = alpha_final >= 1.862;

                end


                limits_ok = A_final >= A_min && A_final <= A_max && w_final >= w_min && w_final <= w_max;


                if region_ok && limits_ok

                    N_final = rle_df(A_final,w_final,R);

                    G_final = squeeze(freqresp(Kp_test*G_base,w_final));

                    LC_error = abs(1 + G_final*N_final);


                    if LC_error < 1e-5

                        duplicate = false;


                        for m = 1:length(Kp_result)

                            same_Kp = abs(Kp_result(m) - Kp_test) < 1e-8;

                            same_A = abs(A_result(m) - A_final) < 1e-3;

                            same_w = abs(w_result(m) - w_final) < 1e-4;


                            if same_Kp && same_A && same_w

                                duplicate = true;

                            end

                        end


                        if ~duplicate

                            w_onset_final = R/A_final;

                            w_crit_final = 1.862*w_onset_final;


                            Kp_result(end+1,1) = Kp_test;

                            A_result(end+1,1) = A_final;

                            w_result(end+1,1) = w_final;

                            w_onset_result(end+1,1) = w_onset_final;

                            w_crit_result(end+1,1) = w_crit_final;

                            error_result(end+1,1) = LC_error;

                        end

                    end

                end

            end

        end

    end

end


%% ==============================================================
%  CREATE RESULT TABLE
% ==============================================================

if isempty(Kp_result)

    Region = strings(0,1);

    Kp_result = zeros(0,1);

    A_result = zeros(0,1);

    w_result = zeros(0,1);

    w_onset_result = zeros(0,1);

    w_crit_result = zeros(0,1);

    error_result = zeros(0,1);


    ResultTable = table(Region,Kp_result,A_result,w_result,w_onset_result,w_crit_result,error_result,'VariableNames',{'Region','Kp','A','w','w_onset','w_crit','Error'});

else

    Region = repmat(RegionName,length(Kp_result),1);


    ResultTable = table(Region,Kp_result,A_result,w_result,w_onset_result,w_crit_result,error_result,'VariableNames',{'Region','Kp','A','w','w_onset','w_crit','Error'});


    ResultTable = sortrows(ResultTable,'Kp');

end


end


%% ==============================================================
%  LOCAL FUNCTION
%  REGION II / REGION III LIMIT-CYCLE ERROR
% ==============================================================

function error = region_error(x,Kp,G_base,R,region_number,A_min,A_max,w_min,w_max)


w = exp(x(1));

alpha = x(2);


if ~isfinite(w) || ~isfinite(alpha)

    error = 1e12;

    return;

end


%% Check region

if region_number == 2

    if alpha <= 1 || alpha >= 1.862

        error = 1e6;

        return;

    end

else

    if alpha < 1.862

        error = 1e6;

        return;

    end

end


%% Calculate amplitude

A = alpha*R/w;


%% Check limits

if A < A_min || A > A_max || w < w_min || w > w_max

    error = 1e6;

    return;

end


%% Describing function

N = rle_df(A,w,R);


%% Frequency response

Gjw = squeeze(freqresp(Kp*G_base,w));


%% Limit-cycle condition

error = abs(1 + Gjw*N)^2;


end


%% ==============================================================
%  LOCAL FUNCTION
%  RLE DESCRIBING FUNCTION
% ==============================================================

function N = rle_df(A,w,R)


%% Onset frequency

w_onset = R/A;


%% Normalized frequency

alpha = w/w_onset;


%% Initialize magnitude and phase

M = ones(size(w));

phi = zeros(size(w));


%% Define regions

region1 = alpha <= 1;

region2 = alpha > 1 & alpha < 1.862;

region3 = alpha >= 1.862;


%% Region I

M(region1) = 1;

phi(region1) = 0;


%% Region II

a = alpha(region2);


M(region2) = 0.2908*a.^3 - 1.4396*a.^2 + 1.9232*a + 0.223;

phi(region2) = 0.528*a.^3 - 2.6213*a.^2 + 3.5056*a - 1.4171;


%% Region III

varpi = w_onset./w(region3);


M(region3) = (4/pi).*varpi;

phi(region3) = -acos((pi/2).*varpi);


%% Complex describing function

N = M.*exp(1j*phi);


end


%% ==============================================================
%  LOCAL FUNCTION
%  RLE DESCRIBING FUNCTION FROM alpha
% ==============================================================

function N = rle_df_alpha(alpha)


%% Region I

if alpha <= 1

    M = 1;

    phi = 0;


%% Region II

elseif alpha < 1.862

    M = 0.2908*alpha^3 - 1.4396*alpha^2 + 1.9232*alpha + 0.223;

    phi = 0.528*alpha^3 - 2.6213*alpha^2 + 3.5056*alpha - 1.4171;


%% Region III

else

    varpi = 1/alpha;

    M = (4/pi)*varpi;

    phi = -acos((pi/2)*varpi);

end


%% Complex describing function

N = M*exp(1j*phi);


end
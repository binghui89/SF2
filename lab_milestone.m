function lab_milestone()
% test_logn_distribution();
% test_discretize_lognormal();
% ramping_requirement_rtd_nl();
ramping_requirement_rtpd_nl();
% cell_history = return_30_days(2019, 5, 31);
end

function cell_history= return_historical_days(YYYY, MM, DD, ndays)
if nargin==3
    ndays = 30;
end
thisday = datetime(YYYY, MM, DD);
if mod(weekday(thisday), 6) == 1
    todayisweekday = false;
else
    todayisweekday = true;
end
cell_history = cell(ndays, 1);

validnumber = 0;
historyday  = thisday;
while validnumber < ndays
    historyday = historyday - days(1);
    if mod(weekday(historyday), 6) == 1
        historyisweekday = false;
    else
        historyisweekday = true;
    end
    if historyisweekday == todayisweekday
        validnumber = validnumber + 1;
        cell_history{validnumber} = historyday;
    end
end

end

function ramping_requirement_rtpd_nl(write_flag)
if nargin == 0
    write_flag = false;
end
oasis = readtable('for_flexiramp_summary.csv');

M_b = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghuinew\binghui_10min.csv', 1, 0); % Binding forecast
M_a = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghuinew\binghui_40min.csv', 1, 0); % Advisory forecast


% Replace 0 with nan
M_b(M_b(:, 6) ==0, 6) = nan;
M_b(M_b(:, 7) ==0, 7) = nan;
M_a(M_a(:, 6) ==0, 6) = nan;
M_a(M_a(:, 7) ==0, 7) = nan;

% Show the actual NL and the deterministically forecasted NL
figure();
subplot(2, 1, 1);
plot(datetime(M_b(:, 1), M_b(:, 2), M_b(:, 3), M_b(:, 4), M_b(:, 5), 0), M_b(:, 7), datetime(M_b(:, 1), M_b(:, 2), M_b(:, 3), M_b(:, 4), M_b(:, 5), 0), M_b(:, 6))
legend('Advisory RTPD', 'Actual 15-min mean');
title('Advisory RTPD');
subplot(2, 1, 2);
plot(datetime(M_a(:, 1), M_a(:, 2), M_a(:, 3), M_a(:, 4), M_a(:, 5), 0), M_a(:, 7), datetime(M_a(:, 1), M_a(:, 2), M_a(:, 3), M_a(:, 4), M_a(:, 5), 0), M_a(:, 6))
legend('Binding RTD', 'Actual 5-min');
title('Binding RTD');

% Since M_a is 15min, we need to repeat it three times.
M_a5 = nan(size(M_b));
M_a5(:, 1:5) = M_b(:, 1:5);
for i = 0:5:55
    rows = (M_a(:, 5) == (i-mod(i, 15)));
    M_a5(M_a5(:, 5)==i, 6:end) = M_a(rows, 6:end);
end
M_a = M_a5;

Mquantiles_b = M_b(:, 8:end);
Mquantiles_a = M_a(:, 8:end);
nlactual_b = M_b(:, 6);
nlactual_a = M_a(:, 6);
nldeterm_b = M_b(:, 7);
nldeterm_a = M_a(:, 7);

p = 1:99; % 1 to 99 quantiles
binsize = 100;
select_month = 5;

for select_day = 27: 31
    
    % Let's calculate the baseline requirement first
    % First, collect historical NL forecast errors
    cell_history = return_historical_days(2019, select_month, select_day, 40);
    ndays = numel(cell_history);
    nlerr_rtpd_max = nan(96, ndays);
    nlerr_rtpd_min = nan(96, ndays);
    nlerr_rtpd = nan(96, ndays, 3); % One 15-min interval includes three 5-min
    for i = 1: ndays
        for j = 0: 2
            frcst_a_rtpd = M_a((M_a(:, 2) ==cell_history{i}.Month)&(M_a(:, 3) ==cell_history{i}.Day)&(mod(M_a(:, 5), 15)==5*j), 7);
            frcst_b_rtd  = M_b((M_b(:, 2) ==cell_history{i}.Month)&(M_b(:, 3) ==cell_history{i}.Day)&(mod(M_b(:, 5), 15)==5*j), 7);
            nlerr_rtpd(:, i, j+1) = frcst_b_rtd - frcst_a_rtpd;
        end
    end
    nlerr_rtpd_max = max(nlerr_rtpd, [], 3);
    nlerr_rtpd_min = min(nlerr_rtpd, [], 3);

    % Next, calculate ramping reserve requirements
    fru_determ = nan(24, 1);
    frd_determ = nan(24, 1);
    for i = 1: 24
        rtpd_col = (i-1)*4+1: i*4;
        samples_fru = nlerr_rtpd_max(rtpd_col, :);
        [f,x] = ecdf(samples_fru(:));
        fru_determ(i) = interp1(f, x, 0.975);

        samples_frd = nlerr_rtpd_min(rtpd_col, :);
        [f,x] = ecdf(samples_frd(:));
        frd_determ(i) = interp1(f, x, 0.025);
    end

    % Now, let's calculate probabilistic requirements
    irows = find((M_b(:, 2) ==select_month)&(M_b(:, 3) == select_day));
    frd_prob = nan(size(irows, 1), 1);
    fru_prob = nan(size(irows, 1), 1);
    for i = irows'
        Mrow_b = Mquantiles_b(i, :);
        Mrow_a = Mquantiles_a(i, :);
        [h_bin_b, binedge_b] = discretize_lognormal(Mrow_b, p, binsize);
        [h_bin_a, binedge_a] = discretize_lognormal(Mrow_a, p, binsize);
        [a_convoluted, t_convoluted] = conv_poly([h_bin_a(:); 0], binedge_a(:), [flipud(h_bin_b(:)); 0], flipud(-binedge_b(:)), binsize); % Distribution of binding - advisory forecast NL

        frd_prob(irows==i) = interp1(cdf_poly(a_convoluted, t_convoluted), t_convoluted, 0.025);  % 5 percentile
        fru_prob(irows==i) = interp1(cdf_poly(a_convoluted, t_convoluted), t_convoluted, 0.975);  % 95 percentile

%         if mod(i, 12) == 1
%             fig = figure();
%             subplot(3, 1, 1);
%             plot_conv_poly(fig, [h_bin_b(:); 0], binedge_b, 'b');
%             subplot(3, 1, 2);
%             plot_conv_poly(fig, [h_bin_a(:); 0], binedge_a, 'b');
%             subplot(3, 1, 3);
%             plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
%         end
    end

    figure();
    hold on;
    xtime = datetime(2019,select_month,select_day,0,5,0):minutes(5): (datetime(2019,select_month,select_day,0,5,0) + days(1) - minutes(5));

    % h1 = plot(xtime, reshape(repmat(min(reshape(frd_prob, 12, 24)), 12, 1), 288, 1), xtime, reshape(repmat(max(reshape(fru_prob, 12, 24)), 12, 1), 288, 1));
    h1 =plot(xtime, reshape(repmat(max(reshape(frd_prob, 12, 24)), 12, 1), 288, 1), xtime, reshape(repmat(min(reshape(fru_prob, 12, 24)), 12, 1), 288, 1));
%     h1 = plot(xtime, reshape(repmat(mean(reshape(frd_prob, 12, 24), 1), 12, 1), 288, 1), xtime, reshape(repmat(mean(reshape(fru_prob, 12, 24), 1), 12, 1), 288, 1));
    set(h1, 'linewidth', 2, 'color', 'k');

    h2 = plot(xtime, reshape(repmat(fru_determ, 1, 12)', 288, 1), xtime, reshape(repmat(frd_determ, 1, 12)', 288, 1));
    set(h2, 'linewidth', 2, 'color', 'r');

%     fru_oasis = oasis.UP_RTPD((oasis.OPR_DT.Year==2019) & (oasis.OPR_DT.Month==select_month) & (oasis.OPR_DT.Day==select_day));
%     frd_oasis = oasis.DOWN_RTPD((oasis.OPR_DT.Year==2019) & (oasis.OPR_DT.Month==select_month) & (oasis.OPR_DT.Day==select_day));
%     h3 = plot(xtime, fru_oasis, xtime, -frd_oasis);
%     set(h3, 'linewidth', 2,  'color', 'b');

    ylabel('MW');
    legend([h1(1); h2(1)], 'Prob', 'Baseline');
%     legend([h1(1); h2(1); h3(1)], 'FRD prob', 'FRD determ', 'FRP OASIS');
%     legend([h2(1); h3(1)], 'Baseline', 'OASIS');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    legend boxoff;
    box on;
    
    figure();
    hold on;
    frd_floor24 = min(max(reshape(frd_prob, 12, 24)), 0);
    fru_floor24 = max(min(reshape(fru_prob, 12, 24)), 0);
    h1 = plot(xtime, reshape(repmat(frd_floor24, 12, 1), 288, 1)./reshape(repmat(frd_determ, 1, 12)', 288, 1));
    h2 = plot(xtime, reshape(repmat(fru_floor24, 12, 1), 288, 1)./reshape(repmat(fru_determ, 1, 12)', 288, 1));
    h3 = plot(xtime, 0.9.*ones(288, 1), 'r');
    set([h1,h2], 'linewidth', 2);
    set(gca, 'YScale', 'log')
    legend([h1, h2], 'FRD', 'FRU');
    ylabel('Ratio of Prob/Baseline');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    legend boxoff;
    box on;
    
%     ratio_up = min(max(reshape(frd_prob, 12, 24)), 0)./frd_determ';
%     ratio_dn = max(min(reshape(fru_prob, 12, 24)), 0)./fru_determ';
    
    % Save results
    if write_flag
        M = [... 
            xtime(:).Year, ...
            xtime(:).Month, ...
            xtime(:).Day, ...
            xtime(:).Hour, ...
            xtime(:).Minute, ...
            xtime(:).Second, ...
            fru_oasis, ...
            -frd_oasis, ...
            reshape(repmat(fru_determ, 1, 12)', 288, 1), ...
            reshape(repmat(frd_determ, 1, 12)', 288, 1), ... 
            reshape(repmat(min(reshape(fru_prob, 12, 24)), 12, 1), 288, 1), ... 
            reshape(repmat(max(reshape(frd_prob, 12, 24)), 12, 1), 288, 1) ...
            ];

        csvname_write = strcat(int2str(select_day), '.csv');
        dirhome = pwd;
        cHeader = {'Year' 'Month' 'Day' 'Hour' 'Minute' 'Second' 'FRU_OASIS' 'FRD_OASIS' 'FRU_BASELINE' 'FRD_BASELINE' 'FRU_PROB' 'FRD_PROB'}; %dummy header
        commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
        commaHeader = commaHeader(:)';
        textHeader = cell2mat(commaHeader); % cHeader in text with commas

        cd('FRP_RTPD');
        fid = fopen(csvname_write,'w'); 
        fprintf(fid,'%s\n',textHeader); % write header to file
        fclose(fid);
        dlmwrite(csvname_write, M, '-append');
        cd(dirhome);
    end

end

end

function ramping_requirement_rtd_nl()
oasis = readtable('for_flexiramp_summary.csv');

M_b = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghuinew\binghui_10min.csv', 1, 0); % Binding forecast
M_a = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghuinew\binghui_15min.csv', 1, 0); % Advisory forecast


% Replace 0 with nan
M_b(M_b(:, 6) ==0, 6) = nan;
M_b(M_b(:, 7) ==0, 7) = nan;
M_a(M_a(:, 6) ==0, 6) = nan;
M_a(M_a(:, 7) ==0, 7) = nan;

% Show the actual NL and the deterministically forecasted NL
figure();
ax1 = subplot(2, 1, 1);
xtime = datetime(M_b(:, 1), M_b(:, 2), M_b(:, 3), M_b(:, 4), M_b(:, 5), 0);
plot(xtime, M_b(:, 7), xtime, M_b(:, 6), xtime, M_a(:, 7))
legend('Binding RTD', 'Actual 5-min', 'Advisory RTD');
ax2 = subplot(2, 1, 2);
plot(xtime, M_b(:, 7)-M_a(:, 7))
title('B - A');
linkaxes([ax1, ax2],'x');

Mquantiles_b = M_b(:, 8:end);
Mquantiles_a = M_a(:, 8:end);
nlactual_b = M_b(:, 6);
nlactual_a = M_a(:, 6);
nldeterm_b = M_b(:, 7);
nldeterm_a = M_a(:, 7);

p = 1:99; % 1 to 99 quantiles
binsize = 100;
select_month = 5;

for select_day = 27: 31
    
    % Let's calculate the baseline requirement first
    % First, collect historical NL forecast errors
    cell_history = return_historical_days(2019, select_month, select_day);
    ndays = numel(cell_history);
    nlerr_rtpd = nan(288, ndays); % One 15-min interval includes three 5-min
    for i = 1: ndays
        frcst_a_rtd = M_a((M_a(:, 2) ==cell_history{i}.Month)&(M_a(:, 3) ==cell_history{i}.Day), 7);
        frcst_b_rtd = M_b((M_b(:, 2) ==cell_history{i}.Month)&(M_b(:, 3) ==cell_history{i}.Day), 7);
        nlerr_rtpd(:, i) = frcst_b_rtd - frcst_a_rtd;
    end

    % Next, calculate ramping reserve requirements
    fru_determ = nan(24, 1);
    frd_determ = nan(24, 1);
    for i = 1: 24
        rtpd_col = (i-1)*12+1: i*12;
        samples_fru = nlerr_rtpd(rtpd_col, :);
        [f,x] = ecdf(samples_fru(:));
        fru_determ(i) = interp1(f, x, 0.975);
        frd_determ(i) = interp1(f, x, 0.025);
    end

    % Now, let's calculate probabilistic requirements
    irows = find((M_b(:, 2) ==select_month)&(M_b(:, 3) == select_day));
    frd_prob = nan(size(irows, 1), 1);
    fru_prob = nan(size(irows, 1), 1);
    for i = irows'
        Mrow_b = Mquantiles_b(i, :);
        Mrow_a = Mquantiles_a(i, :);
        [h_bin_b, binedge_b] = discretize_lognormal(Mrow_b, p, binsize);
        [h_bin_a, binedge_a] = discretize_lognormal(Mrow_a, p, binsize);
        [a_convoluted, t_convoluted] = conv_poly([h_bin_a(:); 0], binedge_a(:), [flipud(h_bin_b(:)); 0], flipud(-binedge_b(:)), binsize); % Distribution of binding - advisory forecast NL

        frd_prob(irows==i) = interp1(cdf_poly(a_convoluted, t_convoluted), t_convoluted, 0.025);  % 5 percentile
        fru_prob(irows==i) = interp1(cdf_poly(a_convoluted, t_convoluted), t_convoluted, 0.975);  % 95 percentile

%         if mod(i, 12) == 1
%             fig = figure();
%             subplot(3, 1, 1);
%             plot_conv_poly(fig, [h_bin_b(:); 0], binedge_b, 'b');
%             subplot(3, 1, 2);
%             plot_conv_poly(fig, [h_bin_a(:); 0], binedge_a, 'b');
%             subplot(3, 1, 3);
%             plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
%         end
    end

    figure(); 
    hold on;
    xtime = datetime(2019,select_month,select_day,0,5,0):minutes(5): (datetime(2019,select_month,select_day,0,5,0) + day(1) - minutes(5));

    % h1 = plot(xtime, reshape(repmat(min(reshape(frd_prob, 12, 24)), 12, 1), 288, 1), xtime, reshape(repmat(max(reshape(fru_prob, 12, 24)), 12, 1), 288, 1));
    % h1 =plot(xtime, reshape(repmat(min(reshape(frd_prob, 12, 24)), 12, 1), 288, 1), xtime, reshape(repmat(min(reshape(fru_prob, 12, 24)), 12, 1), 288, 1))
    h1 = plot(xtime, reshape(repmat(mean(reshape(frd_prob, 12, 24), 1), 12, 1), 288, 1), xtime, reshape(repmat(mean(reshape(fru_prob, 12, 24), 1), 12, 1), 288, 1));
    set(h1, 'linewidth', 2, 'color', 'k');

    h2 = plot(xtime, reshape(repmat(fru_determ, 1, 12)', 288, 1), xtime, reshape(repmat(frd_determ, 1, 12)', 288, 1));
    set(h2, 'linewidth', 2, 'color', 'r');

%     fru_oasis = oasis.UP_RTD((oasis.OPR_DT.Year==2019) & (oasis.OPR_DT.Month==select_month) & (oasis.OPR_DT.Day==select_day));
%     frd_oasis = oasis.DOWN_RTD((oasis.OPR_DT.Year==2019) & (oasis.OPR_DT.Month==select_month) & (oasis.OPR_DT.Day==select_day));
%     h3 = plot(xtime, fru_oasis, xtime, -frd_oasis);
%     set(h3, 'linewidth', 2,  'color', 'b');

    legend([h1(1); h2(1)], 'Prob', 'Baseline');
%     legend([h1(1); h2(1); h3(1)], 'FRD prob', 'FRD determ', 'FRP OASIS');
%     legend([h2(1); h3(1)], 'Baseline', 'OASIS');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    
    figure();
    hold on;
    h1 = plot(xtime, reshape(repmat(max(reshape(frd_prob, 12, 24)), 12, 1), 288, 1)./reshape(repmat(frd_determ, 1, 12)', 288, 1));
    h2 = plot(xtime, reshape(repmat(min(reshape(fru_prob, 12, 24)), 12, 1), 288, 1)./reshape(repmat(fru_determ, 1, 12)', 288, 1));
    h3 = plot(xtime, 0.9.*ones(288, 1), 'r');
    set([h1,h2], 'linewidth', 2);
    set(gca, 'YScale', 'log')
    legend([h1, h2], 'FRU', 'FRD');
    ylabel('Ratio of Prob/Baseline');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);

end

end

function test_logn_distribution()
% Fit Mucun's quantiles to log-normal distribution and plot QQplot
% M = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghui\binghui_10min.csv', 1, 0);
% M = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghui\binghui_15min.csv', 1, 0);
M = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghui\binghui_40min.csv', 1, 0);
Mquantiles = M(:, 8:end);
nlactual = M(:, 6);
nldeterm = M(:, 7);

irows = find(M(:, 3) == 1); % Let's just take a look at May 1.
p = 1:99; % 1 to 99 quantiles

figure();
for i = irows(:)'
    Mrow = Mquantiles(i, :);
    logM = log(Mrow);
    y = logM(:);
    x = sqrt(2).*erfinv(2.*p./100-1);
    x = sqrt(2).*erfinv(2.*p(:)./100-1);
    X = [x ones(size(x))];
    b = (X'*X)\X'*y;
    SIG = b(1);
    MU = b(2);
    plot(p./100, logncdf(Mrow, b(2), b(1)), 'k.');
    hold on;
end
plot([0, 1], [0, 1], 'r', 'linewidth', 2);
title('QQplot');
set(findall(gcf,'-property','FontSize'),'FontSize',14);
end

function test_discretize_lognormal()
% Discretize the fitted log-normal distribution from Mucun's data
M = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghui\binghui_10min.csv', 1, 0);
Mquantiles = M(:, 8:end);
nlactual = M(:, 6);
nldeterm = M(:, 7);

binsize = 100; % 100 MW

irows = find(M(:, 3) == 1); % Let's just take a look at May 1.
p = 1:99; % 1 to 99 quantiles

for i = 1:10
    Mrow = Mquantiles(i, :);
    [h_bin, binedge, MU, SIG] = discretize_lognormal(Mrow, p, binsize);
    
    figure();
    plot(binedge, lognpdf(binedge, MU, SIG));
    hold on;
    stairs(binedge, [h_bin, 0]);
    xlabel('Net load (MW)');
    ylabel('Probability density');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    legend('Log-normal', 'Approximation');
end
end

function [h_bin, binedge, MU, SIG] = discretize_lognormal(Mrow, p, binsize)
% Given percentile Mrow and corresponding CDF p, return disretized
% distribution and the fitted log-normal distribution parameters
logM = log(Mrow);
y = logM(:);
x = sqrt(2).*erfinv(2.*p./100-1);
x = sqrt(2).*erfinv(2.*p(:)./100-1);
X = [x ones(size(x))];
b = (X'*X)\X'*y;
SIG = b(1);
MU = b(2);
binedge = [floor(Mrow(1)/binsize) : ceil(Mrow(end)/binsize)].*binsize;
edge_cdf = logncdf(binedge, MU, SIG);
edge_cdf(1) = 0;
edge_cdf(end) = 1; % Because my convolution cannot handle distributions over (0, inf) so we just assume the distribution is bounded here

edge_cdf(end) = 1; % Because my convolution cannot handle distributions over (0, inf) so we just assume the distribution is bounded here
h_bin = (edge_cdf(2:end) - edge_cdf(1:length(edge_cdf)-1)).*1/binsize;

end
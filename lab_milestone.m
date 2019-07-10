function lab_milestone()
% test_logn_distribution();
% test_discretize_lognormal();
% ramping_requirement_rtd_nl()
ramping_requirement_rtpd_nl()
% cell_history = return_30_days(2019, 5, 31);
end

function cell_history= return_30_days(YYYY, MM, DD)
thisday = datetime(YYYY, MM, DD);
if mod(weekday(thisday), 6) == 1
    todayisweekday = false;
else
    todayisweekday = true;
end
cell_history = cell(30, 1);

validnumber = 0;
historyday  = thisday;
while validnumber < 30
    historyday = historyday - day(1);
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

function ramping_requirement_rtpd_nl()
load frp_20190501;
selected_day = 1;
M_b = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghui\binghui_10min.csv', 1, 0); % Binding forecast
M_a = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghui\binghui_40min.csv', 1, 0); % Advisory forecast

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

for selected_day = 1: 30
    irows = find(M_b(:, 3) == selected_day); % Let's just take a look at May 1.
    frd_rtd = nan(size(irows, 1), 1);
    fru_rtd = nan(size(irows, 1), 1);
    for i = irows'
        Mrow_b = Mquantiles_b(i, :);
        Mrow_a = Mquantiles_a(i, :);
        [h_bin_b, binedge_b] = discretize_lognormal(Mrow_b, p, binsize);
        [h_bin_a, binedge_a] = discretize_lognormal(Mrow_a, p, binsize);
        [a_convoluted, t_convoluted] = conv_poly([h_bin_a(:); 0], binedge_a(:), [flipud(h_bin_b(:)); 0], flipud(-binedge_b(:)), binsize); % Distribution of binding - advisory forecast NL

        frd_rtd(i - (selected_day-1)*288) = interp1(cdf_poly(a_convoluted, t_convoluted), t_convoluted, 0.05);  % 5 percentile
        fru_rtd(i - (selected_day-1)*288) = interp1(cdf_poly(a_convoluted, t_convoluted), t_convoluted, 0.95);  % 95 percentile

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
    xtime = datetime(2019,5,selected_day,0,5,0):minutes(5): datetime(2019,5,selected_day+1,0,0,0);
    h1= plot(xtime, reshape(repmat(min(reshape(frd_rtd, 12, 24)), 12, 1), 288, 1), xtime, reshape(repmat(max(reshape(fru_rtd, 12, 24)), 12, 1), 288, 1));
    hold on;
    % plot(1:288, reshape(repmat(min(reshape(frd_rtd, 12, 24)), 12, 1), 288, 1), 1:288, reshape(repmat(min(reshape(fru_rtd, 12, 24)), 12, 1), 288, 1))
    % plot(1:288, reshape(repmat(mean(reshape(frd_rtd, 12, 24), 1), 12, 1), 288, 1), 1:288, reshape(repmat(mean(reshape(fru_rtd, 12, 24), 1), 12, 1), 288, 1))
    h2 = plot(xtime, fru_rtpd_20190501, xtime,-frd_rtpd_20190501);
    set(h1, 'linewidth', 2, 'color', 'k');
    set(h2, 'linewidth', 2, 'color', 'r');
    legend([h1(1); h2(1)], 'FRD prob', 'FRD OASIS');

    % h3 = plot(xtime, frd_rtd, xtime, fru_rtd);
    % set(h3, 'color', 'k');
    % legend([h1(1); h2(1); h3(1)], 'FRD prob', 'FRD OASIS', 'FRP 5-min');
end

end

function ramping_requirement_rtd_nl()
M_b = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghui\binghui_10min.csv', 1, 0); % Binding forecast
M_a = csvread('C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghui\binghui_15min.csv', 1, 0); % Advisory forecast
Mquantiles_b = M_b(:, 8:end);
Mquantiles_a = M_a(:, 8:end);
nlactual_b = M_b(:, 6);
nlactual_a = M_a(:, 6);
nldeterm_b = M_b(:, 7);
nldeterm_a = M_a(:, 7);

p = 1:99; % 1 to 99 quantiles
binsize = 100;


for selected_day = 1: 30
    irows = find(M_b(:, 3) == selected_day); % Let's just take a look at May 1.
    frd_rtd = nan(size(irows, 1), 1);
    fru_rtd = nan(size(irows, 1), 1);

    for i = irows'
        Mrow_b = Mquantiles_b(i, :);
        Mrow_a = Mquantiles_a(i, :);
        [h_bin_b, binedge_b] = discretize_lognormal(Mrow_b, p, binsize);
        [h_bin_a, binedge_a] = discretize_lognormal(Mrow_a, p, binsize);
        [a_convoluted, t_convoluted] = conv_poly([h_bin_a(:); 0], binedge_a(:), [flipud(h_bin_b(:)); 0], flipud(-binedge_b(:)), binsize); % Distribution of binding - advisory forecast NL

        frd_rtd(i - (selected_day-1)*288) = interp1(cdf_poly(a_convoluted, t_convoluted), t_convoluted, 0.05);  % 5 percentile
        fru_rtd(i - (selected_day-1)*288) = interp1(cdf_poly(a_convoluted, t_convoluted), t_convoluted, 0.95);  % 95 percentile

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
    xtime = datetime(2019,5,selected_day,0,5,0):minutes(5): datetime(2019,5,selected_day+1,0,0,0);
    h1= plot(xtime, reshape(repmat(min(reshape(frd_rtd, 12, 24)), 12, 1), 288, 1), xtime, reshape(repmat(max(reshape(fru_rtd, 12, 24)), 12, 1), 288, 1));
    hold on;
    % plot(1:288, reshape(repmat(min(reshape(frd_rtd, 12, 24)), 12, 1), 288, 1), 1:288, reshape(repmat(min(reshape(fru_rtd, 12, 24)), 12, 1), 288, 1))
    % plot(1:288, reshape(repmat(mean(reshape(frd_rtd, 12, 24), 1), 12, 1), 288, 1), 1:288, reshape(repmat(mean(reshape(fru_rtd, 12, 24), 1), 12, 1), 288, 1))
    load frp_20190501;
    h2 = plot(xtime, fru_rtd_20190501, xtime,-frd_rtd_20190501);
    set(h1, 'linewidth', 2, 'color', 'k');
    set(h2, 'linewidth', 2, 'color', 'r');
    legend([h1(1); h2(1)], 'FRD prob', 'FRD OASIS');

    % h3 = plot(xtime, frd_rtd, xtime, fru_rtd);
    % set(h3, 'color', 'k');
    % legend([h1(1); h2(1); h3(1)], 'FRD prob', 'FRD OASIS', 'FRP 5-min');
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
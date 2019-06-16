function lab_IBM()
ibm_april(5);
% CAISO_10_sites();
end

function ibm_april(m)
if m == 4
    dirwork = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_April\power_frcst';
elseif m == 5
    dirwork = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_May\power_frcst';
end
dirhome = pwd;
% siteforconv = {'gen55', 'gen56', 'gen59', 'gen60', 'gen61'}; % This 5 sites correspond to 5 unique IBM sites
siteforconv = {'gen55', 'gen56', 'gen59', 'gen60', 'gen61', 'gen58', 'gen62', 'gen64'}; % These 8 gens have non-zero capacities
ibm_sites   = {'MNCC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz'};
capacity_gen =  [232.44, 107.58, 116.05, 140.4, 151.32, ];
capacity_site = [372.84, 134.88, 116.05, 338.52, 151.32];
cap_scaler = capacity_site./capacity_gen;

% Read data
capacity = nan(length(siteforconv), 1);
cell_data = cell(length(siteforconv), 1);
cd(dirwork);
for i = 1: length(siteforconv)
    csvname = strcat('power_', siteforconv{i}, '.csv');
    M = csvread(csvname, 1, 0);
%     M(:, 7:9) = M(:, 7:9).*cap_scaler(i);
    M(M(:, 2) ~=m, :) = []; % Remove march data
    tarray = datetime(M(:, 1), M(:, 2), M(:, 3), M(:, 4), M(:, 5), M(:, 6), 'TimeZone', 'America/Los_Angeles');
%     tarray.TimeZone = 'America/Los_Angeles';
%     M(:, 1) = tarray.Year;
%     M(:, 2) = tarray.Month;
%     M(:, 3) = tarray.Day;
%     M(:, 4) = tarray.Hour;
%     M(:, 5) = tarray.Minute;
%     M(:, 6) = tarray.Second;
    cell_data{i} = M;
    capacity(i) = max(M(:, end));
end
cd(dirhome);

nT = size(M, 1);

% Identify valid rows
i_valid   = true(nT, 1);
i_allzero = true(nT, 1);
for i = 1: length(cell_data)
    M = cell_data{i};
    i_valid   = i_valid & (M(:, 7) <= M(:, 8)) & (M(:, 8) <= M(:, 9)); 
    i_allzero = i_allzero & (M(:, 7) <=1E-3) & (M(:, 8) <=1E-3) & (M(:, 9) <=1E-3);
end
t_idx = 1: 1: nT;
% t_idx = t_idx(i_valid);

% Direct sum, if perfectly correlated
M_agg = zeros(size(M(:, 7:9)));
for i = 1: length(cell_data)
    M = cell_data{i};
    i_valid   = i_valid & (M(:, 7) <= M(:, 8)) & (M(:, 8) <= M(:, 9)); 
    i_allzero = i_allzero & (M(:, 7) <=1E-3) & (M(:, 8) <=1E-3) & (M(:, 9) <=1E-3);
    M_agg = M_agg + M(:, 7:9);
end

cell_a = cell(nT, 1);
cell_t = cell(nT, 1);
x_agg_950 = nan(nT, 1);
x_agg_050 = nan(nT, 1);
% Do convolution
% for t = 1: nT
for d = 1: 29 % Day 1 to 29
    x_agg_950_d = nan(24*4, 1);
    x_agg_050_d = nan(24*4, 1);
    for t = t_idx(i_valid&(M(:, 3)==d))
    % for t = 2466-28 % This is the row with negative distribution, 28 is the number of rows in March
        n = 0;
        a_convoluted = [];
        t_convoluted = [];
%         fig = figure();
        for i = 1:length(cell_data)
            M = cell_data{i};
            cap = capacity(i);
            normalized_binwidth = 0.04;

            % Use logit-normal function to fit the distribution
            x500_normalized = M(t, 8)/cap;
            x050_normalized = M(t, 7)/cap;
            x950_normalized = M(t, 9)/cap;
            [normalized_h_bin, normalized_binedge, flag_success] = discretize_dist_logitnormal(x050_normalized, x500_normalized, x950_normalized, normalized_binwidth);
            if flag_success
                h_bin = normalized_h_bin./cap;
                binedge = normalized_binedge.*cap;
                a_fitted = [h_bin(:); 0];
                t_fitted = binedge;

                if (i == 1) || isempty(t_convoluted)
                    a_convoluted = a_fitted;
                    t_convoluted = t_fitted;
                else
                    [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, a_fitted, t_fitted, normalized_binwidth*cap);
                end
                [a_convoluted, t_convoluted] = reformulate_poly(a_convoluted, t_convoluted);
                n = n + 1;
            else
                continue
            end
%             if (x500~=0) && (x050~=0)
%                 param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.050-0.5)) 1]\[log(x500/(1-x500)); log(x050/(1-x050))];
%             elseif (x050~=0) && (x950~=0)
%                 if x950==cap
%                     x950=x950-eps;
%                 end
%                 param = [sqrt(2)*erfinv(2*(0.050-0.5)) 1; sqrt(2)*erfinv(2*(0.950-0.5)) 1]\[log(x050/(1-x050)); log(x950/(1-x950))];
%             elseif (x500~=0) && (x950~=0)
%                 if x950==cap
%                     x950=x950-eps;
%                 end
%                 param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.950-0.5)) 1]\[log(x500/(1-x500)); log(x950/(1-x950))];
%             else
%                 continue % We need at least two non-zero points, otherwise we consider it as zero
%             end
%             n = n + 1;
%             SIGMA = param(1);
%             MU = param(2);
%         %     plot_logit(MU, SIGMA);
%         %     hold on;
% 
%             % Prepare for convolution        
%             normalized_binwidth = 0.04;
%             normalized_binedge = 0:normalized_binwidth:1; % This is normalized bin edges
%             edge_cdf = cdf_logit_normal(normalized_binedge, MU, SIGMA);
%             normalized_h_bin = (edge_cdf(2:end) - edge_cdf(1:length(edge_cdf)-1)).*1/normalized_binwidth;
% 
%             h_bin = normalized_h_bin./cap;
%             binedge = normalized_binedge.*cap;

%             if (i == 1) || isempty(t_convoluted)
%                 a_convoluted = [h_bin(:); 0];
%                 t_convoluted = binedge;
%             else
%                 [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [h_bin(:); 0], binedge, normalized_binwidth*cap);
%             end
    %         a_convoluted(end, :) = 0;
%             ax = subplot(3, 3, i);
%             plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
%             [a_convoluted, t_convoluted] = reformulate_poly(a_convoluted, t_convoluted);
%             title(n);
        end
        cell_a{t} = a_convoluted;
        cell_t{t} = t_convoluted;
        
        t_thatday = t - (d-1)*24*4;
        if ~isempty(a_convoluted)
            cdf_agg = cdf_poly(a_convoluted, t_convoluted);
            [~, imin] = min(abs(cdf_agg - 0.05));
            x_agg_050(t) = t_convoluted(imin);
            x_agg_050_d(t_thatday) = t_convoluted(imin);
            [~, imin] = min(abs(cdf_agg - 0.95));
            x_agg_950(t) = t_convoluted(imin);
            x_agg_950_d(t_thatday) = t_convoluted(imin);
        end
        if i_allzero(t)
            x_agg_050(t) = 0;
            x_agg_950(t) = 0;
            x_agg_050_d(t_thatday) = 0;
            x_agg_950_d(t_thatday) = 0;
        end
        fprintf('%4g/%g Done.\n', t, nT);
    end
    
    % One day's snapshot
%     figure();
%     tarray_d = tarray((M(:, 2)==4)&(M(:, 3)==d)); % T array for the day
%     h1 = plot(tarray_d, x_agg_050_d, 'b');
%     hold on;
%     h2 = plot(tarray_d, x_agg_950_d, 'r');
%     hold on;
%     h3 = plot(tarray_d, M_agg(M(:, 3)==d, :));
%     set(h3, {'linestyle'}, {'--'; '--'; '--'});
%     set(h3, {'color'}, {'b'; 'k'; 'r'});
%     legend([h1; h2; h3], {'CI 5%'; 'CI 95%'; 'SUM 5%'; 'SUM MEAN'; 'SUM 95%'});
%     ylabel('kW');
%     title(d);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix: Assume all nan is zero, i.e., when low 95 CI and mean are 0, and up
% 95 CI is not zero, we assume it is zero
x_agg_050(isnan(x_agg_050)) = 0;
x_agg_950(isnan(x_agg_950)) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d = 1: 29
    figure();
    x_agg_950_d = x_agg_950((M(:, 2)==m)&(M(:, 3)==d));
    x_agg_050_d = x_agg_050((M(:, 2)==m)&(M(:, 3)==d));

    tarray_d = tarray((M(:, 2)==m)&(M(:, 3)==d)); % T array for the day
    h1 = plot(tarray_d, x_agg_050_d, 'b');
    hold on;
    h2 = plot(tarray_d, x_agg_950_d, 'r');
    hold on;
    h3 = plot(tarray_d, M_agg(M(:, 3)==d, :));
    set(h3, {'linestyle'}, {'--'; '--'; '--'});
    set(h3, {'color'}, {'b'; 'k'; 'r'});
    legend([h1; h2; h3], {'CI 5%'; 'CI 95%'; 'SUM 5%'; 'SUM MEAN'; 'SUM 95%'});
    ylabel('kW');
    title(d);

end

figure();
tarray = datetime(M(:, 1), M(:, 2), M(:,3), M(:, 4), M(:, 5), 0);
h1 = plot(tarray, x_agg_050, 'b');
hold on;
h2 = plot(tarray, x_agg_950, 'r');
legend([h1, h2], {'CI 5%', 'CI 95%'});
ylabel('kW');
% ylim([-10, 800]);
% sum(isnan(x_agg_050))/size(x_agg_050, 1);
end

function CAISO_10_sites()
%% Test IBM's data
dirwork = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_old\sample';
dirhome = pwd;
listcsv = {...
    'sample_35.38_-120.18.csv', ...
    'sample_38.47_-122.71.csv', ...
    'sample_40.25_-123.31.csv', ...
    'sample_34.45_-119.7.csv', ...
    'sample_40.71_-123.92.csv', ...
    'sample_37.41_-119.74.csv', ...
    'sample_34.31_-117.5.csv', ...
    'sample_34.12_-117.94.csv', ...
    'sample_35.53_-118.63.csv', ...
    'sample_39.12_-123.07.csv', ...
};

cap = 250; % kW

%% Show fitted distributions
row = 36;
for i = 1:length(listcsv)
    fig = figure();
    csvname = listcsv{i};
    cd(dirwork);
    M = csvread(csvname, 1, 1);
    cd(dirhome);

    % Use logit-normal function to fit the distribution
    hour = mod(M(row, 4) - 8+24, 24);
    x500 = M(row, 6)/cap;
    x025 = M(row, 7)/cap;
    x975 = M(row, 8)/cap;
    if (x500~=0) && (x025~=0)
        param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.025-0.5)) 1]\[log(x500/(1-x500)); log(x025/(1-x025))];
    elseif (x025~=0) && (x975~=0)
        if x975==cap
            x975=x975-eps;
        end
        param = [sqrt(2)*erfinv(2*(0.025-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[log(x025/(1-x025)); log(x975/(1-x975))];
    elseif (x500~=0) && (x975~=0)
        if x975==cap
            x975=x975-eps;
        end
        if x975==cap
            x975=x975-eps;
        end
        param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[log(x500/(1-x500)); log(x975/(1-x975))];
    else
        continue % We need at least two non-zero points.
    end
    SIGMA = param(1);
    MU = param(2);
%     plot_logit(MU, SIGMA);
%     hold on;

    % Prepare for convolution
    binwidth = 10; % kW
    normalized_binedge = 0:binwidth/cap:1; % This is normalized bin edges
    edge_cdf = cdf_logit_normal(normalized_binedge, MU, SIGMA);
    normalized_h_bin = (edge_cdf(2:end) - edge_cdf(1:length(edge_cdf)-1)).*cap/binwidth;
    
    h_bin = normalized_h_bin./250;
    binedge = normalized_binedge.*cap;
    
    % Plot discretized and continuous figures
    stairs_conv_poly(fig, [h_bin(:); 0], binedge, 'k');
    hold on;
    plot(0:binwidth:cap, pdf_logit_normal(normalized_binedge, MU, SIGMA)./cap, 'b', 'LineWidth', 1.5);
    y1=get(gca,'ylim');
    plot([x500 x500].*cap,y1 ,'r', 'LineWidth', 1.5)
%     plot([x025 x025].*cap,y1)
%     plot([x975 x975].*cap,y1)
    patch([x025, x975, x975, x025].*cap, [y1(1), y1(1), y1(2), y1(2)], 'k', 'FaceAlpha', 0.2);
    title(i);
    xlabel('kW');
    set(gca, 'FontSize', 12');
    box on;
end

%% See the distribution of aggretated solar power output
rows_selected = 30:76;
% rows_selected = 36;

for row = rows_selected
    fig = figure;
    n = 0; % Number of plants with non-negative output
    for i = 1:length(listcsv)
        csvname = listcsv{i};
        cd(dirwork);
        M = csvread(csvname, 1, 1);
        cd(dirhome);

        % Use logit-normal function to fit the distribution
        hour = mod(M(row, 4) - 8+24, 24);
        x500 = M(row, 6)/cap;
        x025 = M(row, 7)/cap;
        x975 = M(row, 8)/cap;
        if (x500~=0) && (x025~=0)
%             param = [sqrt(2)*erfinv(2*(0.050-0.5)) 1; sqrt(2)*erfinv(2*(0.025-0.5)) 1]\[x500/(1-x500); x025/(1-x025)];
            param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.025-0.5)) 1]\[log(x500/(1-x500)); log(x025/(1-x025))];
        elseif (x025~=0) && (x975~=0)
            if x975==cap
                x975=x975-eps;
            end
            param = [sqrt(2)*erfinv(2*(0.025-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[log(x025/(1-x025)); log(x975/(1-x975))];
%             param = [sqrt(2)*erfinv(2*(0.025-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[x025/(1-x025); x975/(1-x975)];
        elseif (x500~=0) && (x975~=0)
            if x975==cap
                x975=x975-eps;
            end
            if x975==cap
                x975=x975-eps;
            end
            param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[log(x500/(1-x500)); log(x975/(1-x975))];
%             param = [sqrt(2)*erfinv(2*(0.050-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[x500/(1-x500); x975/(1-x975)];
        else
            continue % We need at least two non-zero points.
        end
        n = n + 1;
        SIGMA = param(1);
        MU = param(2);
    %     plot_logit(MU, SIGMA);
    %     hold on;

        % Prepare for convolution
%         binwidth = 10; % kW
%         binedge = 0:binwidth/cap:1; % This is normalized bin edges
%         edge_cdf = cdf_logit_normal(binedge, MU, SIGMA);
%         binedge = binedge.*cap;
%         bin_pdf = (edge_cdf(2:end) - edge_cdf(1:length(edge_cdf)-1))./binwidth;
        
        binwidth = 10; % kW
        normalized_binedge = 0:binwidth/cap:1; % This is normalized bin edges
        edge_cdf = cdf_logit_normal(normalized_binedge, MU, SIGMA);
        normalized_h_bin = (edge_cdf(2:end) - edge_cdf(1:length(edge_cdf)-1)).*cap/binwidth;
        
        h_bin = normalized_h_bin./250;
        binedge = normalized_binedge.*cap;



        if i == 1
            a_convoluted = [h_bin(:); 0];
            t_convoluted = binedge;
        else
            [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [h_bin(:); 0], binedge, binwidth);
        end
        ax = subplot(3, 4, i);
        plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
        title(n);
    end
end

end

function stairs_conv_poly(fig, a, t, c)
figure(fig); hold on;
x_l_all = nan(size(a, 1)-1, 1);
y_l_all = nan(size(a, 1)-1, 1);
x_r_all = nan(size(a, 1)-1, 1);
y_r_all = nan(size(a, 1)-1, 1);
for i = 1:(size(a, 1)-1)
    x_l = t(i);
    x_r = t(i+1);
    x = linspace(x_l, x_r, 10);
    y = polyval(fliplr(a(i, :)), x);
%     h = stairs(x, y, 'color', c);
    x_l_all(i) = x_l;
    x_r_all(i) = x_r;
    y_l_all(i) = y(1);
    y_r_all(i) = y(end);
end
stairs(x_l_all, y_l_all, 'k')

end

function [normalized_h_bin, normalized_binedge, flag_success] = discretize_dist_logitnormal(x050, x500, x950, normalized_binwidth)
% Fit IBM's forecasted percentiles to distributions
% All inputs should be normalized, i.e., between 0 and 1
flag_lack_data = false;
if (x050~=0) && (x950~=0)
    if x950==1
        x950=x950-eps;
    end
    param = [sqrt(2)*erfinv(2*(0.050-0.5)) 1; sqrt(2)*erfinv(2*(0.950-0.5)) 1]\[log(x050/(1-x050)); log(x950/(1-x950))];
elseif (x500~=0) && (x050~=0)
    param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.050-0.5)) 1]\[log(x500/(1-x500)); log(x050/(1-x050))];
elseif (x500~=0) && (x950~=0)
    if x950==1
        x950=x950-eps;
    end
    param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.950-0.5)) 1]\[log(x500/(1-x500)); log(x950/(1-x950))];
else
    flag_lack_data = true; % We need at least two non-zero points, otherwise we consider it as zero
end

if ~flag_lack_data
    SIGMA = param(1);
    MU = param(2);
    %     plot_logit(MU, SIGMA);
    %     hold on;

    % Prepare for convolution        
    normalized_binedge = 0:normalized_binwidth:1; % This is normalized bin edges
    edge_cdf = cdf_logit_normal(normalized_binedge, MU, SIGMA);
    normalized_h_bin = (edge_cdf(2:end) - edge_cdf(1:length(edge_cdf)-1)).*1/normalized_binwidth;

    flag_success = true;
else
    normalized_binedge = nan;
    normalized_h_bin = nan;
    flag_success = false;
end

end

function y = cdf_logit_normal(x, MU, SIGMA)

y = 0.5.*( 1+erf((logit(x)-MU)./(sqrt(2)*SIGMA)) );
end

function y = pdf_logit_normal(x, MU, SIGMA)
logitx = logit(x);
y = 1/(SIGMA*sqrt(2*pi))*1./(x.*(1-x)).*exp(-(logitx-MU).^2./(2*SIGMA^2));
end

function plot_logit(MU, SIGMA)
x = 0.01:0.01:0.99;
% logitx = logit(x);
% y = 1/(SIGMA*sqrt(2*pi))*1./(x.*(1-x)).*exp(-(logitx-MU).^2./(2*SIGMA^2));
y = pdf_logit_normal(x, MU, SIGMA);
plot(x, y);
end

function y = logit(x)
y = log(x./(1-x));
end
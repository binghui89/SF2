function lab_IBM()
CAISO_10_sites();
end

function CAISO_10_sites()
%% Test IBM's data
dirwork = 'C:\Users\bxl180002\Downloads\RampSolar\IBM';
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
            [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [h_bin(:); 0], binedge);
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
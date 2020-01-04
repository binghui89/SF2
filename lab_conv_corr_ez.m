function lab_conv_corr_ez()

conv_corr_ez_nd()

end

function conv_corr_ez_nd()
addpath('ndhist\histcn'); % Add the path to the N-Dimensional histcount function
load caiso;
all_years = [2016; 2017; 2018; 2019];
bin_width = 200; % MW

% Let's use GW for better visual effects
bin_width = bin_width./1E3;
for i = 1: 4
    btm{i} = btm{i}./1000;
    loads{i} = loads{i}./1000;
    wind{i} = wind{i}./1000;
    if i <= 3
        pv{i} = pv{i}./1E3;
    end
    if i <= 2
        st{i} = st{i}./1E3;
    end
    stpv2019 = stpv2019./1E3;
end

year_index = 2; % Change to 4 if want to run all four years
switch year_index
    case 1 % 2016
        ts = [loads{year_index}, -wind{year_index}, -st{year_index}, -btm{year_index}, -pv{year_index}];
    case 2 % 2017
        ts = [loads{year_index}, -wind{year_index}, -st{year_index}, -btm{year_index}, -pv{year_index}];
    case 3 % 2018
        ts = [loads{year_index}, -wind{year_index}, -btm{year_index}, -pv{year_index}];
    case 4 % 2019
        ts = [loads{year_index}, -wind{year_index}, -btm{year_index}, -stpv2019];
end

N = size(ts, 2); % Size of components
ts_cdf = nan(size(ts, 1), N); % Time series CDF
gr_pdf = cell(N, 1); % Probability density over each grid
gr_p   = cell(N, 1); % Probability of each grid;
gr_bincenter = cell(N, 1); % Bin center of each grid
gr_edges = cell(N, 1); % Bin edges of each grid, length + 1
gr_cdf = cell(N, 1); % Cumulative probability at the edges of grids

for i = 1:N
    [tmp_pdf, tmp_bincenter, tmp_edges] = return_pdf(ts(:, i), bin_width);
    gr_pdf{i} = tmp_pdf(:);
    gr_bincenter{i} = tmp_bincenter(:);
    gr_edges{i} = tmp_edges(:);
    gr_p{i} = gr_pdf{i}.*bin_width;
    gr_cdf{i} = [0; cumsum(gr_p{i})];
    gr_cdf{i}(end) = 1;
    ts_cdf(:, i) = interp1(gr_edges{i}, gr_cdf{i}, ts(:, i));
end

% Obtain the mean copula pdf (empirical copula)
ts_cdf(ts_cdf==1) = 1-eps; % 1 is the upper limit of the last bin, and if any data point falls on that boundary, histcn will add one more bin.

switch N
    case 4
        count = histcn(ts_cdf, gr_cdf{1}, gr_cdf{2}, gr_cdf{3}, gr_cdf{4});
        [X1, X2, X3, X4] = ndgrid(gr_pdf{1}.*bin_width, gr_pdf{2}.*bin_width, gr_pdf{3}.*bin_width, gr_pdf{4}.*bin_width); % Grid size along each dimension
        grid_vol = X1.*X2.*X3.*X4; % Grid volume
        clear X1 X2 X3 X4; % Save some memory
    case 5
        count = histcn(ts_cdf, gr_cdf{1}, gr_cdf{2}, gr_cdf{3}, gr_cdf{4}, gr_cdf{5});
        [X1, X2, X3, X4, X5] = ndgrid(gr_pdf{1}.*bin_width, gr_pdf{2}.*bin_width, gr_pdf{3}.*bin_width, gr_pdf{4}.*bin_width, gr_pdf{5}.*bin_width); % Grid size along each dimension
        grid_vol = X1.*X2.*X3.*X4.*X5; % Grid volume
        clear X1 X2 X3 X4 X5; % Save some memory
end

copula_pdf_nd = count./(sum(count(:)).*grid_vol); % This is the n-dim empirical copula
copula_pdf_nd(grid_vol==0) = 0; % Sometimes the grid volumn is 0 (when the pdf in a grid is 0).
clear count grid_vol;

[x_conv, p_conv] = conv_ez_corr(gr_bincenter, gr_p, bin_width, copula_pdf_nd);

e_conv = [x_conv(1) - bin_width/2; x_conv + bin_width/2]; % Bin edges of convolved results
count_actual = histcounts(sum(ts, 2), e_conv); % Actual counts

figure();
bar(x_conv ,count_actual./sum(count_actual));
hold on;
plot(x_conv, p_conv, 'k');
end


function lab_conv_corr_ez()

conv_corr_ez_nd()

% copula_type = 'empirical';
% conv_corr_ez_2d(copula_type);
end

function conv_corr_ez_nd()
addpath('ndhist\histcn'); % Add the path to the N-Dimensional histcount function
load caiso;
all_years = [2016; 2017; 2018; 2019];
bin_width = 100; % MW
year_index = 1;

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

year_index = 1; % Change to 4 if want to run all four years
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

ts_cdf(ts_cdf==1) = 1-eps; % 1 is the upper limit of the last bin, and if any data point falls on that boundary, histcn will add one more bin.
count = histcn(ts_cdf, gr_cdf{1}, gr_cdf{2}, gr_cdf{3}, gr_cdf{4}, gr_cdf{5});
[X1, X2, X3, X4, X5] = ndgrid(gr_pdf{1}.*bin_width, gr_pdf{2}.*bin_width, gr_pdf{3}.*bin_width, gr_pdf{4}.*bin_width, gr_pdf{5}.*bin_width); % Grid size along each dimension
grid_vol = X1.*X2.*X3.*X4.*X5; % Grid volume
clear X1 X2 X3 X4 X5; % Save some memory

copula_pdf_nd = count./(sum(count(:)).*grid_vol); % This is the n-dim empirical copula
clear count grid_vol;

ind = zeros(numel(copula_pdf_nd), N);

[ind(:, 1), ind(:, 2), ind(:, 3), ind(:, 4), ind(:, 5)] = ind2sub(size(copula_pdf_nd), [1: numel(copula_pdf_nd)]');
sumind = sum(ind, 2);

flat_p_alldim = zeros(numel(copula_pdf_nd), N);
for i = 1: N
    flat_p_alldim(:, i) = gr_p{i}(ind(:, i));
end
flat_copulapdf_p = [copula_pdf_nd(:), flat_p_alldim];
clear ind copula_pdf_nd flat_p_alldim;

ar_delta = [N: max(sumind)]';
x_conv0 = 0;
for i = 1:N
    x_conv0 = x_conv0 + gr_bincenter{i}(1) - bin_width;
end
x_conv = x_conv0 + ar_delta.*bin_width; % Bin center of convolved results
p_conv = zeros(numel(ar_delta), 1); % Probability of each bin in the results
for k = 1: numel(ar_delta)
    disp(k);
    thissum = ar_delta(k);
    p_conv(k) = sum(prod(flat_copulapdf_p(sumind==thissum, :), 2), 1);
end
e_conv = [x_conv(1) - bin_width/2; x_conv + bin_width/2]; % Bin edges of convolved results
count_actual = histcounts(sum(ts, 2), e_conv); % Actual counts
clear flat_p_alldim;

figure();
bar(x_conv ,count_actual./sum(count_actual));
hold on;
plot(x_conv, p_conv, 'k');
end

function conv_corr_ez_2d(copula_type)

% Net load convolution with correlation, using CAISO's flexibility assessment dataset
load caiso;
all_years = [2016; 2017; 2018; 2019];
bin_width = 50; % MW
year_index = 1;

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

year_index = 1; % Change to 4 if want to run all four years

fig = figure();
%     ts = array; % Place holder
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

time_convcorr = 0;
time_conv = 0; 
for i = 2: size(ts, 2)
    ts_1 = sum(ts(:, 1:i-1), 2);
    ts_2 = ts(:, i);
    ts_3 = ts_2 + ts_1; % Actual net load

    % Calculate pdf
    [gr_pdf_3, gr_bincenter_3, gr_edges_3] = return_pdf(ts_3, bin_width);
    [gr_pdf_2, gr_bincenter_2, gr_edges_2] = return_pdf(ts_2, bin_width);
    [gr_pdf_1, gr_bincenter_1, gr_edges_1] = return_pdf(ts_1, bin_width);

    % Find copula
    gr_cdf_1 = [0 cumsum(gr_pdf_1.*bin_width)];
    gr_cdf_1(end) = 1;

    gr_cdf_2 = [0 cumsum(gr_pdf_2.*bin_width)];
    gr_cdf_2(end) = 1;

    ts_cdf_1 = interp1(gr_edges_1, gr_cdf_1, ts_1);
    ts_cdf_2  = interp1(gr_edges_2,  gr_cdf_2,  ts_2);

    i_valid = (ts_cdf_2<1) & (ts_cdf_2>0) & (ts_cdf_1<1) & (ts_cdf_1>0);

    % Determine convolved functions
    if i == 2
%         a_prev = [gr_pdf_1(:); 0];
%         t_prev = gr_edges_1(:);
%         a_prev_corr = a_prev;
%         t_prev_corr = t_prev;

        gr_cdf_forconv_1 = gr_cdf_1;
        gr_cdf_forconv_2 = gr_cdf_2;
    else
%         a_prev = a_conv;
%         t_prev = t_conv;
%         a_prev_corr = a_convcorr;
%         t_prev_corr = t_convcorr;
% 
%         gr_cdf_forconv_1 = cdf_poly(a_prev_corr, t_prev_corr);
        gr_cdf_forconv_1 = cumsum(p_conv_corr);
        gr_cdf_forconv_2 = gr_cdf_2;
    end


    % Calculate discretized copula pdf
    Ix = size(gr_cdf_forconv_1(:), 1);
    Iy = size(gr_cdf_forconv_2(:), 1);

%         rhohat = copulafit('Gaussian', [ts_cdf_2(i_valid) ts_cdf_1(i_valid)]);
%         [copula_x_grid, copula_y_grid] = meshgrid(gr_cdf_forconv_1, gr_cdf_forconv_2);
%         copula_cdf = copulacdf('Gaussian', [copula_x_grid(:) copula_y_grid(:)], rhohat);
%         copula_cdf = reshape(copula_cdf, Iy, Ix);

    switch copula_type
        case 'Gaussian' % Gaussian copula
            rhohat = copulafit('Gaussian', [ts_cdf_2(i_valid) ts_cdf_1(i_valid)]);
            [copula_x_grid, copula_y_grid] = meshgrid(gr_cdf_forconv_1, gr_cdf_forconv_2);
            copula_cdf = copulacdf('Gaussian', [copula_x_grid(:) copula_y_grid(:)], rhohat);
            copula_cdf = reshape(copula_cdf, Iy, Ix);
        case 'empirical' % Empirical copula
            cdfcounts = histcounts2(ts_cdf_1, ts_cdf_2, gr_cdf_forconv_1, gr_cdf_forconv_2);
            cdfcounts = cdfcounts'; % Matlab uses different x and y directions
        otherwise % Archimedean copulas and t copula, param2 is the 2nd parameter
            [rhohat, param2] = copulafit(copula_type, [ts_cdf_2(i_valid) ts_cdf_1(i_valid)]);
            [copula_x_grid, copula_y_grid] = meshgrid(gr_cdf_forconv_1, gr_cdf_forconv_2);
            copula_cdf = copulacdf(copula_type, [copula_x_grid(:) copula_y_grid(:)], rhohat, param2);
            copula_cdf = reshape(copula_cdf, Iy, Ix);
    end

    mean_copula_pdf = nan(Iy-1, Ix-1);
    switch copula_type
        case 'empirical'
            for ii = 1: Iy-1
                for jj = 1: Ix-1
%                         grid_area = (gr_cdf_2(ii+1) - gr_cdf_2(ii))*(gr_cdf_1(jj+1) - gr_cdf_1(jj)); % Note grid area is not homogeneous!
                    grid_area = (gr_cdf_forconv_2(ii+1) - gr_cdf_forconv_2(ii))*(gr_cdf_forconv_1(jj+1) - gr_cdf_forconv_1(jj)); % Note grid area is not homogeneous!
                    mean_copula_pdf(ii,jj) = cdfcounts(ii,jj)/(grid_area*sum(cdfcounts(:)));
                end
            end
        otherwise
        for ii = 1: Iy-1
            for jj = 1: Ix-1
                grid_area = (gr_cdf_forconv_2(ii+1) - gr_cdf_forconv_2(ii))*(gr_cdf_forconv_1(jj+1) - gr_cdf_forconv_1(jj)); % Note grid area is not homogeneous!
                mean_copula_pdf(ii,jj) = (copula_cdf(ii+1, jj+1) + copula_cdf(ii, jj) - copula_cdf(ii+1, jj) - copula_cdf(ii, jj+1))/grid_area;
            end
        end 
    end

    p_1 = gr_pdf_1.*bin_width;
    p_2 = gr_pdf_2.*bin_width;

    flat_mean_copula_pdf = mean_copula_pdf(:);
    [flat_i, flat_j] = ind2sub(size(mean_copula_pdf), [1: numel(mean_copula_pdf)]');
    copula_and_p = [flat_mean_copula_pdf p_2(flat_i)' p_1(flat_j)'];
    sum_sub = sum([flat_i, flat_j], 2);
    conv_sub = [1+1:1:numel(p_1)+numel(p_2)]';
    p_conv_corr = zeros(length(conv_sub), 1);
    for k = 1: length(p_conv_corr)
        copula_and_p_tmp = copula_and_p(sum_sub==conv_sub(k), :);
        p_conv_corr(k) = sum(prod(copula_and_p_tmp, 2));
    end
    x_conv_corr = (gr_bincenter_1(1) + gr_bincenter_2(1) - 2*bin_width) + bin_width.*conv_sub;
    
    x_conv = gr_bincenter_2(1)+gr_bincenter_1(1) : bin_width: gr_bincenter_2(end)+gr_bincenter_1(end);
    p_conv = conv(p_1, p_2);

    plot(x_conv, p_conv); % Convolution no correlation
    hold on;
    plot(gr_bincenter_3, gr_pdf_3.*bin_width); % Actual distribution
    plot(x_conv_corr, p_conv_corr, 'k')
end
end

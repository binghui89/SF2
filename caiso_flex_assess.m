function caiso_flex_assess()
% This script demonstrates the convolution-based method using the CAISO
% flexibility assessment dataset. ISGT paper.

% validate_caiso_data();
% convolution_2component_no_corr();
% convolution_netload_no_corr();
% convolution_2component_no_corr_new();

% copula_type = {'t', 'Clayton', 'Frank', 'Gumbel'};
% copula_type = {'empirical'};
% copula_type = {'Gaussian'};
% for i = 1: length(copula_type)
%     display(copula_type{i});
%     convolution_2component_corr(copula_type{i});
% end

% copula_type = {'Gaussian'};
copula_type = {'empirical'};
for i = 1: length(copula_type)
    display(copula_type{i});
    convolution_netload_corr(copula_type{i});
end
end


function validate_caiso_data()
% The CAISO data: Data validation
load caiso;
all_years = [2016; 2017; 2018; 2019];

figure();
for i = 1: length(loads)
    subplot(2, 2, i);
    plot(loads{i, 1}./1E3);
    if i >= 3
        xlabel('Time (min)');
    end
    if mod(i, 2) == 1
        ylabel('Power (GW)');
    end
    title(all_years(i));
end
suptitle('Load');

figure();
for i = 1: length(wind)
    subplot(2, 2, i);
    plot(wind{i, 1}./1E3);
    if i >= 3
        xlabel('Time (min)');
    end
    if mod(i, 2) == 1
        ylabel('Power (GW)');
    end
    title(all_years(i));
end
suptitle('Wind');

figure();
for i = 1: length(pv)
    subplot(2, 2, i);
    plot(pv{i, 1}./1E3);
    if i >= 3
        xlabel('Time (min)');
    end
    if mod(i, 2) == 1
        ylabel('Power (GW)');
    end
    title(all_years(i));
end
subplot(2, 2, 4);
plot(stpv2019./1E3);
xlabel('Time (min)');
% ylabel('Power (GW)');
title(all_years(i));
suptitle('Solar PV/Solar');

figure();
for i = 1: length(st)
    subplot(2, 2, i);
    plot(st{i, 1}./1E3);
    if i >= 3
        xlabel('Time (min)');
    end
    if mod(i, 2) == 1
        ylabel('Power (GW)');
    end
    title(all_years(i));
end
suptitle('Solar thermal');

figure();
for i = 1: length(btm)
    subplot(2, 2, i);
    plot(btm{i, 1}./1E3);
    if i >= 3
        xlabel('Time (min)');
    end
    if mod(i, 2) == 1
        ylabel('Power (GW)');
    end
    title(all_years(i));
end
suptitle('Behind-the-meter');

end

function convolution_2component_no_corr()
% The CAISO data, demonstration of convolution: just two components and no
% correlation
load caiso;
all_years = [2016; 2017; 2018; 2019];
bin_width = 0.050; % GW
year_index = 1;

all_figures = nan(length(all_years), 1);
all_titles = {'Load - wind', 'Load - solar PV/solar', 'Load - solar thermal', 'Load - BTM'};
for i = 1: 4 % We have 4 components: PV, thermal, wind, and BTM
    all_figures(i) = figure();
end

for year_index = 1: 4

    [l_pdf, l_bincenter, l_edges] = return_pdf(loads{year_index}./1E3, bin_width);
    [w_pdf, w_bincenter, w_edges] = return_pdf(-wind{year_index}./1E3, bin_width);
    [btm_pdf, btm_bincenter, btm_edges] = return_pdf(-btm{year_index}./1E3, bin_width);
    if year_index == 3
        [pv_pdf, pv_bincenter, pv_edges] = return_pdf(-pv{year_index}./1E3, bin_width);
    elseif year_index == 4
        [pv_pdf, pv_bincenter, pv_edges] = return_pdf(-stpv2019./1E3, bin_width);
    else
        [pv_pdf, pv_bincenter, pv_edges] = return_pdf(-pv{year_index}./1E3, bin_width);
        [st_pdf, st_bincenter, st_edges] = return_pdf(-st{year_index}./1E3, bin_width);
    end

    % Just load and wind
    fig = figure(all_figures(1));
    subplot(2, 2, year_index);
    nl = loads{year_index} - wind{year_index};
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
    [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [w_pdf(:); 0], w_edges(:), bin_width);
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');

    % Just load and solar PV
    fig = figure(all_figures(2));
    subplot(2, 2, year_index);
    if year_index <= 3
        nl = loads{year_index} - pv{year_index};
    else
        nl = loads{year_index} - stpv2019;
    end
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
    [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [pv_pdf(:); 0], pv_edges(:), bin_width);
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');

    % Just load and solar thermal
    fig = figure(all_figures(3));
    if year_index <= 2
        subplot(2, 2, year_index);
        nl = loads{year_index} - st{year_index};
        [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
        [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [st_pdf(:); 0], st_edges(:), bin_width);
        plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
        plot(nl_bincenter, nl_pdf, 'r');
        title(all_years(year_index));
        xlabel('GW');
    end

    % Just load and BTM data
    fig = figure(all_figures(4));
    subplot(2, 2, year_index);
    nl = loads{year_index} - btm{year_index};
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
    [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [btm_pdf(:); 0], btm_edges(:), bin_width);
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');
end

for i = 1: length(all_figures)
    figure(all_figures(i));
    suptitle(all_titles{i});
end
end

function convolution_netload_no_corr()
% The CAISO data, demonstration of convolution: net load
load caiso;
all_years = [2016; 2017; 2018; 2019];
bin_width = 0.050; % GW
year_index = 1;

% Net load
fig = figure();
for year_index = 1: 4
    [l_pdf, l_bincenter, l_edges] = return_pdf(loads{year_index}./1E3, bin_width);
    [w_pdf, w_bincenter, w_edges] = return_pdf(-wind{year_index}./1E3, bin_width);
    [btm_pdf, btm_bincenter, btm_edges] = return_pdf(-btm{year_index}./1E3, bin_width);
    if year_index == 3
        [pv_pdf, pv_bincenter, pv_edges] = return_pdf(-pv{year_index}./1E3, bin_width);
    elseif year_index == 4
        [pv_pdf, pv_bincenter, pv_edges] = return_pdf(-stpv2019./1E3, bin_width);
    else
        [pv_pdf, pv_bincenter, pv_edges] = return_pdf(-pv{year_index}./1E3, bin_width);
        [st_pdf, st_bincenter, st_edges] = return_pdf(-st{year_index}./1E3, bin_width);
    end

    if year_index <= 2
        nl = loads{year_index} - wind{year_index} - st{year_index} - btm{year_index} - pv{year_index};
    elseif year_index == 3
        nl = loads{year_index} - wind{year_index} - btm{year_index} - pv{year_index};
    elseif year_index == 4
        nl = loads{year_index} - wind{year_index} - btm{year_index} - stpv2019;
    end
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
    [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [w_pdf(:); 0], w_edges(:), bin_width);
    [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [btm_pdf(:); 0], btm_edges(:), bin_width);
    if year_index <= 2
        [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [st_pdf(:); 0], st_edges(:), bin_width);
        [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [pv_pdf(:); 0], pv_edges(:), bin_width);
    else
        [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [pv_pdf(:); 0], pv_edges(:), bin_width);
    end
    
    ax = subplot(2, 2, year_index);
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');
end
suptitle('Net load');
end

function convolution_2component_no_corr_new()
% Convolution without correlation, but with the conv_poly_corr function, using CAISO's flexibility assessment dataset
load caiso;
all_years = [2016; 2017; 2018; 2019];
bin_width = 0.050; % GW
year_index = 1;

all_figures = nan(length(all_years), 1);
all_titles = {'Load - wind', 'Load - solar PV/solar', 'Load - solar thermal', 'Load - BTM'};
for i = 1: 4 % We have 4 components: PV, thermal, wind, and BTM
    all_figures(i) = figure();
end

for year_index = 1: 4

    [l_pdf, l_bincenter, l_edges] = return_pdf(loads{year_index}./1E3, bin_width);
    [w_pdf, w_bincenter, w_edges] = return_pdf(-wind{year_index}./1E3, bin_width);
    [btm_pdf, btm_bincenter, btm_edges] = return_pdf(-btm{year_index}./1E3, bin_width);
    if year_index == 3
        [pv_pdf, pv_bincenter, pv_edges] = return_pdf(-pv{year_index}./1E3, bin_width);
    elseif year_index == 4
        [pv_pdf, pv_bincenter, pv_edges] = return_pdf(-stpv2019./1E3, bin_width);
    else
        [pv_pdf, pv_bincenter, pv_edges] = return_pdf(-pv{year_index}./1E3, bin_width);
        [st_pdf, st_bincenter, st_edges] = return_pdf(-st{year_index}./1E3, bin_width);
    end

    % Just load and wind
    fig = figure(all_figures(1));
    subplot(2, 2, year_index);
    nl = loads{year_index} - wind{year_index};
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
%     [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [w_pdf(:); 0], w_edges(:));
    Ix = size(w_edges(:), 1);
    Iy = size(l_edges(:), 1);
    [a_convoluted, t_convoluted] = conv_poly_corr([l_pdf(:); 0], l_edges(:), [w_pdf(:); 0], w_edges(:), ones(Ix-1, Iy-1), bin_width); % Discretized copula does not inclue when CDF > 1.
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');

    % Just load and solar PV
    fig = figure(all_figures(2));
    subplot(2, 2, year_index);
    if year_index <= 3
        nl = loads{year_index} - pv{year_index};
    else
        nl = loads{year_index} - stpv2019;
    end
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
    Ix = size(pv_edges(:), 1);
    Iy = size(l_edges(:), 1);
    [a_convoluted, t_convoluted] = conv_poly_corr([l_pdf(:); 0], l_edges(:), [pv_pdf(:); 0], pv_edges(:), ones(Ix-1, Iy-1), bin_width);
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');

    % Just load and solar thermal
    fig = figure(all_figures(3));
    if year_index <= 2
        subplot(2, 2, year_index);
        nl = loads{year_index} - st{year_index};
        [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
        Ix = size(st_edges(:), 1);
        Iy = size(l_edges(:), 1);
        [a_convoluted, t_convoluted] = conv_poly_corr([l_pdf(:); 0], l_edges(:), [st_pdf(:); 0], st_edges(:), ones(Ix-1, Iy-1), bin_width);
        plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
        plot(nl_bincenter, nl_pdf, 'r');
        title(all_years(year_index));
        xlabel('GW');
    end

    % Just load and BTM data
    fig = figure(all_figures(4));
    subplot(2, 2, year_index);
    nl = loads{year_index} - btm{year_index};
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
    Ix = size(btm_edges(:), 1);
    Iy = size(l_edges(:), 1);
    [a_convoluted, t_convoluted] = conv_poly_corr([l_pdf(:); 0], l_edges(:), [btm_pdf(:); 0], btm_edges(:), ones(Ix-1, Iy-1), bin_width);
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');
end

for i = 1: length(all_figures)
    figure(all_figures(i));
    suptitle(all_titles{i});
end

end

function convolution_2component_corr(copula_type)
% Convolution with correlation, using CAISO's flexibility assessment dataset
load caiso;
all_years = [2016; 2017; 2018; 2019];
bin_width = 100; % MW
flag_subplot = 0; % 1 means panel graph
% copula_type = 'Gaussian';

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

all_figures = nan(length(all_years), 1);
all_titles = {'Load - wind', 'Load - solar PV/solar', 'Load - solar thermal', 'Load - BTM'};
% for i = 1: 4 % We have 4 components: PV, thermal, wind, and BTM
%     all_figures(i) = figure();
% end

for i = 1:4
%     if flag_subplot
%         fig = figure(i);
%     end
    for year_index = 1:1 % For all four years, change the 2nd 1 to 4
        % Select component 1 (x axis component) based on i
        switch i
            case 1 % Wind
                ts_1 = -wind{year_index};
            case 2 % Solar PV
                if year_index == 4
                    ts_1 = -stpv2019;
                else
                    ts_1 = -pv{year_index};
                end
            case 3 % Solar thermal
                if year_index > 2
                    continue;
                else
                    ts_1 = -st{year_index};
                end
            case 4 % BTM
                ts_1 = -btm{year_index};
        end
        ts_2 = loads{year_index}; % Y axis component
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

        i_valid = (ts_cdf_2<1) & (ts_cdf_2>0) & (ts_cdf_1<1) & (ts_cdf_2>0);
        Ix = size(gr_cdf_1(:), 1);
        Iy = size(gr_cdf_2(:), 1);

        switch copula_type
            case 'Gaussian' % Gaussian copula
                rhohat = copulafit('Gaussian', [ts_cdf_2(i_valid) ts_cdf_1(i_valid)]);
                [copula_x_grid, copula_y_grid] = meshgrid(gr_cdf_1, gr_cdf_2);
                copula_cdf = copulacdf('Gaussian', [copula_x_grid(:) copula_y_grid(:)], rhohat);
                copula_cdf = reshape(copula_cdf, Iy, Ix);
            case 'empirical' % Empirical copula
                cdfcounts = histcounts2(ts_cdf_1, ts_cdf_2, gr_cdf_1, gr_cdf_2);
                cdfcounts = cdfcounts'; % Matlab uses different x and y directions
%             case 't'
%                 [rhohat, nu] = copulafit('t', [ts_cdf_2(i_valid) ts_cdf_1(i_valid)]);
%                 [copula_x_grid, copula_y_grid] = meshgrid(gr_cdf_1, gr_cdf_2);
%                 copula_cdf = copulacdf('t', [copula_x_grid(:) copula_y_grid(:)], rhohat, nu);
%                 copula_cdf = reshape(copula_cdf, Iy, Ix);
            otherwise % Archimedean copulas and t copula, param2 is the 2nd parameter
                [rhohat, param2] = copulafit(copula_type, [ts_cdf_2(i_valid) ts_cdf_1(i_valid)]);
                [copula_x_grid, copula_y_grid] = meshgrid(gr_cdf_1, gr_cdf_2);
                copula_cdf = copulacdf(copula_type, [copula_x_grid(:) copula_y_grid(:)], rhohat, param2);
                copula_cdf = reshape(copula_cdf, Iy, Ix);
        end
        
        mean_copula_pdf = nan(Iy-1, Ix-1);
        switch copula_type
            case 'empirical'
                for ii = 1: Iy-1
                    for jj = 1: Ix-1
                        grid_area = (gr_cdf_2(ii+1) - gr_cdf_2(ii))*(gr_cdf_1(jj+1) - gr_cdf_1(jj)); % Note grid area is not homogeneous!
                        mean_copula_pdf(ii,jj) = cdfcounts(ii,jj)/(grid_area*sum(cdfcounts(:)));
                    end
                end
            otherwise
                for ii = 1: Iy-1
                    for jj = 1: Ix-1
                        grid_area = (gr_cdf_2(ii+1) - gr_cdf_2(ii))*(gr_cdf_1(jj+1) - gr_cdf_1(jj)); % Note grid area is not homogeneous!
                        mean_copula_pdf(ii,jj) = (copula_cdf(ii+1, jj+1) + copula_cdf(ii, jj) - copula_cdf(ii+1, jj) - copula_cdf(ii, jj+1))/grid_area;
                    end
                end
        end

        [a_convcorr, t_convcorr] = conv_poly_corr([gr_pdf_1(:); 0], gr_edges_1(:), [gr_pdf_2(:); 0], gr_edges_2(:), mean_copula_pdf, bin_width);
        [a_conv, t_conv] = conv_poly([gr_pdf_1(:); 0], gr_edges_1(:), [gr_pdf_2(:); 0], gr_edges_2(:), bin_width);
        
%         if flag_subplot
%             subplot(2, 2, year_index);
%         else
%             fig = figure();
%         end
        fig_pdf = figure();
%         h3 = plot(gr_bincenter_3, gr_pdf_3, 'color', 'r', 'linewidth', 1);
        h3 = bar(gr_bincenter_3, gr_pdf_3, 'FaceColor', [.8, .8, .8], 'EdgeColor', [.8, .8, .8], 'FaceAlpha', .5);
        h1 = plot_conv_poly(fig_pdf, a_conv, t_conv, 'b');
        h2 = plot_conv_poly(fig_pdf, a_convcorr, t_convcorr, 'k');
        set(h1, 'LineWidth', 1);
        set(h2, 'LineWidth', 1);
        legend([h1, h2, h3], 'No copula', 'Copula', 'Actual');
        legend boxoff;
%         title(all_years(year_index));
        xlabel('GW');
        ylabel('Probability density (GW^{-1})');
        box on;
        set(findall(fig_pdf,'-property','FontSize'),'FontSize',20);

        % Prepare for PP and CoD
        gr_cdf_3 = [0 cumsum(gr_pdf_3.*bin_width)];
        cdf_conv = cdf_poly(a_conv,t_conv);
        cdf_convcorr = cdf_poly(a_convcorr, t_convcorr);

        % PP plot
        fig_pp = figure();
        i_conv_start = find(abs(t_conv(:) - gr_edges_3(1)) < bin_width/10 );
        i_conv_end   = find(abs(t_conv(:) - gr_edges_3(end)) < bin_width/10 );
        h12 = plot(gr_cdf_3, cdf_conv(i_conv_start:i_conv_end), '-b.', gr_cdf_3, cdf_convcorr(i_conv_start:i_conv_end), '-k.');
        hold on;
        h3 = plot([0, 1], [0, 1], '-r', 'LineWidth', 1);
        legend([h12; h3], 'No copula', 'Copula', 'Diagonal');
        legend boxoff;
        xlabel('Empirical CDF');
        ylabel('Convolution CDF');
        box on;
        set(findall(fig_pp,'-property','FontSize'),'FontSize',20);

        % Coefficients of determination
        mdl_conv = fitlm(gr_cdf_3, cdf_conv(i_conv_start:i_conv_end));
        mdl_convcorr = fitlm(gr_cdf_3, cdf_convcorr(i_conv_start:i_conv_end));
        
        % Chi squared calculation
        ndegree = numel(ts_3) - 1; % Degree of the chi2 distribution
        f_conv = (cdf_conv(2:end) - cdf_conv(1:end-1)).*numel(ts_3); % Frequency
        f_convcorr = (cdf_convcorr(2:end) - cdf_convcorr(1:end-1)).*numel(ts_3); % Frequency
        f_observe = (gr_cdf_3(2:end) - gr_cdf_3(1:end-1)).*numel(ts_3);
        f_observe_complete = [zeros(i_conv_start-1, 1); f_observe(:); zeros(numel(f_conv) - i_conv_end + 1 ,1)]; % Frequency, including missing points
        chi2_conv = sum((f_conv - f_observe_complete).^2./f_conv);
        chi2_convcorr = nansum((f_convcorr - f_observe_complete).^2./f_convcorr); % Ignore nan due to 0 in the denominator
        chi2_critical = chi2inv(0.95, ndegree);
        fprintf('Component: (%g, %g), chi2: %.1f (n-1: %g, 0.05 reject: %.1f), R^2   No correlation: %f\n', i,  all_years(year_index), chi2_conv, ndegree, chi2_critical, mdl_conv.Rsquared.Adjusted);
        fprintf('Component: (%g, %g), chi2: %.1f (n-1: %g, 0.05 reject: %.1f), R^2 With correlation: %f\n', i,  all_years(year_index), chi2_convcorr, ndegree, chi2_critical, mdl_convcorr.Rsquared.Adjusted);
    end
%     if flag_subplot
%         figure(i);
%         suptitle(all_titles{i});
%     end
end
end

function convolution_netload_corr(copula_type)
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

for year_index = 1:1 % Change to 4 if want to run all four years
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
            a_prev = [gr_pdf_1(:); 0];
            t_prev = gr_edges_1(:);
            a_prev_corr = a_prev;
            t_prev_corr = t_prev;
            
            gr_cdf_forconv_1 = gr_cdf_1;
            gr_cdf_forconv_2 = gr_cdf_2;
        else
            a_prev = a_conv;
            t_prev = t_conv;
            a_prev_corr = a_convcorr;
            t_prev_corr = t_convcorr;
            
            gr_cdf_forconv_1 = cdf_poly(a_prev_corr, t_prev_corr);
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
        
        tic; 
        [a_convcorr, t_convcorr] = conv_poly_corr(a_prev_corr, t_prev_corr, [gr_pdf_2(:); 0], gr_edges_2(:), mean_copula_pdf, bin_width);
        [a_convcorr, t_convcorr] = cleanup_poly(a_convcorr, t_convcorr); % Clean up leading (0) and trailing (1) constant terms
        time_convcorr = time_convcorr + toc; 
        tic;
        [a_conv, t_conv] = conv_poly(a_prev, t_prev, [gr_pdf_2(:); 0], gr_edges_2(:), bin_width);
        time_conv = time_conv + toc; 
    end
%     h3 = plot(gr_bincenter_3, gr_pdf_3, 'color', 'r', 'linewidth', 1);
    h3 = bar(gr_bincenter_3, gr_pdf_3, 'FaceColor', [.8, .8, .8], 'EdgeColor', [.8, .8, .8], 'FaceAlpha', .5);
    h1 = plot_conv_poly(fig, a_conv, t_conv, 'b');
    h2 = plot_conv_poly(fig, a_convcorr, t_convcorr, 'k');
    set(h1, 'LineWidth', 1);
    set(h2, 'LineWidth', 1);
    legend([h1, h2, h3], 'No copula', 'Copula', 'Actual');
%     title(all_years(year_index));
    xlabel('GW');
    legend boxoff;
    ylabel('PDF (GW^{-1})');
    box on;
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    ax = gca;
%     ax.YAxis.Exponent = -3;
    
    % Prepare results for graphs
    gr_cdf_3 = [0 cumsum(gr_pdf_3.*bin_width)];
    cdf_conv = cdf_poly(a_conv,t_conv);
    cdf_convcorr = cdf_poly(a_convcorr, t_convcorr);

    % CDF comparison
    figure();
    plot(t_conv, cdf_conv, t_convcorr, cdf_convcorr, gr_edges_3, gr_cdf_3);
    legend('No correlation', 'With correlation', 'Raw');

    % PP plot
    figure();
%     i_conv_start = find(t_conv(:) == gr_edges_3(1) );
%     i_conv_end   = find(t_conv(:) == gr_edges_3(end) );
    i_conv_start = find(abs(t_conv(:) - gr_edges_3(1)) < bin_width/10 );
    i_conv_end   = find(abs(t_conv(:) - gr_edges_3(end)) < bin_width/10 );
    i_convcorr_start = find(abs(t_convcorr(:) - gr_edges_3(1)) < bin_width/10 );
    i_convcorr_end   = find(abs(t_convcorr(:) - gr_edges_3(end)) < bin_width/10 );

    h12 = plot(gr_cdf_3, cdf_conv(i_conv_start:i_conv_end), '-b.', gr_cdf_3, cdf_convcorr(i_convcorr_start:i_convcorr_end), '-k.');
    hold on;
    h3 = plot([0, 1], [0, 1], '-r', 'LineWidth', 1);
    legend([h12; h3], 'No correlation', 'With correlation', 'Diagonal');
    xlabel('Empirical CDF');
    ylabel('Convolution CDF');
    legend boxoff;
    box on;
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    
    % Coefficients of determination
    mdl_conv = fitlm(gr_cdf_3, cdf_conv(i_conv_start:i_conv_end));
    mdl_convcorr = fitlm(gr_cdf_3, cdf_convcorr(i_convcorr_start:i_convcorr_end));
%     fprintf('R^2   No correlation: %f\n', mdl_conv.Rsquared.Adjusted);
%     fprintf('R^2 With correlation: %f\n', mdl_convcorr.Rsquared.Adjusted);
    
    % Chi squared calculation
    ndegree = numel(ts_3) - 1; % Degree of the chi2 distribution
    f_conv = (cdf_conv(2:end) - cdf_conv(1:end-1)).*numel(ts_3); % Frequency
    f_convcorr = (cdf_convcorr(2:end) - cdf_convcorr(1:end-1)).*numel(ts_3); % Frequency
    f_observe = (gr_cdf_3(2:end) - gr_cdf_3(1:end-1)).*numel(ts_3);
    f_observe_complete = [zeros(i_conv_start-1, 1); f_observe(:); zeros(numel(f_conv) - i_conv_end + 1 ,1)]; % Frequency, including missing points
    chi2_conv = sum((f_conv - f_observe_complete).^2./f_conv, 'omitnan');
    f_observe_complete_corr = [zeros(i_convcorr_start-1, 1); f_observe(:); zeros(numel(f_convcorr) - i_convcorr_end + 1 ,1)]; % Frequency, including missing points
    chi2_convcorr = sum((f_convcorr - f_observe_complete_corr).^2./f_convcorr, 'omitnan');
    chi2_critical = chi2inv(0.95, ndegree);
    fprintf('Component: (%g, %g), time: %5.1f, chi2: %.1f (n-1: %g, 0.05 reject: %.1f), R^2   No correlation: %f\n', i,  all_years(year_index), time_conv, chi2_conv, ndegree, chi2_critical, mdl_conv.Rsquared.Adjusted);
    fprintf('Component: (%g, %g), time: %5.1f, chi2: %.1f (n-1: %g, 0.05 reject: %.1f), R^2 With correlation: %f\n', i,  all_years(year_index), time_convcorr, chi2_convcorr, ndegree, chi2_critical, mdl_convcorr.Rsquared.Adjusted);
end
end
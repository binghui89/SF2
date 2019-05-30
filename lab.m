%% The Irwin–Hall distribution when n = 1 to 5
a0 = [1; 0;]; % Size: I by M, I is number of segments + 1 (last segment is always 0), M is the highest order + 1 (lowest order is always 0)
t0 = [0; 1;]; % Size: I by 1
a = a0;
t = t0;

fig = figure();
for i = 1: 5
    plot_conv_poly(fig, a, t, rand(1,3));
    [a_convoluted, t_convoluted] = conv_poly(a0, t0, a, t);
    a = a_convoluted; t = t_convoluted;
end
%% The paper example
c1 = [1; -1;]; t1 = [-3;-2;];
c2 = [0 2; 0 -4; 0 2; 0 2; 0 -4; 0 2; 0 2; 0 -4; 0 2;]; 
t2 = [0;0.5;1;2;2.5;3;4;4.5;5];

a1 = a_by_c_and_t(c1, t1); a2 = a_by_c_and_t(c2, t2);
[a_convoluted, t_convoluted] = conv_poly(a1, t1, a2, t2);
% a_convoluted = a_by_c_and_t(c_convoluted, t_convoluted);

fig = figure(); hold on;
plot_conv_poly(fig, a1, t1, 'b'); plot_conv_poly(fig, a2, t2, 'b');
plot_conv_poly(fig, a_convoluted, t_convoluted, 'k');

%% The CAISO data: Data validation
all_years = [2016; 2017; 2018; 2019];

figure();
for i = 1: length(load)
    subplot(2, 2, i);
    plot(load{i, 1}./1E3);
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

%% The CAISO data, demonstration of convolution
all_years = [2016; 2017; 2018; 2019];
bin_width = 0.050; % GW
year_index = 1;

all_figures = nan(length(all_years), 1);
all_titles = {'Load - wind', 'Load - solar PV/solar', 'Load - solar thermal', 'Load - BTM'};
for i = 1: 4 % We have 4 components: PV, thermal, wind, and BTM
    all_figures(i) = figure();
end

for year_index = 1: 4

    [l_pdf, l_bincenter, l_edges] = return_pdf(load{year_index}./1E3, bin_width);
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
    nl = load{year_index} - wind{year_index};
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
    [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [w_pdf(:); 0], w_edges(:));
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');

    % Just load and solar PV
    fig = figure(all_figures(2));
    subplot(2, 2, year_index);
    if year_index <= 3
        nl = load{year_index} - pv{year_index};
    else
        nl = load{year_index} - stpv2019;
    end
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
    [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [pv_pdf(:); 0], pv_edges(:));
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');

    % Just load and solar thermal
    fig = figure(all_figures(3));
    if year_index <= 2
        subplot(2, 2, year_index);
        nl = load{year_index} - st{year_index};
        [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
        [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [st_pdf(:); 0], st_edges(:));
        plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
        plot(nl_bincenter, nl_pdf, 'r');
        title(all_years(year_index));
        xlabel('GW');
    end

    % Just load and BTM data
    fig = figure(all_figures(4));
    subplot(2, 2, year_index);
    nl = load{year_index} - btm{year_index};
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
    [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [btm_pdf(:); 0], btm_edges(:));
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');
end

for i = 1: length(all_figures)
    figure(all_figures(i));
    suptitle(all_titles{i});
end

% Net load
fig = figure();
for year_index = 1: 4
    [l_pdf, l_bincenter, l_edges] = return_pdf(load{year_index}./1E3, bin_width);
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
        nl = load{year_index} - wind{year_index} - st{year_index} - btm{year_index} - pv{year_index};
    elseif year_index == 3
        nl = load{year_index} - wind{year_index} - btm{year_index} - pv{year_index};
    elseif year_index == 4
        nl = load{year_index} - wind{year_index} - btm{year_index} - stpv2019;
    end
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
    [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [w_pdf(:); 0], w_edges(:));
    [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [btm_pdf(:); 0], btm_edges(:));
    if year_index <= 2
        [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [st_pdf(:); 0], st_edges(:));
        [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [pv_pdf(:); 0], pv_edges(:));
    else
        [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [pv_pdf(:); 0], pv_edges(:));
    end
    
    ax = subplot(2, 2, year_index);
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');
end
suptitle('Net load');

%% "Convolution" with correlation using triangular distributions
copula_discrete = [0.3 1; 1 1.7;];

a = [0 1; 2 -1; 0 0;];
t = [0; 1; 2;];

% No correlation
[a_convoluted, t_convoluted] = conv_poly(a, t, a, t);
fig = figure();
plot_conv_poly(fig, a_convoluted, t_convoluted, 'r');

% With correlation
% deltat = unique(diff(t)); % Samples must be equally placed
% tc_1 = t(1)+t(1);
% tc_e = t(end)+t(end);
% K = round((tc_e - tc_1)/deltat);
% 
% t_convcorr = nan(K+1, 1);
% a_convcorr = nan(K+1, size(a, 2) + size(a, 2));
% t_convcorr(1) = tc_1;
% a_convcorr(end, :) = 0;
% 
% for k = 1: K
%     [cx, cy] = c_vector(copula_discrete, k);
%     ax = repmat([cx; 0], 1, size(a, 2)).*a;
%     ay = repmat([cy; 0], 1, size(a, 2)).*a;
%     [a_tmp, t_tmp] = conv_poly(ax, t, ay, t);
%     a_convcorr(k, :) = a_tmp(t_tmp==t_convcorr(k), :);
%     t_convcorr(k+1) = tc_1+k*deltat;
% end
[a_convcorr, t_convcorr] = conv_poly_corr(a, t, a, t, copula_discrete);
% fig = figure();
plot_conv_poly(fig, a_convcorr, t_convcorr, 'k');

%% "Convolution" with correlation using CAISO's flexibility assessment dataset
all_years = [2016; 2017; 2018; 2019];
bin_width = 0.050; % GW
year_index = 1;

all_figures = nan(length(all_years), 1);
all_titles = {'Load - wind', 'Load - solar PV/solar', 'Load - solar thermal', 'Load - BTM'};
for i = 1: 4 % We have 4 components: PV, thermal, wind, and BTM
    all_figures(i) = figure();
end

for year_index = 1: 4

    [l_pdf, l_bincenter, l_edges] = return_pdf(load{year_index}./1E3, bin_width);
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
    nl = load{year_index} - wind{year_index};
    [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
%     [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [w_pdf(:); 0], w_edges(:));
    Ix = size(w_edges(:), 1);
    Iy = size(l_edges(:), 1);
    [a_convoluted, t_convoluted] = conv_poly_corr([l_pdf(:); 0], l_edges(:), [w_pdf(:); 0], w_edges(:), ones(Ix-1, Iy-1)); % Discretized copula does not inclue when CDF > 1.
    plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
    plot(nl_bincenter, nl_pdf, 'r');
    title(all_years(year_index));
    xlabel('GW');

%     % Just load and solar PV
%     fig = figure(all_figures(2));
%     subplot(2, 2, year_index);
%     if year_index <= 3
%         nl = load{year_index} - pv{year_index};
%     else
%         nl = load{year_index} - stpv2019;
%     end
%     [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
%     [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [pv_pdf(:); 0], pv_edges(:));
%     plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
%     plot(nl_bincenter, nl_pdf, 'r');
%     title(all_years(year_index));
%     xlabel('GW');
% 
%     % Just load and solar thermal
%     fig = figure(all_figures(3));
%     if year_index <= 2
%         subplot(2, 2, year_index);
%         nl = load{year_index} - st{year_index};
%         [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
%         [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [st_pdf(:); 0], st_edges(:));
%         plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
%         plot(nl_bincenter, nl_pdf, 'r');
%         title(all_years(year_index));
%         xlabel('GW');
%     end
% 
%     % Just load and BTM data
%     fig = figure(all_figures(4));
%     subplot(2, 2, year_index);
%     nl = load{year_index} - btm{year_index};
%     [nl_pdf, nl_bincenter, nl_edges] = return_pdf(nl./1E3, bin_width);
%     [a_convoluted, t_convoluted] = conv_poly([l_pdf(:); 0], l_edges(:), [btm_pdf(:); 0], btm_edges(:));
%     plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
%     plot(nl_bincenter, nl_pdf, 'r');
%     title(all_years(year_index));
%     xlabel('GW');
end

for i = 1: length(all_figures)
    figure(all_figures(i));
    suptitle(all_titles{i});
end

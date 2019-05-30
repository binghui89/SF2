%% Basic data
all_years = [2016; 2017; 2018; 2019];

%% Load, 2016, 2017, 2018, 2019
figure(1);
hold on;
for i = 1:4
    x = load{i};
    [counts, binedges] = histcounts(x./1000, 50);
    bincenter = (binedges(1: end-1) + binedges(2: end))/2;
    binwidth = binedges(2)-binedges(1);
    plot(bincenter, counts./(sum(counts)*binwidth), 'LineWidth', 1)
end
hold off;
box on;
xlabel('Load (GW)');
ylabel('Probability density');
legend('2016', '2017', '2018', '2019');

%% Wind, 2016, 2017, 2018, 2019
figure(2);
hold on;
for i = 1:4
    x = wind{i};
    [counts, binedges] = histcounts(x, 50);
    bincenter = (binedges(1: end-1) + binedges(2: end))/2;
    binwidth = binedges(2)-binedges(1);
    plot(bincenter, counts./(sum(counts)*binwidth), 'LineWidth', 1)
end
hold off;
box on;
xlabel('Wind power (MW)');
ylabel('Probability density');
legend('2016', '2017', '2018', '2019');
set(gca, 'YScale', 'log')

%% solar pv, 2016, 2017, 2018, 2019(?)
figure(3);
hold on;
for i = 1:4
    if i == 4
        x = stpv2019;
    else
        x = pv{i};
    end
    [counts, binedges] = histcounts(x, 50);
    bincenter = (binedges(1: end-1) + binedges(2: end))/2;
    binwidth = binedges(2)-binedges(1);
    plot(bincenter, counts./(sum(counts)*binwidth), 'LineWidth', 2)
end
hold off;
box on;
xlabel('Solar PV (MW)');
ylabel('Probability density');
legend('2016', '2017', '2018', '2019');
set(gca, 'YScale', 'log')

%% Solar thermal, 2016, 2017, 2019(?)
figure(5);
hold on;
for i = 1:2
    if i == 4
        x = stpv2019;
    elseif i == 3
        continue
    else
        x = st{i};
    end
    [counts, binedges] = histcounts(x, 50);
    bincenter = (binedges(1: end-1) + binedges(2: end))/2;
    binwidth = binedges(2)-binedges(1);
    plot(bincenter, counts./(sum(counts)*binwidth), 'LineWidth', 2)
end
hold off;
box on;
xlabel('Solar thermal (MW)');
ylabel('Probability density');
legend('2016', '2017', '2019'); % 2018 solar thermal data is missing
set(gca, 'YScale', 'log')

%% Behind-the-meter, 2016, 2017, 2018, 2019
figure(4);
hold on;
for i = 1:4
    x = btm{i};
    [counts, binedges] = histcounts(x, 50);
    bincenter = (binedges(1: end-1) + binedges(2: end))/2;
    binwidth = binedges(2)-binedges(1);
    plot(bincenter, counts./(sum(counts)*binwidth), 'LineWidth', 2)
end
hold off;
box on;
xlabel('Behind-the-meter (MW)');
ylabel('Probability density');
legend('2016', '2017', '2018', '2019');
set(gca, 'YScale', 'log')

%% Try convolution, using 2016 data only, load and wind

bin_width = 0.050; % GW
year_index = 2;

figure()
for year_index = 1: 4
    l_min = floor(min(load{year_index}./1E3)./bin_width)*bin_width;
    l_max = ceil(max(load{year_index}./1E3)./bin_width)*bin_width;
    [l_counts, l_edges] = histcounts(load{year_index}./1E3, l_min:bin_width:l_max);
    l_bincenter = (l_edges(1: end-1) + l_edges(2: end))/2;

    neg_wind = -wind{year_index};
    w_min = floor(min(neg_wind./1E3)./bin_width)*bin_width;
    w_max = ceil(max(neg_wind./1E3)./bin_width)*bin_width;
    [w_counts, w_edges] = histcounts(neg_wind./1E3, w_min:bin_width:w_max);
    w_bincenter = (w_edges(1: end-1) + w_edges(2: end))/2;

    lmw = load{year_index}./1E3 - wind{year_index}./1E3;
    lmw_min = floor(min(lmw)./bin_width)*bin_width;
    lmw_max = ceil(max(lmw)./bin_width)*bin_width;
    [lmw_counts, lmw_edges] = histcounts(lmw, lmw_min:bin_width:lmw_max);
    lmw_bincenter = (lmw_edges(1: end-1) + lmw_edges(2: end))/2;

    [z, z_pdf] = my_convolution(load{year_index}./1E3, -wind{year_index}./1E3, bin_width);

    % figure();
    % plot(l_bincenter, l_counts./sum(l_counts), -w_bincenter, flipud(w_counts./sum(w_counts)), 'LineWidth', 1);
    % legend('Load', 'Wind');
    % xlabel('Power (GW)');
    % ylabel('Probability');

%     figure();
    subplot(2, 2, year_index);
    plot(z, z_pdf, 'Color', 'k', 'LineWidth', 1);
    hold on;
    plot(lmw_bincenter, lmw_counts./sum(lmw_counts), 'Color', 'r', 'LineWidth', 1);
    legend('Convolution', 'Actual');
    xlim([-10, 50]);
    xlabel('Load - wind (GW)');
    ylabel('Probability');
    title(all_years(year_index));
end
%% Try convolution, using 2016 data only, load and solar PV

bin_width = 0.050; % GW
year_index = 1;

figure();
for year_index = 1: 4
    l_min = floor(min(load{year_index}./1E3)./bin_width)*bin_width;
    l_max = ceil(max(load{year_index}./1E3)./bin_width)*bin_width;
    [l_counts, l_edges] = histcounts(load{year_index}./1E3, l_min:bin_width:l_max);
    l_bincenter = (l_edges(1: end-1) + l_edges(2: end))/2;

    if year_index == 4
        neg_pv = -stpv2019;
    else
        neg_pv = -pv{year_index};
    end
    w_min = floor(min(neg_pv./1E3)./bin_width)*bin_width;
    w_max = ceil(max(neg_pv./1E3)./bin_width)*bin_width;
    [w_counts, w_edges] = histcounts(neg_pv./1E3, w_min:bin_width:w_max);
    w_bincenter = (w_edges(1: end-1) + w_edges(2: end))/2;

    lmw = load{year_index}./1E3 + neg_pv./1E3;
    lmw_min = floor(min(lmw)./bin_width)*bin_width;
    lmw_max = ceil(max(lmw)./bin_width)*bin_width;
    [lmw_counts, lmw_edges] = histcounts(lmw, lmw_min:bin_width:lmw_max);
    lmw_bincenter = (lmw_edges(1: end-1) + lmw_edges(2: end))/2;

    [z, z_pdf] = my_convolution(load{year_index}./1E3, neg_pv./1E3, bin_width);

    % figure();
    % plot(l_bincenter, l_counts./sum(l_counts), -w_bincenter, flipud(w_counts./sum(w_counts)), 'LineWidth', 1);
    % legend('Load', 'PV');
    % xlabel('Power (GW)');
    % ylabel('Probability');

%     figure()
    subplot(2, 2, year_index);
    plot(z, z_pdf, 'Color', 'k', 'LineWidth', 1);
    hold on;
    plot(lmw_bincenter, lmw_counts./sum(lmw_counts), 'Color', 'r', 'LineWidth', 1);
    legend('Convolution', 'Actual');
    xlim([0, 60]);
    xlabel('Load - PV (GW)');
    ylabel('Probability');
    title(all_years(year_index));
end
% Load data and pre-processing
add_pvlib();

load breakpoint; % Output from raster_process_data.m
Time.UTCOffset = zeros(size(ar_datetime, 1), 1); % Because IBM uses UTC time, so utc offset is zero
Time.year   = ar_datetime.Year;
Time.month  = ar_datetime.Month;
Time.day    = ar_datetime.Day;
Time.hour   = ar_datetime.Hour;
Time.minute = ar_datetime.Minute;
Time.second = ar_datetime.Second;

load_caiso_data;
ar_datetime_hour = datetime(ar_datetime.Year, ar_datetime.Month, ar_datetime.Day, ar_datetime.Hour, 0, 0, 'TimeZone', 'UTC');
T_rtpd = T_rtpd((T_rtpd{:, 'HOUR_START'}>=min(ar_datetime))&(T_rtpd{:, 'HOUR_START'}<=max(ar_datetime)), :); % Only look at the hours with raster forecast
T_rtpd_forcong = T_rtpd(:, {'HOUR_START', 'error_max', 'error_min'});
T_rtpd_forcong = grpstats(T_rtpd_forcong, 'HOUR_START', {'max', 'min'}, 'DataVars', {'error_max', 'error_min'});
T_rtpd_forcong = T_rtpd_forcong(:, {'HOUR_START', 'max_error_max', 'min_error_min'});
T_rtpd_forcong.Properties.VariableNames{'max_error_max'} = 'FRU';
T_rtpd_forcong.Properties.VariableNames{'min_error_min'} = 'FRD';

%% Run baseline and kNN using sum of PVlib power
% Baseline run
karray = 5:5:40; % the number of nearest neighbors
ar_datetime_hour_unique = unique(ar_datetime_hour);
testhour_end = datetime(2020, 4, 1, 2, 0, 0, 'TimeZone', 'UTC');
testhour_start = datetime(2020, 3, 4, 9, 0, 0, 'TimeZone', 'UTC');
T_baseline_rtpd = array2table(...
    ar_datetime_hour_unique((ar_datetime_hour_unique>=testhour_start)&(ar_datetime_hour_unique<=testhour_end)), ...
    'VariableNames', ...
    {'HOUR_START'}); % Result container

fprintf('Baseline\n');
for k = karray
    T_baseline_rtpd{:, strcat('FRU_k', num2str(k))} = nan;
    T_baseline_rtpd{:, strcat('FRD_k', num2str(k))} = nan;
    for i = 1: size(T_baseline_rtpd, 1)
        this_date = datetime(T_baseline_rtpd.HOUR_START.Year(i), T_baseline_rtpd.HOUR_START.Month(i), T_baseline_rtpd.HOUR_START.Day(i), 'TimeZone', 'UTC');
        this_hour = T_baseline_rtpd.HOUR_START.Hour(i);
        selected_days = (T_rtpd.TIME_START<this_date)&(T_rtpd.TIME_START>=this_date-days(k))&(T_rtpd.HOUR_START.Hour==this_hour); % We use k previous days
        sample_error_max = T_rtpd{selected_days, 'error_max'};
        sample_error_min = T_rtpd{selected_days, 'error_min'};
        [f,x] = ecdf(sample_error_max(:));
        T_baseline_rtpd{i, strcat('FRU_k', num2str(k))} = interp1(f, x, 0.975);
        [f,x] = ecdf(sample_error_min(:));
        T_baseline_rtpd{i, strcat('FRD_k', num2str(k))} = interp1(f, x, 0.025);
    end

    fprintf('k = %g\n', k);
end

%% Calculate classifier based on kcs and normalize
cell_pvcellkcs = cell(size(cell_pvcellghi));
for i = 1: numel(cell_pvcellkcs)
    cell_pvcellkcs{i} = cell_pvcellghi{i}./ghi_cs_cell;
end
ar_hourstart = unique(datetime(ar_datetime.Year, ar_datetime.Month, ar_datetime.Day, ar_datetime.Hour, 0, 0, 'TimeZone', 'UTC'));

ar_pvcellkcs_p25 = cell_pvcellkcs{2};
ar_pvcellkcs_p50 = cell_pvcellkcs{3};
ar_pvcellkcs_p75 = cell_pvcellkcs{4};
for j = 1: numel(ar_hourstart)
    istart = 6*(j-1) + 1;
    iend   = 6*j;
    ar_pvcellkcs_mean(j, :) = mean(ar_pvcellkcs_p50(istart:iend, :), 1, 'omitnan');
    ar_pvcellkcs_std(j, :) = std(ar_pvcellkcs_p50(istart:iend, :), 0, 1, 'omitnan');
    ar_pvcellkcs_vrb(j, :) = sqrt(1/5.*sum(diff(ar_pvcellkcs_p50(istart:iend, :), 1, 1).^2, 1));
    ar_pvcellkcs_wmean(j, :) = mean(ar_pvcellkcs_p75(istart:iend, :) - ar_pvcellkcs_p25(istart:iend, :), 1, 'omitnan');
end

T_classifier_cellkcs = table(ar_hourstart, ar_pvcellkcs_mean, ar_pvcellkcs_std, ar_pvcellkcs_vrb, ar_pvcellkcs_wmean, 'VariableNames', {'HOUR_START', 'mean_kcs', 'std_kcs', 'vrb_kcs', 'wmean_kcs'});
T_classifier_cellkcs.DATE = datetime(T_classifier_cellkcs.HOUR_START.Year, T_classifier_cellkcs.HOUR_START.Month, T_classifier_cellkcs.HOUR_START.Day, 'TimeZone', 'UTC');

T_classifier_cellkcs_norm = T_classifier_cellkcs(:, {'HOUR_START', 'mean_kcs', 'std_kcs', 'vrb_kcs', 'wmean_kcs', 'DATE'});
T_classifier_cellkcs_norm{:, 'mean_kcs'}  = T_classifier_cellkcs{:, 'mean_kcs'}./repmat(max(T_classifier_cellkcs{:, 'mean_kcs'}, [], 1), height(T_classifier_cellkcs), 1);
T_classifier_cellkcs_norm{:, 'std_kcs'}   = T_classifier_cellkcs{:, 'std_kcs'}./repmat(max(T_classifier_cellkcs{:, 'std_kcs'}, [], 1), height(T_classifier_cellkcs), 1);
T_classifier_cellkcs_norm{:, 'vrb_kcs'}   = T_classifier_cellkcs{:, 'vrb_kcs'}./repmat(max(T_classifier_cellkcs{:, 'vrb_kcs'}, [], 1), height(T_classifier_cellkcs), 1);
T_classifier_cellkcs_norm{:, 'wmean_kcs'} = T_classifier_cellkcs{:, 'wmean_kcs'}./repmat(max(T_classifier_cellkcs{:, 'wmean_kcs'}, [], 1), height(T_classifier_cellkcs), 1);

wcap = T_pvcell.TotalCapacity(T_pvcell.within_boundary==1)./sum(T_pvcell.TotalCapacity(T_pvcell.within_boundary==1));

%% kNN 1-dim
fprintf('kNN\n');
cell_knnkcs_rtpd = cell(4, 1);
for c = 1:4
    T_knnkcs_rtpd = T_classifier_cellkcs_norm((ar_datetime_hour_unique>=testhour_start)&(ar_datetime_hour_unique<=testhour_end), :); % Result container
    switch c
        case 1
            T_classifier_cellkcs_norm{:, 'c'} = T_classifier_cellkcs_norm{:, 'mean_kcs'};
            T_knnkcs_rtpd{:, 'c'}    = T_knnkcs_rtpd{:, 'mean_kcs'};
        case 2
            T_classifier_cellkcs_norm{:, 'c'} = T_classifier_cellkcs_norm{:, 'std_kcs'};
            T_knnkcs_rtpd{:, 'c'}    = T_knnkcs_rtpd{:, 'std_kcs'};
        case 3
            T_classifier_cellkcs_norm{:, 'c'} = T_classifier_cellkcs_norm{:, 'vrb_kcs'};
            T_knnkcs_rtpd{:, 'c'}    = T_knnkcs_rtpd{:, 'vrb_kcs'};
        case 4
            T_classifier_cellkcs_norm{:, 'c'} = T_classifier_cellkcs_norm{:, 'wmean_kcs'};
            T_knnkcs_rtpd{:, 'c'}    = T_knnkcs_rtpd{:, 'wmean_kcs'};
    end
    for k = karray
        T_knnkcs_rtpd{:, strcat('FRU_k', num2str(k))} = nan;
        T_knnkcs_rtpd{:, strcat('FRD_k', num2str(k))} = nan;

        for i = 1: size(T_knnkcs_rtpd, 1)
            this_date = datetime(T_knnkcs_rtpd.HOUR_START.Year(i), T_knnkcs_rtpd.HOUR_START.Month(i), T_knnkcs_rtpd.HOUR_START.Day(i), 'TimeZone', 'UTC');
            this_hour = T_knnkcs_rtpd.HOUR_START.Hour(i);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Can use the days AFTER this_date if results are not ideal
            T_sample = T_classifier_cellkcs_norm((T_classifier_cellkcs_norm.HOUR_START.Hour==this_hour), :);
%             T_sample = T_classifier_cellkcs_norm((T_classifier_cellkcs_norm.DATE<this_date)&(T_classifier_cellkcs_norm.HOUR_START.Hour==this_hour), :);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if any(isnan(T_knnkcs_rtpd{i, 'c'}))
                T_sample_sorted = T_sample(T_sample.DATE>=this_date-days(k), :); % We use 30 previous days, i.e., baseline, if data is nan
                selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); % 30 the nearest days for baseline
            else
                T_sample.dist = abs(T_sample{:, 'c'}-T_knnkcs_rtpd{i, 'c'})*wcap; % Manhattan distance weighted by cell capacity
                T_sample_sorted = sortrows(T_sample, 'dist');
                selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); 
            end

            sample_error_max = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_max'};
            sample_error_min = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_min'};
            [f,x] = ecdf(sample_error_max(:));
            T_knnkcs_rtpd{i, strcat('FRU_k', num2str(k))} = interp1(f, x, 0.975);
            [f,x] = ecdf(sample_error_min(:));
            T_knnkcs_rtpd{i, strcat('FRD_k', num2str(k))} = interp1(f, x, 0.025);
        end
        fprintf('c = %g, k = %g\n', c, k);
    end
    cell_knnkcs_rtpd{c} = T_knnkcs_rtpd;
end

%% Process results and visualize
f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});
T_actual_error = [T_baseline_rtpd(:, 'HOUR_START') T_errormax_rtpd T_errormin_rtpd];
ar_actual_error_fru = T_actual_error{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}};
ar_actual_error_frd = T_actual_error{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}};

%% Process results and visualize: Calculate oversupply and p_shortage
% Baseline
for k = karray
    [oversupply_baseline_fru(karray==k), pshort_baseline_fru(karray==k)] = calculate_oversupply_and_pshort(ar_actual_error_fru, T_baseline_rtpd{:, strcat('FRU_k', num2str(k))});
    [oversupply_baseline_frd(karray==k), pshort_baseline_frd(karray==k)] = calculate_oversupply_and_pshort(-ar_actual_error_frd, -T_baseline_rtpd{:, strcat('FRD_k', num2str(k))});
end

% kNN, 1-dim
for c = 1:4
    T_knnkcs_rtpd = cell_knnkcs_rtpd{c};
    for k = karray
        [oversupply_knnkcs_fru(karray==k, c), pshort_knnkcs_fru(karray==k, c)] = calculate_oversupply_and_pshort(ar_actual_error_fru, T_knnkcs_rtpd{:, strcat('FRU_k', num2str(k))});
        [oversupply_knnkcs_frd(karray==k, c), pshort_knnkcs_frd(karray==k, c)] = calculate_oversupply_and_pshort(-ar_actual_error_frd, -T_knnkcs_rtpd{:, strcat('FRD_k', num2str(k))});
    end
end

%% Produce figures: Trends as a function of K
figure(); 
h = plot(karray, oversupply_knnkcs_frd./1E3, karray, oversupply_baseline_frd./1E3);
set(h(1), 'Color', 'r', 'Marker', '.');
set(h(2), 'Color', 'g', 'Marker', '.');
set(h(3), 'Color', 'm', 'Marker', '.');
set(h(4), 'Color', 'b', 'Marker', '.');
set(h(5), 'Color', 'k', 'Marker', '.');
title('Oversupply (FRD)');
xlabel('K');
ylabel('Oversupply (GWh)');
legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');

figure(); 
h = plot(karray, pshort_knnkcs_frd, karray, pshort_baseline_frd); 
set(h(1), 'Color', 'r', 'Marker', '.');
set(h(2), 'Color', 'g', 'Marker', '.');
set(h(3), 'Color', 'm', 'Marker', '.');
set(h(4), 'Color', 'b', 'Marker', '.');
set(h(5), 'Color', 'k', 'Marker', '.');
title('Pshort (FRD)');
xlabel('K');
ylabel('Probability of shortage');
legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');

figure(); 
h = plot(karray, oversupply_knnkcs_fru./1E3, karray, oversupply_baseline_fru./1E3); 
set(h(1), 'Color', 'r', 'Marker', '.');
set(h(2), 'Color', 'g', 'Marker', '.');
set(h(3), 'Color', 'm', 'Marker', '.');
set(h(4), 'Color', 'b', 'Marker', '.');
set(h(5), 'Color', 'k', 'Marker', '.');
title('Oversupply (FRU)');
xlabel('K');
ylabel('Oversupply (GWh)');
legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');

figure(); 
h = plot(karray, pshort_knnkcs_fru, karray, pshort_baseline_fru); 
set(h(1), 'Color', 'r', 'Marker', '.');
set(h(2), 'Color', 'g', 'Marker', '.');
set(h(3), 'Color', 'm', 'Marker', '.');
set(h(4), 'Color', 'b', 'Marker', '.');
set(h(5), 'Color', 'k', 'Marker', '.');
title('Pshort (FRU)');
xlabel('K');
ylabel('Probability of shortage');
legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');

%% Tile layout
figure(); 
t = tiledlayout(2, 2);
ax1 = nexttile;
h = plot(karray, oversupply_knnkcs_frd./1E3, karray, oversupply_baseline_frd./1E3);
set(h(1), 'Color', 'r', 'Marker', '.');
set(h(2), 'Color', 'g', 'Marker', '.');
set(h(3), 'Color', 'm', 'Marker', '.');
set(h(4), 'Color', 'b', 'Marker', '.');
set(h(5), 'Color', 'k', 'Marker', '.');
% title('Oversupply (FRD)');
% xlabel('K');
ylabel('Oversupply (GWh)');
% legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
xline(30, 'Color', 'r', 'LineStyle', '--')
ylim([50, 250]);
grid on;

% figure(); 
ax2 = nexttile;
h = plot(karray, oversupply_knnkcs_fru./1E3, karray, oversupply_baseline_fru./1E3); 
set(h(1), 'Color', 'r', 'Marker', '.');
set(h(2), 'Color', 'g', 'Marker', '.');
set(h(3), 'Color', 'm', 'Marker', '.');
set(h(4), 'Color', 'b', 'Marker', '.');
set(h(5), 'Color', 'k', 'Marker', '.');
% title('Oversupply (FRU)');
% xlabel('K');
% ylabel('Oversupply (GWh)');
% legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
xline(30, 'Color', 'r', 'LineStyle', '--')
ylim([50, 250]);
grid on;

% figure(); 
ax3 = nexttile;
h = plot(karray, pshort_knnkcs_frd, karray, pshort_baseline_frd); 
set(h(1), 'Color', 'r', 'Marker', '.');
set(h(2), 'Color', 'g', 'Marker', '.');
set(h(3), 'Color', 'm', 'Marker', '.');
set(h(4), 'Color', 'b', 'Marker', '.');
set(h(5), 'Color', 'k', 'Marker', '.');
% title('Pshort (FRD)');
% xlabel('K');
ylabel('Prob of shortage');
% legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
ylim([0, .1]);
xline(30, 'Color', 'r', 'LineStyle', '--')
grid on;

% figure(); 
ax4 = nexttile;
h = plot(karray, pshort_knnkcs_fru, karray, pshort_baseline_fru); 
set(h(1), 'Color', 'r', 'Marker', '.');
set(h(2), 'Color', 'g', 'Marker', '.');
set(h(3), 'Color', 'm', 'Marker', '.');
set(h(4), 'Color', 'b', 'Marker', '.');
set(h(5), 'Color', 'k', 'Marker', '.');
% title('Pshort (FRU)');
% xlabel('K');
% ylabel('Probability of shortage');
% legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
ylim([0, .1]);
xline(30, 'Color', 'r', 'LineStyle', '--')
grid on;

xticklabels(ax1,{});
xticklabels(ax2,{});
yticklabels(ax2,{})
yticklabels(ax4,{})

xlabel(t,'K')
t.TileSpacing = 'compact';
t.Padding = 'none';
%% Produce figures: Pareto frontier, 1-dim

% This is the same figure scatter plots
% figure(); 
% h(5)= scatter(pshort_baseline_fru, oversupply_baseline_fru./1E3, 60, 'ko');
% hold on;
% for c = 1: 4
%     h(c) = scatter(pshort_knnkcs_fru(:, c), oversupply_knnkcs_fru(:, c)./1E3, 30, '^');
% end
% % h(6) = scatter(pshort_knnkpvndim_fru(:), oversupply_knnkpvndim_fru(:)./1E3, 30, '^');
% set(h(1), 'MarkerFaceColor', 'r');
% set(h(2), 'MarkerFaceColor', 'g');
% set(h(3), 'MarkerFaceColor', 'm');
% set(h(4), 'MarkerFaceColor', 'b');
% set(h(5), 'MarkerFaceColor', 'k');
% % set(h(6), 'MarkerFaceColor', 'y');
% xline(pshort_baseline_fru(karray==30), 'Color', 'r', 'LineStyle', '--');
% yline(oversupply_baseline_fru(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
% title('FRU');
% xlabel('Probability of shortage');
% ylabel('Oversupply (GWh)');
% legend(h, '\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
% box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the same figure with lines, not scatter plots
figure();
h(5)= line(pshort_baseline_fru, oversupply_baseline_fru./1E3);
hold on;
for c = 1: 4
    h(c) = line(pshort_knnkcs_fru(:, c), oversupply_knnkcs_fru(:, c)./1E3);
end
set(h(1), 'Color', 'r', 'MarkerFaceColor', 'r', 'Marker', '^', 'MarkerSize', 6);
set(h(2), 'Color', 'g', 'MarkerFaceColor', 'g', 'Marker', '^', 'MarkerSize', 6);
set(h(3), 'Color', 'm', 'MarkerFaceColor', 'm', 'Marker', '^', 'MarkerSize', 6);
set(h(4), 'Color', 'b', 'MarkerFaceColor', 'b', 'Marker', '^', 'MarkerSize', 6);
set(h(5), 'Color', 'k', 'MarkerFaceColor', 'k', 'Marker', 'o', 'MarkerSize', 6);
xline(pshort_baseline_fru(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(oversupply_baseline_fru(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
title('FRU');
xlabel('Probability of shortage');
ylabel('Oversupply (GWh)');
legend(h, '\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the same figure with scatter plots
% figure(); 
% h(5)= scatter(pshort_baseline_frd, oversupply_baseline_frd./1E3, 60, 'ko');
% hold on;
% for c = 1: 4
%     h(c) = scatter(pshort_knnkcs_frd(:, c), oversupply_knnkcs_frd(:, c)./1E3, 30, '^');
% end
% % h(6) = scatter(pshort_knnkpvndim_frd(:), oversupply_knnkpvndim_frd(:)./1E3, 30, '^');
% set(h(1), 'MarkerFaceColor', 'r');
% set(h(2), 'MarkerFaceColor', 'g');
% set(h(3), 'MarkerFaceColor', 'm');
% set(h(4), 'MarkerFaceColor', 'b');
% set(h(5), 'MarkerFaceColor', 'k');
% % set(h(6), 'MarkerFaceColor', 'y');
% xline(pshort_baseline_frd(karray==30), 'Color', 'r', 'LineStyle', '--');
% yline(oversupply_baseline_frd(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
% title('FRD');
% xlabel('Probability of shortage');
% ylabel('Oversupply (GWh)');
% legend(h, '\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
% box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the same figure with lines, not scatter plots
figure();
h(5)= line(pshort_baseline_frd, oversupply_baseline_frd./1E3);
hold on;
for c = 1: 4
    h(c) = line(pshort_knnkcs_frd(:, c), oversupply_knnkcs_frd(:, c)./1E3);
end
set(h(1), 'Color', 'r', 'MarkerFaceColor', 'r', 'Marker', '^', 'MarkerSize', 6);
set(h(2), 'Color', 'g', 'MarkerFaceColor', 'g', 'Marker', '^', 'MarkerSize', 6);
set(h(3), 'Color', 'm', 'MarkerFaceColor', 'm', 'Marker', '^', 'MarkerSize', 6);
set(h(4), 'Color', 'b', 'MarkerFaceColor', 'b', 'Marker', '^', 'MarkerSize', 6);
set(h(5), 'Color', 'k', 'MarkerFaceColor', 'k', 'Marker', 'o', 'MarkerSize', 6);
xline(pshort_baseline_frd(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(oversupply_baseline_frd(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
title('FRD');
xlabel('Probability of shortage');
ylabel('Oversupply (GWh)');
legend(h, '\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% knn 4-dim
fprintf('kNN (n-dim)\n');

w = []; % Weight, increment of 0.1 at each dimension
for i = 0: 10
    for j = 0:10-i
        for t = 0: 10-i-j
            k=10-i-j-t;
            w = [w, [i/10.*wcap; j/10.*wcap; t/10.*wcap; k/10.*wcap]];
        end
    end
end

cell_knnkcs_rtpd_ndim = cell(size(w, 2), 1);
name_classifiers = {'mean_kcs', 'std_kcs', 'vrb_kcs', 'wmean_kcs'};
for j = 1: size(w, 2)
    w_this = w(:, j);
    lw_this = true(size(w_this));
    lw_this(w_this==0) = false;
    T_knnkcs_rtpd_ndim = T_classifier_cellkcs_norm((ar_datetime_hour_unique>=testhour_start)&(ar_datetime_hour_unique<=testhour_end), :); % Result container
    for k = karray
        T_knnkcs_rtpd_ndim{:, strcat('FRU_k', num2str(k))} = nan;
        T_knnkcs_rtpd_ndim{:, strcat('FRD_k', num2str(k))} = nan;

        for i = 1: size(T_knnkcs_rtpd_ndim, 1)
            this_date = datetime(T_knnkcs_rtpd_ndim.HOUR_START.Year(i), T_knnkcs_rtpd_ndim.HOUR_START.Month(i), T_knnkcs_rtpd_ndim.HOUR_START.Day(i), 'TimeZone', 'UTC');
            this_hour = T_knnkcs_rtpd_ndim.HOUR_START.Hour(i);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Can use the days AFTER this_date if results are not ideal
            T_sample = T_classifier_cellkcs_norm((T_classifier_cellkcs_norm.HOUR_START.Hour==this_hour), :);
    %             T_sample = T_classifier_kpv((T_classifier_kpv.DATE<this_date)&(T_classifier_kpv.HOUR_START.Hour==this_hour), :);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if any(isnan(T_knnkcs_rtpd_ndim{i, {'mean_kpv', 'std_kpv', 'vrb_kpv'}}), 2)
            if any(isnan(T_knnkcs_rtpd_ndim{i, name_classifiers}(:, lw_this)), 2)
                T_sample_sorted = T_sample(T_sample.DATE>=this_date-days(k), :); % We use 30 previous days, i.e., baseline, if data is nan
                selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); % 30 the nearest days for baseline
            else
                T_sample.dist = abs(T_sample{:, name_classifiers}(:, lw_this)-T_knnkcs_rtpd_ndim{i, name_classifiers}(:, lw_this))*w(lw_this, j); % Manhattan distance
%                 T_sample.dist = abs(T_sample{:, {'mean_kpv', 'std_kpv', 'vrb_kpv'}}-T_knnkcs_rtpd_ndim{i, {'mean_kpv', 'std_kpv', 'vrb_kpv'}})*w(:, j); % Manhattan distance
                T_sample_sorted = sortrows(T_sample, 'dist');
                selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); 
            end

            sample_error_max = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_max'};
            sample_error_min = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_min'};
            [f,x] = ecdf(sample_error_max(:));
            T_knnkcs_rtpd_ndim{i, strcat('FRU_k', num2str(k))} = interp1(f, x, 0.975);
            [f,x] = ecdf(sample_error_min(:));
            T_knnkcs_rtpd_ndim{i, strcat('FRD_k', num2str(k))} = interp1(f, x, 0.025);
        end
        fprintf('j = %g, k = %g\n', j, k);
    end
    cell_knnkcs_rtpd_ndim{j, 1} = T_knnkcs_rtpd_ndim;
end

%% Process results and visualize: Calculate oversupply and p_shortage
% kNN, 4-dim
% Calculate kNN exhaustive search oversupply and pshort 
for j = 1: size(w, 2)
    T_knnkcs_rtpd_ndim = cell_knnkcs_rtpd_ndim{j, 1};
    for k = karray
        [oversupply_knnkcsndimopt_fru(j, karray==k), pshort_knnkcsndimopt_fru(j, karray==k)] = calculate_oversupply_and_pshort(ar_actual_error_fru, T_knnkcs_rtpd_ndim{:, strcat('FRU_k', num2str(k))});
        [oversupply_knnkcsndimopt_frd(j, karray==k), pshort_knnkcsndimopt_frd(j, karray==k)] = calculate_oversupply_and_pshort(-ar_actual_error_frd, -T_knnkcs_rtpd_ndim{:, strcat('FRD_k', num2str(k))});
    end
end


%% Color by dimension of classifiers
figure();
hold on;
h = [];
for ndim = 1:4
    pshort_frd_tmp = pshort_knnkcsndimopt_frd(sum(w==0, 1)/340==4-ndim, :);
    oversupply_frd_tmp = oversupply_knnkcsndimopt_frd(sum(w==0, 1)/340==4-ndim, :);
    h(ndim) = scatter(pshort_frd_tmp(:), oversupply_frd_tmp(:)./1e3, 'o');
end
h(numel(h)+1)= line(pshort_baseline_frd, oversupply_baseline_frd./1E3);
set(h(1), 'MarkerFaceColor', 'r', 'Marker', '^');
set(h(2), 'MarkerEdgeColor', 'b', 'Marker', '+');
set(h(3), 'MarkerFaceColor', 'm', 'Marker', '.');
set(h(4), 'MarkerFaceColor', 'g', 'Marker', '.');
set(h(5), 'Color', 'k', 'MarkerFaceColor', 'k', 'Marker', 'o', 'MarkerSize', 6);
xline(pshort_baseline_frd(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(oversupply_baseline_frd(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
legend(h, '1-dim', '2-dim', '3-dim', '4-dim', 'Baseline');
title('FRD');

figure();
hold on;
h = [];
for ndim = 1:4
    pshort_fru_tmp = pshort_knnkcsndimopt_fru(sum(w==0, 1)/340==4-ndim, :);
    oversupply_fru_tmp = oversupply_knnkcsndimopt_fru(sum(w==0, 1)/340==4-ndim, :);
    h(ndim) = scatter(pshort_fru_tmp(:), oversupply_fru_tmp(:)./1e3);
end
h(numel(h)+1)= line(pshort_baseline_fru, oversupply_baseline_fru./1E3);
xline(pshort_baseline_fru(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(oversupply_baseline_fru(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
legend(h, '1-dim', '2-dim', '3-dim', 'Baseline');
set(h(1), 'MarkerFaceColor', 'r', 'Marker', '^');
set(h(2), 'MarkerEdgeColor', 'b', 'Marker', '+');
set(h(3), 'MarkerFaceColor', 'm', 'Marker', '.');
set(h(4), 'MarkerFaceColor', 'g', 'Marker', '.');
% set(h(5), 'MarkerFaceColor', 'k');
set(h(5), 'Color', 'k', 'MarkerFaceColor', 'k', 'Marker', 'o', 'MarkerSize', 6);
legend(h, '1-dim', '2-dim', '3-dim', '4-dim', 'Baseline');
title('FRU');
%%
trainhour_end = datetime(2020, 3, 25, 2, 0, 0, 'TimeZone', 'UTC');
trainhour_start = datetime(2020, 3, 4, 9, 0, 0, 'TimeZone', 'UTC');

T_actual_error{:, 'istrain'} = 0;
T_actual_error{(T_actual_error.HOUR_START<=trainhour_end)&(T_actual_error.HOUR_START>=trainhour_start), 'istrain'} = 1;

for j = 1: size(w, 2)
    T_knnkcs_rtpd_ndim = cell_knnkcs_rtpd_ndim{j, 1};
    T_knnkcs_rtpd_ndim{:, 'istrain'} = 0;
    T_knnkcs_rtpd_ndim{(T_knnkcs_rtpd_ndim.HOUR_START<=trainhour_end)&(T_knnkcs_rtpd_ndim.HOUR_START>=trainhour_start), 'istrain'} = 1;
    cell_knnkcs_rtpd_ndim{j, 1} = T_knnkcs_rtpd_ndim;
end

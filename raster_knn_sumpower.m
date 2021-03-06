% Load data and pre-processing
add_pvlib();

load breakpoint;
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

% Calculate PV power using the simple way
t = 25 + ghi_cs_cell./800.*(45.5-20);
p1_solar_cs = ghi_cs_cell./1000.*(1+0.004.*(t - 25));
p1_solar = nan(size(ar_datetime, 1), size(cell_pvcellghi{1}, 2), 5);
for i = 1:5
    t = 25 + cell_pvcellghi{i}./800.*(45.5-20);
    p1_solar(:, :, i) = cell_pvcellghi{i}./1000.*(1+0.004.*(t - 25));
end
pvlib_flag = input('PVlib? (1 - Yes/0 - No): ');
if ~pvlib_flag
    p_solar_cs = p1_solar_cs;
    p_solar = p1_solar;
end

%% Run baseline and kNN using sum of PVlib power
% Baseline run
karray = 5:5:40; % the number of nearest neighbors
ar_datetime_hour_unique = unique(ar_datetime_hour);
testhour_end = datetime(2020, 4, 1, 2, 0, 0, 'TimeZone', 'UTC');
testhour_start = datetime(2020, 3, 25, 9, 0, 0, 'TimeZone', 'UTC');
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
        selected_days = (T_rtpd.TIME_START<this_date)&(T_rtpd.TIME_START>=this_date-days(k))&(T_rtpd.HOUR_START.Hour==this_hour); % We use 30 previous days
        sample_error_max = T_rtpd{selected_days, 'error_max'};
        sample_error_min = T_rtpd{selected_days, 'error_min'};
        [f,x] = ecdf(sample_error_max(:));
        T_baseline_rtpd{i, strcat('FRU_k', num2str(k))} = interp1(f, x, 0.975);
        [f,x] = ecdf(sample_error_min(:));
        T_baseline_rtpd{i, strcat('FRD_k', num2str(k))} = interp1(f, x, 0.025);
    end

    fprintf('k = %g\n', k);
end
plot(T_baseline_rtpd{:, 'HOUR_START'}, T_baseline_rtpd{:, 2:end}, 'b')
hold on;
plot(...
    T_rtpd_forcong{(T_rtpd_forcong{:, 'HOUR_START'}>=testhour_start)&(T_rtpd_forcong{:, 'HOUR_START'}<=testhour_end), 'HOUR_START'}, ...
    T_rtpd_forcong{(T_rtpd_forcong{:, 'HOUR_START'}>=testhour_start)&(T_rtpd_forcong{:, 'HOUR_START'}<=testhour_end), 2:end}, ...
    'r' ...
    );

%% Calculate classifier based on kpv
sigma_psolar = nan(size(p_solar, 1), 5);
% plot(sum(squeeze(p_solar(:, T_eia860_caiso{:, 'within_boundary'}==1, 3)*diag(T_eia860_caiso{T_eia860_caiso{:, 'within_boundary'}==1, 'NameplateCapacity_MW_'})), 2));
if pvlib_flag
    for q = 1:5
        sigma_psolar(:, q) = sum(squeeze(p_solar(:, T_eia860_caiso{:, 'within_boundary'}==1, 3)*diag(T_eia860_caiso{T_eia860_caiso{:, 'within_boundary'}==1, 'NameplateCapacity_MW_'})), 2);
    end
    sigma_psolarcs = sum(p_solar_cs(:, T_eia860_caiso{:, 'within_boundary'}==1)*diag(T_eia860_caiso{T_eia860_caiso{:, 'within_boundary'}==1, 'NameplateCapacity_MW_'}), 2);
else
    for q = 1:5
        sigma_psolar(:, q) = sum(squeeze(p_solar(:, :, 3)*diag(T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'TotalCapacity'})), 2);
    end
    sigma_psolarcs = sum(p_solar_cs*diag(T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'within_boundary'}), 2);
end

kpv = sigma_psolar./repmat(sigma_psolarcs, 1, 5);

T_classifier_kpv = array2table(ar_datetime_hour_unique, 'VariableNames', {'HOUR_START'});
T_classifier_kpv{:, 'mean_kpv'} = nan; % Mean of p50 kpv
T_classifier_kpv{:, 'std_kpv'} = nan; % Standard deviation of p50 kpv
T_classifier_kpv{:, 'vrb_kpv'} = nan; % Variability of p50 kpv
T_classifier_kpv{:, 'wmean_kpv'} = nan; % Mean of p75 - p25
for i = 1: height(T_classifier_kpv)
    kpv_1h_p25 = kpv(ar_datetime_hour==T_classifier_kpv{i, 'HOUR_START'}, 2);
    kpv_1h_p50 = kpv(ar_datetime_hour==T_classifier_kpv{i, 'HOUR_START'}, 3);
    kpv_1h_p75 = kpv(ar_datetime_hour==T_classifier_kpv{i, 'HOUR_START'}, 4);

    T_classifier_kpv{i, 'mean_kpv'} = mean(kpv_1h_p50, 'omitnan');
    T_classifier_kpv{i, 'std_kpv'} = std(kpv_1h_p50, 0, 'omitnan');
    T_classifier_kpv{i, 'vrb_kpv'} = sqrt(1/5.*sum(diff(kpv_1h_p50).^2));
    T_classifier_kpv{i, 'wmean_kpv'} = mean(kpv_1h_p75 - kpv_1h_p25, 'omitnan');
end
T_classifier_kpv.DATE = datetime(T_classifier_kpv.HOUR_START.Year, T_classifier_kpv.HOUR_START.Month, T_classifier_kpv.HOUR_START.Day, 'TimeZone', 'UTC');

fprintf('kNN\n');
cell_knnkpv_rtpd = cell(4, 1);
for c = 1:4
    T_knnkpv_rtpd = T_classifier_kpv((ar_datetime_hour_unique>=testhour_start)&(ar_datetime_hour_unique<=testhour_end), :); % Result container
    switch c
        case 1
            T_classifier_kpv{:, 'c'} = T_classifier_kpv{:, 'mean_kpv'};
            T_knnkpv_rtpd{:, 'c'}    = T_knnkpv_rtpd{:, 'mean_kpv'};
        case 2
            T_classifier_kpv{:, 'c'} = T_classifier_kpv{:, 'std_kpv'};
            T_knnkpv_rtpd{:, 'c'}    = T_knnkpv_rtpd{:, 'std_kpv'};
        case 3
            T_classifier_kpv{:, 'c'} = T_classifier_kpv{:, 'vrb_kpv'};
            T_knnkpv_rtpd{:, 'c'}    = T_knnkpv_rtpd{:, 'vrb_kpv'};
        case 4
            T_classifier_kpv{:, 'c'} = T_classifier_kpv{:, 'wmean_kpv'};
            T_knnkpv_rtpd{:, 'c'}    = T_knnkpv_rtpd{:, 'wmean_kpv'};
    end
    for k = karray
        T_knnkpv_rtpd{:, strcat('FRU_k', num2str(k))} = nan;
        T_knnkpv_rtpd{:, strcat('FRD_k', num2str(k))} = nan;

        for i = 1: size(T_knnkpv_rtpd, 1)
            this_date = datetime(T_knnkpv_rtpd.HOUR_START.Year(i), T_knnkpv_rtpd.HOUR_START.Month(i), T_knnkpv_rtpd.HOUR_START.Day(i), 'TimeZone', 'UTC');
            this_hour = T_knnkpv_rtpd.HOUR_START.Hour(i);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Can use the days AFTER this_date if results are not ideal
            T_sample = T_classifier_kpv((T_classifier_kpv.HOUR_START.Hour==this_hour), :);
%             T_sample = T_classifier_kpv((T_classifier_kpv.DATE<this_date)&(T_classifier_kpv.HOUR_START.Hour==this_hour), :);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if any(isnan(T_knnkpv_rtpd{i, 'c'}))
                T_sample_sorted = T_sample(T_sample.DATE>=this_date-days(k), :); % We use 30 previous days, i.e., baseline, if data is nan
                selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); % 30 the nearest days for baseline
            else
                T_sample.dist = sqrt(sum((T_sample{:, 'c'}-T_knnkpv_rtpd{i, 'c'}).^2, 2)); % Euclidean distance
                T_sample_sorted = sortrows(T_sample, 'dist');
                selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); 
            end

            sample_error_max = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_max'};
            sample_error_min = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_min'};
            [f,x] = ecdf(sample_error_max(:));
            T_knnkpv_rtpd{i, strcat('FRU_k', num2str(k))} = interp1(f, x, 0.975);
            [f,x] = ecdf(sample_error_min(:));
            T_knnkpv_rtpd{i, strcat('FRD_k', num2str(k))} = interp1(f, x, 0.025);
        end
        fprintf('c = %g, k = %g\n', c, k);
    end
    cell_knnkpv_rtpd{c} = T_knnkpv_rtpd;
end

%% Process results and visualize
f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});
T_actual_error = [T_baseline_rtpd(:, 'HOUR_START') T_errormax_rtpd T_errormin_rtpd];
ar_actual_error_fru = T_actual_error{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}};
ar_actual_error_frd = T_actual_error{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}};

% Baseline
for k = karray
    [oversupply_baseline_fru(karray==k), pshort_baseline_fru(karray==k)] = calculate_oversupply_and_pshort(ar_actual_error_fru, T_baseline_rtpd{:, strcat('FRU_k', num2str(k))});
    [oversupply_baseline_frd(karray==k), pshort_baseline_frd(karray==k)] = calculate_oversupply_and_pshort(-ar_actual_error_frd, -T_baseline_rtpd{:, strcat('FRD_k', num2str(k))});
end

% kNN
for c = 1:4
    T_knnkpv_rtpd = cell_knnkpv_rtpd{c};
    for k = karray
        [oversupply_knnkpv_fru(karray==k, c), pshort_knnkpv_fru(karray==k, c)] = calculate_oversupply_and_pshort(ar_actual_error_fru, T_knnkpv_rtpd{:, strcat('FRU_k', num2str(k))});
        [oversupply_knnkpv_frd(karray==k, c), pshort_knnkpv_frd(karray==k, c)] = calculate_oversupply_and_pshort(-ar_actual_error_frd, -T_knnkpv_rtpd{:, strcat('FRD_k', num2str(k))});
    end
end

%% Produce figures
figure(); 
h = plot(karray, oversupply_knnkpv_frd./1E3, karray, oversupply_baseline_frd./1E3, 'k');
set(h(1), 'Color', 'r');
set(h(2), 'Color', 'g');
set(h(3), 'Color', 'm');
set(h(4), 'Color', 'b');
set(h(5), 'Color', 'k');
title('Oversupply (FRD)');
xlabel('K');
ylabel('Oversupply (GWh)');
legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');

figure(); 
h = plot(karray, pshort_knnkpv_frd, karray, pshort_baseline_frd, 'k'); 
set(h(1), 'Color', 'r');
set(h(2), 'Color', 'g');
set(h(3), 'Color', 'm');
set(h(4), 'Color', 'b');
set(h(5), 'Color', 'k');
title('Pshort (FRD)');
xlabel('K');
ylabel('Probability of shortage');
legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');

figure(); 
h = plot(karray, oversupply_knnkpv_fru./1E3, karray, oversupply_baseline_fru./1E3, 'k'); 
set(h(1), 'Color', 'r');
set(h(2), 'Color', 'g');
set(h(3), 'Color', 'm');
set(h(4), 'Color', 'b');
set(h(5), 'Color', 'k');
title('Oversupply (FRU)');
xlabel('K');
ylabel('Oversupply (GWh)');
legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');

figure(); 
h = plot(karray, pshort_knnkpv_fru, karray, pshort_baseline_fru, 'k'); 
set(h(1), 'Color', 'r');
set(h(2), 'Color', 'g');
set(h(3), 'Color', 'm');
set(h(4), 'Color', 'b');
set(h(5), 'Color', 'k');
title('Pshort (FRU)');
xlabel('K');
ylabel('Probability of shortage');
legend('\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');

% Pareto frontier
figure(); 
h(5)= scatter(pshort_baseline_fru, oversupply_baseline_fru./1E3, 60, 'ko');
hold on;
for c = 1: 4
    h(c) = scatter(pshort_knnkpv_fru(:, c), oversupply_knnkpv_fru(:, c)./1E3, 30, '^');
end
set(h(1), 'MarkerFaceColor', 'r');
set(h(2), 'MarkerFaceColor', 'g');
set(h(3), 'MarkerFaceColor', 'm');
set(h(4), 'MarkerFaceColor', 'b');
set(h(5), 'MarkerFaceColor', 'k');
xline(pshort_baseline_fru(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(oversupply_baseline_fru(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
title('FRU');
xlabel('Probability of shortage');
ylabel('Oversupply (GWh)');
legend(h, '\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
box on;

figure(); 
h(5)= scatter(pshort_baseline_frd, oversupply_baseline_frd./1E3, 60, 'ko');
hold on;
for c = 1: 4
    h(c) = scatter(pshort_knnkpv_frd(:, c), oversupply_knnkpv_frd(:, c)./1E3, 30, '^');
end
set(h(1), 'MarkerFaceColor', 'r');
set(h(2), 'MarkerFaceColor', 'g');
set(h(3), 'MarkerFaceColor', 'm');
set(h(4), 'MarkerFaceColor', 'b');
set(h(5), 'MarkerFaceColor', 'k');
xline(pshort_baseline_frd(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(oversupply_baseline_frd(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
title('FRD');
xlabel('Probability of shortage');
ylabel('Oversupply (GWh)');
legend(h, '\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
box on;

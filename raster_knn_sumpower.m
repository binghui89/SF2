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
    sigma_psolarcs = sum(p_solar_cs*diag(T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'TotalCapacity'}), 2);
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

% Normalize all classifiers
T_classifier_kpv_norm = T_classifier_kpv(:, {'HOUR_START', 'mean_kpv', 'std_kpv', 'vrb_kpv', 'DATE'});
T_classifier_kpv_norm{:, 'mean_kpv'} = T_classifier_kpv{:, 'mean_kpv'}./max(T_classifier_kpv{:, 'mean_kpv'});
T_classifier_kpv_norm{:, 'std_kpv'} = T_classifier_kpv{:, 'std_kpv'}./max(T_classifier_kpv{:, 'std_kpv'});
T_classifier_kpv_norm{:, 'vrb_kpv'} = T_classifier_kpv{:, 'vrb_kpv'}./max(T_classifier_kpv{:, 'vrb_kpv'});

%% kNN 1-dim
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
                T_sample.dist = sqrt(sum((T_sample{:, 'c'}-T_knnkpv_rtpd{i, 'c'}).^2, 2)); % Manhattan distance
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

%% knn 3-dim
fprintf('kNN (n-dim)\n');

w = []; % Weight, increment of 0.1 at each dimension
for i = 0: 10
    for j = 0:10-i
        k=10-i-j;
        w = [w, [i/10; j/10; k/10]];
    end
end

cell_knnkpv_rtpd_ndim = cell(size(w, 2), 1);
name_classifiers = {'mean_kpv', 'std_kpv', 'vrb_kpv'};
for j = 1: size(w, 2)
    w_this = w(:, j);
    lw_this = true(size(w_this));
    lw_this(w_this==0) = false;
    T_knnkpv_rtpd_ndim = T_classifier_kpv_norm((ar_datetime_hour_unique>=testhour_start)&(ar_datetime_hour_unique<=testhour_end), :); % Result container
    for k = karray
        T_knnkpv_rtpd_ndim{:, strcat('FRU_k', num2str(k))} = nan;
        T_knnkpv_rtpd_ndim{:, strcat('FRD_k', num2str(k))} = nan;

        for i = 1: size(T_knnkpv_rtpd_ndim, 1)
            this_date = datetime(T_knnkpv_rtpd_ndim.HOUR_START.Year(i), T_knnkpv_rtpd_ndim.HOUR_START.Month(i), T_knnkpv_rtpd_ndim.HOUR_START.Day(i), 'TimeZone', 'UTC');
            this_hour = T_knnkpv_rtpd_ndim.HOUR_START.Hour(i);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Can use the days AFTER this_date if results are not ideal
            T_sample = T_classifier_kpv_norm((T_classifier_kpv_norm.HOUR_START.Hour==this_hour), :);
    %             T_sample = T_classifier_kpv((T_classifier_kpv.DATE<this_date)&(T_classifier_kpv.HOUR_START.Hour==this_hour), :);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if any(isnan(T_knnkpv_rtpd_ndim{i, {'mean_kpv', 'std_kpv', 'vrb_kpv'}}), 2)
            if any(isnan(T_knnkpv_rtpd_ndim{i, name_classifiers(lw_this)}), 2)
                T_sample_sorted = T_sample(T_sample.DATE>=this_date-days(k), :); % We use 30 previous days, i.e., baseline, if data is nan
                selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); % 30 the nearest days for baseline
            else
                T_sample.dist = abs(T_sample{:, name_classifiers(lw_this)}-T_knnkpv_rtpd_ndim{i, name_classifiers(lw_this)})*w(lw_this, j); % Manhattan distance
%                 T_sample.dist = abs(T_sample{:, {'mean_kpv', 'std_kpv', 'vrb_kpv'}}-T_knnkpv_rtpd_ndim{i, {'mean_kpv', 'std_kpv', 'vrb_kpv'}})*w(:, j); % Manhattan distance
                T_sample_sorted = sortrows(T_sample, 'dist');
                selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); 
            end

            sample_error_max = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_max'};
            sample_error_min = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_min'};
            [f,x] = ecdf(sample_error_max(:));
            T_knnkpv_rtpd_ndim{i, strcat('FRU_k', num2str(k))} = interp1(f, x, 0.975);
            [f,x] = ecdf(sample_error_min(:));
            T_knnkpv_rtpd_ndim{i, strcat('FRD_k', num2str(k))} = interp1(f, x, 0.025);
        end
        fprintf('j = %g, k = %g\n', j, k);
    end
    cell_knnkpv_rtpd_ndim{j, 1} = T_knnkpv_rtpd_ndim;
end

%% Process results and visualize
f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});
T_actual_error = [T_baseline_rtpd(:, 'HOUR_START') T_errormax_rtpd T_errormin_rtpd];
ar_actual_error_fru = T_actual_error{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}};
ar_actual_error_frd = T_actual_error{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}};

% Calculate sum of absolute errors: unfloor
ar_actual_fru_unfloor = max(ar_actual_error_fru, [], 2);
ar_actual_frd_unfloor = min(ar_actual_error_frd, [], 2);
ar_sae_fru_unfloor = nan(size(w, 2), numel(karray));
ar_sae_frd_unfloor = nan(size(w, 2), numel(karray));

for j = 1: size(w, 2)
    T_knnkpv_rtpd_ndim = cell_knnkpv_rtpd_ndim{j, 1};
    for k = karray
        ar_sae_fru_unfloor(j, karray==k) = sum(abs(ar_actual_fru_unfloor - T_knnkpv_rtpd_ndim{:, strcat('FRU_k', num2str(k))}));
        ar_sae_frd_unfloor(j, karray==k) = sum(abs(ar_actual_frd_unfloor - T_knnkpv_rtpd_ndim{:, strcat('FRD_k', num2str(k))}));
    end
end

% Calculate sum of absolute errors: floored
ar_actual_fru_floor = max(0, max(ar_actual_error_fru, [], 2));
ar_actual_frd_floor = min(0, min(ar_actual_error_frd, [], 2));
ar_sae_fru_floor = nan(size(w, 2), numel(karray));
ar_sae_frd_floor = nan(size(w, 2), numel(karray));

for j = 1: size(w, 2)
    T_knnkpv_rtpd_ndim = cell_knnkpv_rtpd_ndim{j, 1};
    for k = karray
        ar_sae_fru_floor(j, karray==k) = sum(abs(ar_actual_fru_floor - max(0, T_knnkpv_rtpd_ndim{:, strcat('FRU_k', num2str(k))})));
        ar_sae_frd_floor(j, karray==k) = sum(abs(ar_actual_frd_floor - min(0, T_knnkpv_rtpd_ndim{:, strcat('FRD_k', num2str(k))})));
    end
end


%% Process results and visualize: Calculate oversupply and p_shortage
% Baseline
for k = karray
    [oversupply_baseline_fru(karray==k), pshort_baseline_fru(karray==k)] = calculate_oversupply_and_pshort(ar_actual_error_fru, T_baseline_rtpd{:, strcat('FRU_k', num2str(k))});
    [oversupply_baseline_frd(karray==k), pshort_baseline_frd(karray==k)] = calculate_oversupply_and_pshort(-ar_actual_error_frd, -T_baseline_rtpd{:, strcat('FRD_k', num2str(k))});
end

% kNN, 1-dim
for c = 1:4
    T_knnkpv_rtpd = cell_knnkpv_rtpd{c};
    for k = karray
        [oversupply_knnkpv_fru(karray==k, c), pshort_knnkpv_fru(karray==k, c)] = calculate_oversupply_and_pshort(ar_actual_error_fru, T_knnkpv_rtpd{:, strcat('FRU_k', num2str(k))});
        [oversupply_knnkpv_frd(karray==k, c), pshort_knnkpv_frd(karray==k, c)] = calculate_oversupply_and_pshort(-ar_actual_error_frd, -T_knnkpv_rtpd{:, strcat('FRD_k', num2str(k))});
    end
end

% kNN, 3-dim
% for k = karray
%     [oversupply_knnkpvndim_fru(karray==k), pshort_knnkpvndim_fru(karray==k)] = calculate_oversupply_and_pshort(ar_actual_error_fru, T_knnkpv_rtpd_ndim{:, strcat('FRU_k', num2str(k))});
%     [oversupply_knnkpvndim_frd(karray==k), pshort_knnkpvndim_frd(karray==k)] = calculate_oversupply_and_pshort(-ar_actual_error_frd, -T_knnkpv_rtpd_ndim{:, strcat('FRD_k', num2str(k))});
% end

% kNN, 3-dim
% Calculate kNN exhaustive search oversupply and pshort 
for j = 1: size(w, 2)
    T_knnkpv_rtpd_ndim = cell_knnkpv_rtpd_ndim{j, 1};
    for k = karray
        [oversupply_knnkpvndimopt_fru(j, karray==k), pshort_knnkpvndimopt_fru(j, karray==k)] = calculate_oversupply_and_pshort(ar_actual_error_fru, T_knnkpv_rtpd_ndim{:, strcat('FRU_k', num2str(k))});
        [oversupply_knnkpvndimopt_frd(j, karray==k), pshort_knnkpvndimopt_frd(j, karray==k)] = calculate_oversupply_and_pshort(-ar_actual_error_frd, -T_knnkpv_rtpd_ndim{:, strcat('FRD_k', num2str(k))});
    end
end


%% Produce figures: Trends as a function of K
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

%% Produce figures: Pareto frontier, 1-dim
figure(); 
h(5)= scatter(pshort_baseline_fru, oversupply_baseline_fru./1E3, 60, 'ko');
hold on;
for c = 1: 4
    h(c) = scatter(pshort_knnkpv_fru(:, c), oversupply_knnkpv_fru(:, c)./1E3, 30, '^');
end
% h(6) = scatter(pshort_knnkpvndim_fru(:), oversupply_knnkpvndim_fru(:)./1E3, 30, '^');
set(h(1), 'MarkerFaceColor', 'r');
set(h(2), 'MarkerFaceColor', 'g');
set(h(3), 'MarkerFaceColor', 'm');
set(h(4), 'MarkerFaceColor', 'b');
set(h(5), 'MarkerFaceColor', 'k');
% set(h(6), 'MarkerFaceColor', 'y');
xline(pshort_baseline_fru(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(oversupply_baseline_fru(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
title('FRU');
xlabel('Probability of shortage');
ylabel('Oversupply (GWh)');
legend(h, '\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the same figure with lines, not scatter plots
% figure();
% h(5)= line(pshort_baseline_fru, oversupply_baseline_fru./1E3);
% hold on;
% for c = 1: 4
%     h(c) = line(pshort_knnkpv_fru(:, c), oversupply_knnkpv_fru(:, c)./1E3);
% end
% set(h(1), 'Color', 'r', 'MarkerFaceColor', 'r', 'Marker', '^', 'MarkerSize', 6);
% set(h(2), 'Color', 'g', 'MarkerFaceColor', 'g', 'Marker', '^', 'MarkerSize', 6);
% set(h(3), 'Color', 'm', 'MarkerFaceColor', 'm', 'Marker', '^', 'MarkerSize', 6);
% set(h(4), 'Color', 'b', 'MarkerFaceColor', 'b', 'Marker', '^', 'MarkerSize', 6);
% set(h(5), 'Color', 'k', 'MarkerFaceColor', 'k', 'Marker', 'o', 'MarkerSize', 6);
% xline(pshort_baseline_fru(karray==30), 'Color', 'r', 'LineStyle', '--');
% yline(oversupply_baseline_fru(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
% title('FRU');
% xlabel('Probability of shortage');
% ylabel('Oversupply (GWh)');
% legend(h, '\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
% box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(); 
h(5)= scatter(pshort_baseline_frd, oversupply_baseline_frd./1E3, 60, 'ko');
hold on;
for c = 1: 4
    h(c) = scatter(pshort_knnkpv_frd(:, c), oversupply_knnkpv_frd(:, c)./1E3, 30, '^');
end
% h(6) = scatter(pshort_knnkpvndim_frd(:), oversupply_knnkpvndim_frd(:)./1E3, 30, '^');
set(h(1), 'MarkerFaceColor', 'r');
set(h(2), 'MarkerFaceColor', 'g');
set(h(3), 'MarkerFaceColor', 'm');
set(h(4), 'MarkerFaceColor', 'b');
set(h(5), 'MarkerFaceColor', 'k');
% set(h(6), 'MarkerFaceColor', 'y');
xline(pshort_baseline_frd(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(oversupply_baseline_frd(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
title('FRD');
xlabel('Probability of shortage');
ylabel('Oversupply (GWh)');
legend(h, '\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the same figure with lines, not scatter plots
% figure();
% h(5)= line(pshort_baseline_frd, oversupply_baseline_frd./1E3);
% hold on;
% for c = 1: 4
%     h(c) = line(pshort_knnkpv_frd(:, c), oversupply_knnkpv_frd(:, c)./1E3);
% end
% set(h(1), 'Color', 'r', 'MarkerFaceColor', 'r', 'Marker', '^', 'MarkerSize', 6);
% set(h(2), 'Color', 'g', 'MarkerFaceColor', 'g', 'Marker', '^', 'MarkerSize', 6);
% set(h(3), 'Color', 'm', 'MarkerFaceColor', 'm', 'Marker', '^', 'MarkerSize', 6);
% set(h(4), 'Color', 'b', 'MarkerFaceColor', 'b', 'Marker', '^', 'MarkerSize', 6);
% set(h(5), 'Color', 'k', 'MarkerFaceColor', 'k', 'Marker', 'o', 'MarkerSize', 6);
% xline(pshort_baseline_frd(karray==30), 'Color', 'r', 'LineStyle', '--');
% yline(oversupply_baseline_frd(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
% title('FRD');
% xlabel('Probability of shortage');
% ylabel('Oversupply (GWh)');
% legend(h, '\mu(k_{pv})', '\sigma(k_{pv})', 'v(k_{pv})', 'w(k_{pv})' ,'Baseline');
% box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Color by dimension of classifiers
figure();
hold on;
h = [];
for ndim = 1:3
    pshort_frd_tmp = pshort_knnkpvndimopt_frd(sum(w==0, 1)==3-ndim, :);
    oversupply_frd_tmp = oversupply_knnkpvndimopt_frd(sum(w==0, 1)==3-ndim, :);
    h(ndim) = scatter(pshort_frd_tmp(:), oversupply_frd_tmp(:)./1e3, 'o');
end
h(numel(h)+1)= scatter(pshort_baseline_frd, oversupply_baseline_frd./1E3, 60, 'ko');
set(h(1), 'MarkerFaceColor', 'r', 'Marker', '^');
set(h(2), 'MarkerEdgeColor', 'b', 'Marker', '+');
set(h(3), 'MarkerFaceColor', 'm', 'Marker', '.');
set(h(4), 'MarkerFaceColor', 'k');
xline(pshort_baseline_frd(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(oversupply_baseline_frd(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
legend(h, '1-dim', '2-dim', '3-dim', 'Baseline');
title('FRD');

figure();
hold on;
h = [];
for ndim = 1:3
    pshort_fru_tmp = pshort_knnkpvndimopt_fru(sum(w==0, 1)==3-ndim, :);
    oversupply_fru_tmp = oversupply_knnkpvndimopt_fru(sum(w==0, 1)==3-ndim, :);
    h(ndim) = scatter(pshort_fru_tmp(:), oversupply_fru_tmp(:)./1e3);
end
h(numel(h)+1)= scatter(pshort_baseline_fru, oversupply_baseline_fru./1E3, 60, 'ko');
xline(pshort_baseline_fru(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(oversupply_baseline_fru(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
legend(h, '1-dim', '2-dim', '3-dim', 'Baseline');
set(h(1), 'MarkerFaceColor', 'r', 'Marker', '^');
set(h(2), 'MarkerEdgeColor', 'b', 'Marker', '+');
set(h(3), 'MarkerFaceColor', 'm', 'Marker', '.');
set(h(4), 'MarkerFaceColor', 'k');
title('FRU');

trainhour_end = datetime(2020, 3, 25, 2, 0, 0, 'TimeZone', 'UTC');
trainhour_start = datetime(2020, 3, 4, 9, 0, 0, 'TimeZone', 'UTC');

T_actual_error{:, 'istrain'} = 0;
T_actual_error{(T_actual_error.HOUR_START<=trainhour_end)&(T_actual_error.HOUR_START>=trainhour_start), 'istrain'} = 1;

for j = 1: size(w, 2)
    T_knnkpv_rtpd_ndim = cell_knnkpv_rtpd_ndim{j, 1};
    T_knnkpv_rtpd_ndim{:, 'istrain'} = 0;
    T_knnkpv_rtpd_ndim{(T_knnkpv_rtpd_ndim.HOUR_START<=trainhour_end)&(T_knnkpv_rtpd_ndim.HOUR_START>=trainhour_start), 'istrain'} = 1;
    cell_knnkpv_rtpd_ndim{j, 1} = T_knnkpv_rtpd_ndim;
end

all_h_local = 5:1:18; % Hourly knn: UTC 1300 to 0200 (+1), or PST 5am to 6pm
all_h_utc = all_h_local + 8;
all_h_utc(all_h_utc >=24) = all_h_utc(all_h_utc >=24) - 24;


%% Multi-objective optimization

% Step I: Calculate evaluating metrics of the training dataset by hour
o_knnkpvnd_byh_fru_train = nan(size(w, 2), numel(karray), numel(all_h_local));
o_knnkpvnd_byh_frd_train = nan(size(w, 2), numel(karray), numel(all_h_local));
pr_knnkpvnd_byh_fru_train = nan(size(w, 2), numel(karray), numel(all_h_local));
pr_knnkpvnd_byh_frd_train = nan(size(w, 2), numel(karray), numel(all_h_local));

for ih = 1:length(all_h_local)
    h_utc = all_h_utc(ih);
    h_local = all_h_local(ih);
    for j = 1: size(w, 2)
        T_knnkpv_rtpd_ndim = cell_knnkpv_rtpd_ndim{j, 1};
        for k = karray
            [o_knnkpvnd_byh_fru_train(j, karray==k, all_h_local==h_local), pr_knnkpvnd_byh_fru_train(j, karray==k, all_h_local==h_local)] = calculate_oversupply_and_pshort(...
                ar_actual_error_fru((T_actual_error.HOUR_START.Hour==h_utc)&(T_actual_error.istrain==1)), ...
                T_knnkpv_rtpd_ndim{(T_knnkpv_rtpd_ndim.HOUR_START.Hour==h_utc)&(T_knnkpv_rtpd_ndim.istrain==1), strcat('FRU_k', num2str(k))} ...
                );
            [o_knnkpvnd_byh_frd_train(j, karray==k, all_h_local==h_local), pr_knnkpvnd_byh_frd_train(j, karray==k, all_h_local==h_local)] = calculate_oversupply_and_pshort(...
                -ar_actual_error_frd((T_actual_error.HOUR_START.Hour==h_utc)&(T_actual_error.istrain==1)), ...
                -T_knnkpv_rtpd_ndim{(T_knnkpv_rtpd_ndim.HOUR_START.Hour==h_utc)&(T_knnkpv_rtpd_ndim.istrain==1), strcat('FRD_k', num2str(k))} ...
                );
        end
    end
end

% Calculate the evaluating metrics of the test dataset by hour
o_knnkpvnd_byh_fru_test = nan(size(w, 2), numel(karray), numel(all_h_local));
o_knnkpvnd_byh_frd_test = nan(size(w, 2), numel(karray), numel(all_h_local));
pr_knnkpvnd_byh_fru_test = nan(size(w, 2), numel(karray), numel(all_h_local));
pr_knnkpvnd_byh_frd_test = nan(size(w, 2), numel(karray), numel(all_h_local));

% Calculate
for ih = 1:length(all_h_local)
    h_utc = all_h_utc(ih);
    h_local = all_h_local(ih);
    for j = 1: size(w, 2)
        T_knnkpv_rtpd_ndim = cell_knnkpv_rtpd_ndim{j, 1};
        for k = karray
            [o_knnkpvnd_byh_fru_test(j, karray==k, all_h_local==h_local), pr_knnkpvnd_byh_fru_test(j, karray==k, all_h_local==h_local)] = calculate_oversupply_and_pshort(...
                ar_actual_error_fru((T_actual_error.HOUR_START.Hour==h_utc)&(T_actual_error.istrain==0)), ...
                T_knnkpv_rtpd_ndim{(T_knnkpv_rtpd_ndim.HOUR_START.Hour==h_utc)&(T_knnkpv_rtpd_ndim.istrain==0), strcat('FRU_k', num2str(k))} ...
                );
            [o_knnkpvnd_byh_frd_test(j, karray==k, all_h_local==h_local), pr_knnkpvnd_byh_frd_test(j, karray==k, all_h_local==h_local)] = calculate_oversupply_and_pshort(...
                -ar_actual_error_frd((T_actual_error.HOUR_START.Hour==h_utc)&(T_actual_error.istrain==0)), ...
                -T_knnkpv_rtpd_ndim{(T_knnkpv_rtpd_ndim.HOUR_START.Hour==h_utc)&(T_knnkpv_rtpd_ndim.istrain==0), strcat('FRD_k', num2str(k))} ...
                );
        end
    end
end

%% Start Multi-objective optimization
alpha = [0.01:0.02:0.50];

% THe real optimal
T_knnkpv_rtpd_ndimopt = T_actual_error(T_actual_error.istrain==0, {'HOUR_START'});
T_knnkpv_rtpd_ndimopt{:, 'FRU'} = 0;
T_knnkpv_rtpd_ndimopt{:, 'FRD'} = 0;

for ialpha = 1: length(alpha)
    % The kNN optimal
    T_knnkpv_rtpd_ndimopt_ = T_actual_error(T_actual_error.istrain==0, {'HOUR_START'});
    T_knnkpv_rtpd_ndimopt_{:, 'FRU'} = 0;
    T_knnkpv_rtpd_ndimopt_{:, 'FRD'} = 0;

    % Reliability first
    % FRD
    for ih = 1:length(all_h_local)

    %     ih = 13;
        % Minimum using training data
        [min_pr_, ~] = min(pr_knnkpvnd_byh_frd_train(:,:, ih), [], 'all','linear');
        imin_all_ = find((pr_knnkpvnd_byh_frd_train(:,:, ih)-min_pr_)<=alpha(ialpha));
        tmp = o_knnkpvnd_byh_frd_train(:, :, ih);
        [min_o_c_minpr_, ~] = min(tmp(imin_all_));
        imin_o_c_minpr_ = find(abs((o_knnkpvnd_byh_frd_train(:,:, ih)-min_o_c_minpr_))<=0.0001);
        if isempty(imin_o_c_minpr_)
            error('No optimal value found!');
        else
            [r_frdopt_(ih), c_frdopt_(ih)] = ind2sub(size(tmp), imin_o_c_minpr_(1));
        end

        % Real minimum
        [min_pr, ~] = min(pr_knnkpvnd_byh_frd_test(:,:, ih), [], 'all','linear');
        imin_all = find((pr_knnkpvnd_byh_frd_test(:,:, ih)-min_pr)<=alpha(ialpha));
        tmp = o_knnkpvnd_byh_frd_test(:, :, ih);
        [min_o_c_minpr, ~] = min(tmp(imin_all));
        imin_o_c_minpr = find(abs((o_knnkpvnd_byh_frd_test(:,:, ih)-min_o_c_minpr))<=0.0001);
        if isempty(imin_o_c_minpr)
            error('No optimal value found!');
        else
            [r_frdopt(ih), c_frdopt(ih)] = ind2sub(size(tmp), imin_o_c_minpr(1));
        end

        error_frd_pr(ih) = pr_knnkpvnd_byh_frd_test(r_frdopt_(ih), c_frdopt_(ih), ih) - pr_knnkpvnd_byh_frd_test(r_frdopt(ih), c_frdopt(ih), ih);
        error_frd_o(ih)  = o_knnkpvnd_byh_frd_test(r_frdopt_(ih), c_frdopt_(ih), ih) - o_knnkpvnd_byh_frd_test(r_frdopt(ih), c_frdopt(ih), ih);

        % Get the optimal FRD and FRU
        tmp = cell_knnkpv_rtpd_ndim{r_frdopt_(ih)};
        T_knnkpv_rtpd_ndimopt_{T_knnkpv_rtpd_ndimopt_.HOUR_START.Hour==all_h_utc(ih), 'FRD'} = tmp{(tmp.HOUR_START.Hour==all_h_utc(ih))&(tmp.istrain==0), strcat( 'FRD_k', num2str(karray( c_frdopt_(ih) )) )};

        tmp = cell_knnkpv_rtpd_ndim{r_frdopt(ih)};
        T_knnkpv_rtpd_ndimopt{T_knnkpv_rtpd_ndimopt.HOUR_START.Hour==all_h_utc(ih), 'FRD'} = tmp{(tmp.HOUR_START.Hour==all_h_utc(ih))&(tmp.istrain==0), strcat( 'FRD_k', num2str(karray( c_frdopt(ih) )) )};    
    end

    % FRU
    for ih = 1:length(all_h_local)
        % Minimum using training data
        [min_pr_, ~] = min(pr_knnkpvnd_byh_fru_train(:,:, ih), [], 'all','linear');
        imin_all_ = find((pr_knnkpvnd_byh_fru_train(:,:, ih)-min_pr_)<=0.0001);
        tmp = o_knnkpvnd_byh_fru_train(:, :, ih);
        [min_o_c_minpr_, ~] = min(tmp(imin_all_));
        imin_o_c_minpr_ = find(abs((o_knnkpvnd_byh_fru_train(:,:, ih)-min_o_c_minpr_))<=0.0001);
        if isempty(imin_o_c_minpr_)
            error('No optimal value found!');
        else
            [r_fruopt_(ih), c_fruopt_(ih)] = ind2sub(size(tmp), imin_o_c_minpr_(1));
        end

        % Real minimum
        [min_pr, ~] = min(pr_knnkpvnd_byh_fru_test(:,:, ih), [], 'all','linear');
        imin_all = find((pr_knnkpvnd_byh_fru_test(:,:, ih)-min_pr)<=0.0001);
        tmp = o_knnkpvnd_byh_fru_test(:, :, ih);
        [min_o_c_minpr, ~] = min(tmp(imin_all));
        imin_o_c_minpr = find(abs((o_knnkpvnd_byh_fru_test(:,:, ih)-min_o_c_minpr))<=0.0001);
        if isempty(imin_o_c_minpr)
            error('No optimal value found!');
        else
            [r_fruopt(ih), c_fruopt(ih)] = ind2sub(size(tmp), imin_o_c_minpr(1));
        end

        error_fru_pr(ih) = pr_knnkpvnd_byh_fru_test(r_fruopt_(ih), c_fruopt_(ih), ih) - pr_knnkpvnd_byh_fru_test(r_fruopt(ih), c_fruopt(ih), ih);
        error_fru_o(ih)  = o_knnkpvnd_byh_fru_test(r_fruopt_(ih), c_fruopt_(ih), ih) - o_knnkpvnd_byh_fru_test(r_fruopt(ih), c_fruopt(ih), ih);

        % Get the optimal FRD and FRU
        tmp = cell_knnkpv_rtpd_ndim{r_fruopt_(ih)};
        T_knnkpv_rtpd_ndimopt_{T_knnkpv_rtpd_ndimopt_.HOUR_START.Hour==all_h_utc(ih), 'FRU'} = tmp{(tmp.HOUR_START.Hour==all_h_utc(ih))&(tmp.istrain==0), strcat( 'FRU_k', num2str(karray( c_fruopt_(ih) )) )};

        tmp = cell_knnkpv_rtpd_ndim{r_fruopt(ih)};
        T_knnkpv_rtpd_ndimopt{T_knnkpv_rtpd_ndimopt.HOUR_START.Hour==all_h_utc(ih), 'FRU'} = tmp{(tmp.HOUR_START.Hour==all_h_utc(ih))&(tmp.istrain==0), strcat( 'FRU_k', num2str(karray( c_fruopt(ih) )) )};
    end

    % Calculate oversupply and pshortage of the estimated optimal kNN requirements
    [o_knnkpvnd_opt_fru_(ialpha), pr_knnkpvnd_opt_fru_(ialpha)] = calculate_oversupply_and_pshort(ar_actual_error_fru((T_actual_error.istrain==0)), T_knnkpv_rtpd_ndimopt_{:, 'FRU'});
    [o_knnkpvnd_opt_frd_(ialpha), pr_knnkpvnd_opt_frd_(ialpha)] = calculate_oversupply_and_pshort(-ar_actual_error_frd((T_actual_error.istrain==0)), -T_knnkpv_rtpd_ndimopt_{:, 'FRD'});

    % Calculate oversupply and pshortage of the real optimal kNN requirements
    [o_knnkpvnd_opt_fru(ialpha), pr_knnkpvnd_opt_fru(ialpha)] = calculate_oversupply_and_pshort(ar_actual_error_fru((T_actual_error.istrain==0)), T_knnkpv_rtpd_ndimopt{:, 'FRU'});
    [o_knnkpvnd_opt_frd(ialpha), pr_knnkpvnd_opt_frd(ialpha)] = calculate_oversupply_and_pshort(-ar_actual_error_frd((T_actual_error.istrain==0)), -T_knnkpv_rtpd_ndimopt{:, 'FRD'});
end

% figure();line(pr_knnkpvnd_opt_frd_, o_knnkpvnd_opt_frd_);
% hold on;
% line(pr_knnkpvnd_opt_frd, o_knnkpvnd_opt_frd);

%% Save 
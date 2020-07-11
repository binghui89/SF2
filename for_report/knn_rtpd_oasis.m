dt_rtd = 5; % min
dt_rtpd = 15; % min

%% Load CAISO data RTD
T_rtd = readtable('C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\tmp_excels\df_rtd.csv', 'Delimiter',',');
T_rtd.TIME = datetime(T_rtd.Var1, 'InputFormat', 'yyyy-MM-dd'' ''HH:mm:ssXXX', 'TimeZone', 'UTC'); % This is the time of the end of an interval
T_rtd.TIME_START = T_rtd.TIME - duration(0, dt_rtd, 0); % This is the time of the start of an interval
T_rtd.HOUR_START = datetime(T_rtd.TIME_START.Year, T_rtd.TIME_START.Month, T_rtd.TIME_START.Day, T_rtd.TIME_START.Hour, 0, 0, 'TimeZone', 'UTC');
% T_rtd.local_time = datetime(T_rtd.TIME, 'TimeZone', 'America/Los_Angeles');
T_rtd.local_time = datetime(T_rtd.TIME, 'TimeZone', '-08:00'); % We use PST as local time

T_rtd.FORECAST_ERROR_Brtd_Artd_solar = (-T_rtd.Solar_NP15_B_RTD-T_rtd.Solar_ZP26_B_RTD-T_rtd.Solar_SP15_B_RTD) - (-T_rtd.Solar_NP15_A_RTD-T_rtd.Solar_ZP26_A_RTD-T_rtd.Solar_SP15_A_RTD); % Pure solar forecasting errors

% Demonstrate FRP requirements and binding interval forecasts, RTD
% ax1 = subplot(2, 1, 1);
% plot(T_rtd.TIME, T_rtd.UP_RTD, '-b', T_rtd.TIME, -1.*T_rtd.DOWN_RTD, '-b', T_rtd.TIME, T_rtd.FORECAST_ERROR_Brtd_Artd, '-r');
% 
% ax2 = subplot(2, 1, 2);
% plot(T_rtd.TIME, sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2), 'b');
% linkaxes([ax1(1), ax2],'x');

%% Load CAISO data RTPD
T_rtpd = readtable('C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\tmp_excels\df_rtpd.csv', 'Delimiter',',');
T_rtpd.TIME = datetime(T_rtpd.Var1, 'InputFormat', 'yyyy-MM-dd'' ''HH:mm:ssXXX', 'TimeZone', 'UTC');
T_rtpd.TIME_START = T_rtpd.TIME - duration(0, dt_rtpd, 0); % This is the time of the start of an interval
T_rtpd.HOUR_START = datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, T_rtpd.TIME_START.Hour, 0, 0, 'TimeZone', 'UTC');
T_rtpd.HOUR_START_local = datetime(T_rtpd.HOUR_START, 'TimeZone', 'America/Los_Angeles');
% T_rtpd.local_time = datetime(T_rtpd.TIME, 'TimeZone', 'America/Los_Angeles');

T_rtpd.TIME_local = datetime(T_rtpd.TIME, 'TimeZone', 'America/Los_Angeles');
T_rtpd.TIME_START_local = datetime(T_rtpd.TIME_START, 'TimeZone', 'America/Los_Angeles');
T_rtpd.DATE_START = datetime(T_rtpd.HOUR_START.Year, T_rtpd.HOUR_START.Month, T_rtpd.HOUR_START.Day, 'TimeZone', 'UTC');
T_rtpd.DATE_START_local = datetime(T_rtpd.HOUR_START_local.Year, T_rtpd.HOUR_START_local.Month, T_rtpd.HOUR_START_local.Day, 'TimeZone', 'America/Los_Angeles');

% calculate pure solar forecast error
tmp_rtp_solar  = sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2);
tmp_rtpd_solar = sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2);
tmp_error_solar = -reshape(tmp_rtp_solar, 3, size(tmp_rtp_solar, 1)/3)' - (-repmat(tmp_rtpd_solar, 1, 3));

% Calculate net load forecast error
tmp_rtd_nl_b  = T_rtd.LOAD_B_RTD - sum(T_rtd{:, {'Wind_NP15_B_RTD', 'Wind_SP15_B_RTD', 'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2);
tmp_rtpd_nl_a = T_rtpd.LOAD_B_RTPD - sum(T_rtpd{:, {'Wind_NP15_RTPD', 'Wind_SP15_RTPD', 'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2); % Note we use binding forecast as advisory forecast since CAISO does not publish advisory forecast
% tmp_rtd_nl_b  = - sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2);
% tmp_rtpd_nl_a = - sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2); % Note we use binding forecast as advisory forecast since CAISO does not publish advisory forecast

tmp_error_nl = reshape(tmp_rtd_nl_b, 3, size(tmp_rtd_nl_b, 1)/3)' - (repmat(tmp_rtpd_nl_a, 1, 3));

T_rtpd.error_max = max(tmp_error_nl, [], 2);
T_rtpd.error_min = min(tmp_error_nl, [], 2);


% Demonstrate FRP requirements and binding interval forecasts, RTPD
% figure();
% ax1 = subplot(2, 1, 1);
% plot(T_rtpd.TIME, T_rtpd.UP_RTPD, '-b', T_rtpd.TIME, -1.*T_rtpd.DOWN_RTPD, '-b', T_rtpd.TIME, T_rtpd.error_max, '-r', T_rtpd.TIME, T_rtpd.error_min, '-r');
% 
% ax2 = subplot(2, 1, 2);
% plot(T_rtpd.TIME, sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2), 'b');
% linkaxes([ax1(1), ax2],'x');

%% Load IBM data, forecast 15-min
month_ibm = [
    201908;
    201909;
    201910;
    201911;
    201912;
    202001;
    202002;
    202003;
    202004;
];

uniquegen     = {'gen55', 'gen56', 'gen59', 'gen60', 'gen61'};
uniquegensite = {'MNCC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz'};
uniquegenlat  = [34.31, 34.12, 37.41, 35.53, 35.38]; 
uniquegenlon  = [-117.50,-117.94,-119.74,-118.63,-120.18];

allgen  = {'gen55', 'gen56', 'gen57', 'gen58', 'gen59', 'gen60', 'gen61',    'gen62', 'gen63', 'gen64'};
allgensite = {'MNCC1', 'STFC1', 'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};
allgenlat = [34.31,34.12,34.12,34.12,37.41,35.53,35.38,34.31,34.31,35.53];
allgenlon = [-117.5,-117.94, -117.94, -117.94,-119.74, -118.63, -120.18,-117.5,-117.5,-118.63];

dirhome = pwd;
cell_pwr = cell(numel(uniquegen), 1);
dirwork = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\IBM\F_pwr.201908.15min';
for i = 1: numel(uniquegen)
    gen = uniquegen{i};
    for m = 1: numel(month_ibm)
        mstr = num2str(month_ibm(m));
        dirwork = strcat('C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\IBM\F_pwr.', mstr, '.15min');
        csvname = strcat(dirwork, '\', 'frcst_', gen, '.csv');
        cd(dirwork);
        if m == 1
            T_pwr = readtable(csvname, 'Delimiter', ',');
        else
            T_pwr = vertcat(T_pwr, readtable(csvname, 'Delimiter', ','));
        end
    end
    [~, ia, ~] = unique(T_pwr{:, {'Year', 'Month', 'Day', 'Hour', 'Minute'}}, 'rows');
    T_pwr = T_pwr(ia, :);
    T_pwr.TIME = datetime(T_pwr.Year, T_pwr.Month, T_pwr.Day, T_pwr.Hour, T_pwr.Minute, 0, 'TimeZone', 'UTC');
    
    % Calculate solar elevation to filter invalid clear-sky index
    Location = pvl_makelocationstruct(uniquegenlat(i), uniquegenlon(i));
    Time.year   = T_pwr.TIME.Year;
    Time.month  = T_pwr.TIME.Month;
    Time.day    = T_pwr.TIME.Day;
    Time.hour   = T_pwr.TIME.Hour;
    Time.minute = T_pwr.TIME.Minute;
    Time.second = T_pwr.TIME.Second;
    Time.UTCOffset = zeros(size(T_pwr, 1), 1);
    [SunAz, SunEl, ApparentSunEl, SolarTime] = pvl_ephemeris(Time, Location);
    T_pwr.ApparentSunEl = ApparentSunEl;
    T_pwr.SunEl = SunEl;
    T_pwr{T_pwr.ApparentSunEl<=3, 'ghi_cs'} = 0; % Consider only solar elevation > 3 degree
    T_pwr{T_pwr.ApparentSunEl<=3, 'pwr_cs'} = 0; % Consider only solar elevation > 3 degree
    
    % Calculate uncertainty bandwidth and variability of clear-sky index (k)
    T_pwr.k_p025 = T_pwr.ghi_p025./T_pwr.ghi_cs;
    T_pwr.k_p050 = T_pwr.ghi_p050./T_pwr.ghi_cs;
    T_pwr.k_p075 = T_pwr.ghi_p075./T_pwr.ghi_cs;
    T_pwr.k_width = T_pwr.k_p075 - T_pwr.k_p025;
    
    % Calculate uncertainty bandwidth and variability of clear-sky index of PV (k_PV)
    T_pwr.kpv_p025 = T_pwr.pwr_p025./T_pwr.pwr_cs;
    T_pwr.kpv_p050 = T_pwr.pwr_p050./T_pwr.pwr_cs;
    T_pwr.kpv_p075 = T_pwr.pwr_p075./T_pwr.pwr_cs;
    T_pwr.kpv_width = T_pwr.kpv_p075 - T_pwr.kpv_p025;
    
%     figure();
%     subplot(2, 1, 1);
%     hist(T_pwr.k_p050);
%     title('k');
%     subplot(2, 1, 2);
%     hist(T_pwr.kpv_p050);
%     title('kpv');
    
    cell_pwr{i} = T_pwr;
end
cd(dirhome);

%% Baseline

% Select month and k
% this_year = 2019;
% this_month = 10;
% karray = 5:5:60;

fprintf('Baseline\n');
% for k = 30
T_baseline_rtpd = T_rtpd((T_rtpd.HOUR_START_local.Month==this_month)&(T_rtpd.HOUR_START_local.Year==this_year), :);
[~, IA, IC ] = unique(T_baseline_rtpd.HOUR_START);
T_baseline_rtpd = T_baseline_rtpd(IA, {'HOUR_START', 'HOUR_START_local', 'DATE_START', 'DATE_START_local', 'UP_RTPD', 'DOWN_RTPD'});
for i = 1: size(T_baseline_rtpd, 1)
    this_date = T_baseline_rtpd.DATE_START(i);
    this_date_local = T_baseline_rtpd.DATE_START_local(i);
    this_hour_local = T_baseline_rtpd.HOUR_START_local.Hour(i);
    isweekday = (mod(weekday(this_date), 6) ~= 1);
    if isweekday
        k = 40;
    else
        k = 20;
    end
    selected_days = ismember(T_rtpd.DATE_START_local, return_history_days(this_date_local, k))&(T_rtpd.HOUR_START_local.Hour==this_hour_local);
    sample_error_max = T_rtpd{selected_days, 'error_max'};
    sample_error_min = T_rtpd{selected_days, 'error_min'};
    [f,x] = ecdf(sample_error_max(:));
    T_baseline_rtpd.FRU(i) = interp1(f, x, 0.975);
    [f,x] = ecdf(sample_error_min(:));
    T_baseline_rtpd.FRD(i) = interp1(f, x, 0.025);
end

% This is the actual need of FRP
f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
fru_need_rtpd = max(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4), [], 1)';
f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min	
frd_need_rtpd = min(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4), [], 1)';	

% Calculate baseline FRP imbalance
T_baseline_rtpd.FRU_error = T_baseline_rtpd.FRU - fru_need_rtpd;
T_baseline_rtpd.FRD_error = T_baseline_rtpd.FRD - frd_need_rtpd;

%% KNN Two-dim classifier, multi-site average, RTPD
karray = 5:5:40;
cell_results_rtpd  = cell(numel(karray), 1, 4); % k, site, classifier
cell_Tpwrhourly = cell(5, 1);
for s = 1: 5
    fprintf('s = %g\n', s);
    
    T_pwr = cell_pwr{s}; % The 5th is CA_Topaz site
    T_pwr.TIME_START  = T_pwr.TIME - duration(0, dt_rtpd, 0);
    T_pwr.HOUR_START  = datetime(T_pwr.TIME_START.Year, T_pwr.TIME_START.Month, T_pwr.TIME_START.Day, T_pwr.TIME_START.Hour, 0, 0, 'TimeZone', 'UTC');

    T_tmp1 = grpstats(T_pwr(:, {'HOUR_START', 'k_p050', 'kpv_p050', 'k_width', 'kpv_width'}), {'HOUR_START'}, 'mean');
    T_tmp2 = grpstats(T_pwr(:, {'HOUR_START', 'k_p050', 'kpv_p050', 'k_width', 'kpv_width'}), {'HOUR_START'}, 'std');
    T_pwr.dk = [nan; diff(T_pwr.k_p050)];
    T_pwr.dk_sq = [nan; diff(T_pwr.k_p050)].^2;
    T_pwr.dkpv = [nan; diff(T_pwr.kpv_p050)];
    T_pwr.dkpv_sq = [nan; diff(T_pwr.kpv_p050)].^2;
    T_pwr.dw = [nan; diff(T_pwr.k_width)];
    T_pwr.dw_sq = [nan; diff(T_pwr.dw)].^2;
    T_pwr.dwpv = [nan; diff(T_pwr.kpv_width)];
    T_pwr.dwpv_sq = [nan; diff(T_pwr.dwpv)].^2;
    T_tmp3 = grpstats(T_pwr(:, {'HOUR_START', 'dk_sq', 'dkpv_sq', 'dw_sq', 'dwpv_sq'}), {'HOUR_START'}, 'mean');
    T_tmp3.vk = sqrt(T_tmp3.mean_dk_sq);
    T_tmp3.vpv = sqrt(T_tmp3.mean_dkpv_sq);
    T_tmp3.vw = sqrt(T_tmp3.mean_dw_sq);
    T_tmp3.vwpv = sqrt(T_tmp3.mean_dwpv_sq);
    T_pwr_hourly = [T_tmp1(:, {'mean_k_p050', 'mean_kpv_p050', 'mean_k_width', 'mean_kpv_width'}) T_tmp2(:, {'std_k_p050', 'std_kpv_p050', 'std_k_width', 'std_kpv_width'}) T_tmp3(:, {'vk', 'vpv', 'vw', 'vwpv'})];
    T_pwr_hourly.HOUR_START = T_tmp1.HOUR_START;
    T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
    cell_Tpwrhourly{s} = T_pwr_hourly;
end
% Calculate 5-site mean
T_pwr_hourly_mean = T_pwr_hourly(:, {'HOUR_START', 'DATE'}); % Five site average
all_column = {'mean_k_p050', 'mean_kpv_p050', 'mean_k_width', 'mean_kpv_width', 'std_k_p050', 'std_kpv_p050', 'std_k_width', 'std_kpv_width', 'vk', 'vpv', 'vw', 'vwpv'};
array_pwr_hourly = nan(size(T_pwr_hourly_mean, 1), 12, 5);
for j = 1:5
    T_pwr_hourly = cell_Tpwrhourly{j};
    array_pwr_hourly(:, :, j) = T_pwr_hourly{:, all_column};
end
T_pwr_hourly_mean = [T_pwr_hourly_mean array2table(mean(array_pwr_hourly, 3), 'VariableNames', all_column)];
T_pwr_hourly_mean.HOUR_START_local = datetime(T_pwr_hourly_mean.HOUR_START, 'TimeZone', 'America/Los_Angeles');
T_pwr_hourly_mean.DATE_START_local = datetime(T_pwr_hourly_mean.HOUR_START_local.Year, T_pwr_hourly_mean.HOUR_START_local.Month, T_pwr_hourly_mean.HOUR_START_local.Day, 'TimeZone', 'America/Los_Angeles');

fprintf('kNN\n');

for s = 1

    % Select classifier
    for classifier = 1: 4
        T_pwr_hourly = T_pwr_hourly_mean;
        T_pwr_hourly.isweekday = (mod(weekday(T_pwr_hourly.DATE_START_local), 6) ~= 1);
        switch classifier
            case 1
                % Classifier 1: k (50 percentile), mean
                T_pwr_hourly.Properties.VariableNames{'mean_k_p050'} = 'classifier_1';
                T_pwr_hourly.Properties.VariableNames{'vk'} = 'classifier_2';
            case 2
                % Classifier 2: k (50 percentile), std.
                T_pwr_hourly.Properties.VariableNames{'mean_kpv_p050'} = 'classifier_1';
                T_pwr_hourly.Properties.VariableNames{'vpv'} = 'classifier_2';
            case 3
                % Classifier 3: k (50 percentile), variability
                T_pwr_hourly.Properties.VariableNames{'mean_k_p050'} = 'classifier_1';
                T_pwr_hourly.Properties.VariableNames{'mean_k_width'} = 'classifier_2';
            case 4
                % Classifier 4: k_pv (50 percentile), mean
                T_pwr_hourly.Properties.VariableNames{'mean_kpv_p050'} = 'classifier_1';
                T_pwr_hourly.Properties.VariableNames{'mean_kpv_width'} = 'classifier_2';
        end
        T_pwr_hourly.classifier_1(isinf(T_pwr_hourly.classifier_1)) = nan;
        T_pwr_hourly.classifier_2(isinf(T_pwr_hourly.classifier_2)) = nan;

        fprintf('classifier = %g\n', classifier);
        for k = karray
            % Test one month knn
            T_results_rtpd = T_pwr_hourly(T_pwr_hourly.DATE_START_local.Month==this_month, :); % Result container
            T_results_rtpd.use_knn = false(size(T_results_rtpd, 1), 1);
            % T_results_rtd = array2table(unique(T_rtd.HOUR_START(T_rtd.HOUR_START.Month==5)), 'VariableNames', {'HOUR_START'}); % Result container
            for i = 1: size(T_results_rtpd, 1)
                this_date = T_results_rtpd.DATE(i);
                this_date_local = T_results_rtpd.DATE_START_local(i);
                this_hour_local = T_results_rtpd.HOUR_START_local.Hour(i);
                isweekday = T_results_rtpd.isweekday(i);
                if any(isnan(T_results_rtpd{i, {'classifier_1', 'classifier_2'}}), 2)
                    if isweekday
                        ndays = 40;
                    else
                        ndays = 20;
                    end
                    selected_days = ismember(T_rtpd.DATE_START_local, return_history_days(this_date_local, ndays))&(T_rtpd.HOUR_START_local.Hour==this_hour_local); % We use 30 previous days
                else
                    T_sample = T_pwr_hourly((T_pwr_hourly.DATE_START_local<this_date_local)&(T_pwr_hourly.HOUR_START_local.Hour==this_hour_local)&(T_pwr_hourly.isweekday==isweekday), :);
                    T_sample.dist = sqrt(sum((T_sample{:, {'classifier_1', 'classifier_2'}}-T_results_rtpd{i, {'classifier_1', 'classifier_2'}}).^2, 2)); % Euclidean distance
                    T_sample_sorted = sortrows(T_sample, 'dist');
                    selected_days = ismember(T_rtpd.DATE_START_local, T_sample_sorted.DATE_START_local(1:k));
                    T_results_rtpd.use_knn(i) = true;
                end

                sample_error_max = T_rtpd{selected_days&(T_rtpd.HOUR_START_local.Hour==this_hour_local), 'error_max'};
                sample_error_min = T_rtpd{selected_days&(T_rtpd.HOUR_START_local.Hour==this_hour_local), 'error_min'};
                [f,x] = ecdf(sample_error_max(:));
                T_results_rtpd.FRU(i) = interp1(f, x, 0.975);
                [f,x] = ecdf(sample_error_min(:));
                T_results_rtpd.FRD(i) = interp1(f, x, 0.025);
            end

            T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
            T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

            cell_results_rtpd{karray==k, s, classifier} = T_results_rtpd;
            fprintf('k = %g\n', k);

        end
    end

end
%% Visualize the results, RTPD (Copied from RTD)

figure();
T_results_rtpd.HOUR  = T_results_rtpd.HOUR_START + duration(1, 0, 0); % We always use the end of an interval as time stamp
T_baseline_rtpd.HOUR = T_baseline_rtpd.HOUR_START + duration(1, 0, 0); % We always use the end of an interval as time stamp

h(1) = stairs(T_rtpd.TIME_START_local, T_rtpd.UP_RTPD, '-b');
hold on;
stairs(T_rtpd.TIME_START_local, -1.*T_rtpd.DOWN_RTPD, '-b');
h(2) = stairs(T_rtpd.TIME_START_local, T_rtpd.error_max, '-r');
stairs(T_rtpd.TIME_START_local, T_rtpd.error_min, '-r')
h(3) = stairs(T_results_rtpd.HOUR_START_local, T_results_rtpd.FRU, '-g');
stairs(T_results_rtpd.HOUR_START_local, T_results_rtpd.FRD, '-g');
h(4) = stairs(T_baseline_rtpd.HOUR_START_local, T_baseline_rtpd.FRU, '-k');
stairs(T_baseline_rtpd.HOUR_START_local, T_baseline_rtpd.FRD, '-k');
legend(h, {'OASIS', 'NL error', 'kNN', 'Baseline'});
ylabel('MW');

% Use box plot to show changes of FRP, RTD
fru_compare_rtpd = [T_baseline_rtpd{(T_baseline_rtpd.HOUR_START.Hour>=16)&(T_baseline_rtpd.HOUR_START.Hour<=24), {'FRU'}} T_results_rtpd{(T_results_rtpd.HOUR_START.Hour>=16)&(T_results_rtpd.HOUR_START.Hour<=24), {'FRU'}}];
delta_fru_percent = (fru_compare_rtpd(:, 2)-fru_compare_rtpd(:, 1))./fru_compare_rtpd(:, 1);

figure();
boxplot(reshape(delta_fru_percent, 8, size(delta_fru_percent, 1)/8)', 'Label', {'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
xlabel('Time (PST)');
ylabel('Relative change to baseline');
title('FRU');

frd_compare_rtpd = [T_baseline_rtpd{(T_baseline_rtpd.HOUR_START.Hour>=16)&(T_baseline_rtpd.HOUR_START.Hour<=24), {'FRD'}} T_results_rtpd{(T_results_rtpd.HOUR_START.Hour>=16)&(T_results_rtpd.HOUR_START.Hour<=24), {'FRD'}}];
delta_frd_percent = -(frd_compare_rtpd(:, 2)-frd_compare_rtpd(:, 1))./frd_compare_rtpd(:, 1);

figure();
% boxplot(reshape(delta_frd_percent.*100, 8, size(delta_frd_percent, 1)/8)', 'Label', {'16', '17', '18', '19', '20', '21', '22', '23'});
boxplot(reshape(delta_frd_percent, 8, size(delta_frd_percent, 1)/8)', 'Label', {'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
xlabel('Time (PST)');
ylabel('Relative change to baseline');
title('FRD');

fru_need_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START)&(T_rtpd.HOUR_START.Hour>=16)&(T_rtpd.HOUR_START.Hour<=24), 'error_max'};
fru_need_rtpd = reshape(fru_need_rtpd, 4, size(fru_need_rtpd, 1)/4)';

fru_shortage_baseline_rtpd = fru_need_rtpd - fru_compare_rtpd(:, 1);
fru_shortage_results_rtpd  = fru_need_rtpd - fru_compare_rtpd(:, 2);
fru_n_shortage_hourly_results_rtpd = sum(reshape(sum(fru_shortage_results_rtpd>0, 2), 8, size(fru_shortage_results_rtpd, 1)/8), 2);
fru_n_not_nan_hourly_results_rtpd  = sum(reshape(~isnan(fru_shortage_results_rtpd), 8, numel(fru_shortage_results_rtpd)/8), 2);
fru_n_shortage_hourly_baseline_rtpd = sum(reshape(sum(fru_shortage_baseline_rtpd>0, 2), 8, size(fru_shortage_baseline_rtpd, 1)/8), 2);
fru_n_not_nan_hourly_baseline_rtpd  = sum(reshape(~isnan(fru_shortage_baseline_rtpd), 8, numel(fru_shortage_baseline_rtpd)/8), 2);
figure(); 
bar([fru_n_shortage_hourly_results_rtpd./fru_n_not_nan_hourly_results_rtpd fru_n_shortage_hourly_baseline_rtpd./fru_n_not_nan_hourly_baseline_rtpd]); % Percentage of violations each hour
legend({'New', 'Baseline'});
set(gca,'xticklabel',{'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
ylabel('Time of FRP shortage');
title('FRU');

frd_need_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START)&(T_rtpd.HOUR_START.Hour>=16)&(T_rtpd.HOUR_START.Hour<=24), 'error_min'};
frd_need_rtpd = reshape(frd_need_rtpd, 4, size(frd_need_rtpd, 1)/4)';

frd_shortage_baseline_rtpd = frd_need_rtpd - frd_compare_rtpd(:, 1);
frd_shortage_results_rtpd  = frd_need_rtpd - frd_compare_rtpd(:, 2);
frd_n_shortage_hourly_results_rtpd = sum(reshape(sum(frd_shortage_results_rtpd<0, 2), 8, size(frd_shortage_results_rtpd, 1)/8), 2);
frd_n_not_nan_hourly_results_rtpd  = sum(reshape(~isnan(frd_shortage_results_rtpd), 8, numel(frd_shortage_results_rtpd)/8), 2);
frd_n_shortage_hourly_baseline_rtpd = sum(reshape(sum(frd_shortage_baseline_rtpd<0, 2), 8, size(frd_shortage_baseline_rtpd, 1)/8), 2);
frd_n_not_nan_hourly_baseline_rtpd  = sum(reshape(~isnan(frd_shortage_baseline_rtpd), 8, numel(frd_shortage_baseline_rtpd)/8), 2);
figure(); 
bar([frd_n_shortage_hourly_results_rtpd./frd_n_not_nan_hourly_results_rtpd frd_n_shortage_hourly_baseline_rtpd./frd_n_not_nan_hourly_baseline_rtpd]); % Percentage of violations each hour
legend({'New', 'Baseline'});
set(gca,'xticklabel',{'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
ylabel('Time of FRP shortage');
title('FRD');

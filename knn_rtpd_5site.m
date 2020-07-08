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
T_rtpd.local_time = datetime(T_rtpd.TIME, 'TimeZone', 'America/Los_Angeles');

% calculate pure solar forecast error
tmp_rtp_solar  = sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2);
tmp_rtpd_solar = sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2);
tmp_error_solar = -reshape(tmp_rtp_solar, 3, size(tmp_rtp_solar, 1)/3)' - (-repmat(tmp_rtpd_solar, 1, 3));

% Calculate net load forecast error
% tmp_rtd_nl_b  = T_rtd.LOAD_B_RTD - sum(T_rtd{:, {'Wind_NP15_B_RTD', 'Wind_SP15_B_RTD', 'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2);
% tmp_rtpd_nl_a = T_rtpd.LOAD_B_RTPD - sum(T_rtpd{:, {'Wind_NP15_RTPD', 'Wind_SP15_RTPD', 'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2); % Note we use binding forecast as advisory forecast since CAISO does not publish advisory forecast
tmp_rtd_nl_b  = - sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2);
tmp_rtpd_nl_a = - sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2); % Note we use binding forecast as advisory forecast since CAISO does not publish advisory forecast

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

%% One-dim classifier, multi-site average, RTPD

% Select month and k
% this_year = 2019;
% this_month = 10;
karray = 5:5:60;

% Result container
cell_baseline_rtpd = cell(numel(karray), 1); % site, k, classifier
cell_results_rtpd  = cell(numel(karray), 1, 12); % site, k, classifier

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

for s = 1
    fprintf('Baseline\n');
    for k = karray
        T_baseline_rtpd = array2table(unique(T_rtpd.HOUR_START((T_rtpd.HOUR_START.Month==this_month)&(T_rtpd.HOUR_START.Year==this_year))), 'VariableNames', {'HOUR_START'}); % Result container
        for i = 1: size(T_baseline_rtpd, 1)
            this_date = datetime(T_baseline_rtpd.HOUR_START.Year(i), T_baseline_rtpd.HOUR_START.Month(i), T_baseline_rtpd.HOUR_START.Day(i), 'TimeZone', 'UTC');
            this_hour = T_baseline_rtpd.HOUR_START.Hour(i);
            selected_days = (T_rtpd.TIME_START<this_date)&(T_rtpd.TIME_START>=this_date-days(k))&(T_rtpd.HOUR_START.Hour==this_hour); % We use 30 previous days
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
        
        cell_baseline_rtpd{karray==k, s} = T_baseline_rtpd;

        fprintf('k = %g\n', k);
    end
    
    fprintf('kNN\n');
    % Select classifier
    for classifier = 1: 12
        T_pwr_hourly = T_pwr_hourly_mean;
        switch classifier
            case 1
                T_pwr_hourly.Properties.VariableNames{'mean_k_p050'} = 'classifier_1';
            case 2
                T_pwr_hourly.Properties.VariableNames{'std_k_p050'} = 'classifier_1';
            case 3
                T_pwr_hourly.Properties.VariableNames{'vk'} = 'classifier_1';
            case 4
                T_pwr_hourly.Properties.VariableNames{'mean_kpv_p050'} = 'classifier_1';
            case 5
                T_pwr_hourly.Properties.VariableNames{'std_kpv_p050'} = 'classifier_1';
            case 6
                T_pwr_hourly.Properties.VariableNames{'vpv'} = 'classifier_1';
            case 7
                T_pwr_hourly.Properties.VariableNames{'mean_k_width'} = 'classifier_1';
            case 8
                T_pwr_hourly.Properties.VariableNames{'std_k_width'} = 'classifier_1';
            case 9
                T_pwr_hourly.Properties.VariableNames{'vw'} = 'classifier_1';
            case 10
                T_pwr_hourly.Properties.VariableNames{'mean_kpv_width'} = 'classifier_1';
            case 11
                T_pwr_hourly.Properties.VariableNames{'std_kpv_width'} = 'classifier_1';
            case 12
                T_pwr_hourly.Properties.VariableNames{'vwpv'} = 'classifier_1';
        end

        fprintf('classifier = %g\n', classifier);
        for k = karray
            % Test one month knn
            T_results_rtpd = T_pwr_hourly(T_pwr_hourly.DATE.Month==this_month, :); % Result container
            % T_results_rtd = array2table(unique(T_rtd.HOUR_START(T_rtd.HOUR_START.Month==5)), 'VariableNames', {'HOUR_START'}); % Result container
            for i = 1: size(T_results_rtpd, 1)
                this_date = datetime(T_results_rtpd.HOUR_START.Year(i), T_results_rtpd.HOUR_START.Month(i), T_results_rtpd.HOUR_START.Day(i), 'TimeZone', 'UTC');
                this_hour = T_results_rtpd.HOUR_START.Hour(i);
                T_sample = T_pwr_hourly((T_pwr_hourly.DATE<this_date)&(T_pwr_hourly.HOUR_START.Hour==this_hour), :);
                if any(isnan(T_results_rtpd{i, 'classifier_1'}))	
                    T_sample_sorted = T_sample(T_sample.DATE>=this_date-days(k), :); % We use 30 previous days, i.e., baseline, if data is nan
                    selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); % 30 the nearest days for baseline
                else
                    T_sample.dist = sqrt(sum((T_sample{:, 'classifier_1'}-T_results_rtpd{i, 'classifier_1'}).^2, 2)); % Euclidean distance
                    T_sample_sorted = sortrows(T_sample, 'dist');
                    selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); 
                end

                sample_error_max = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_max'};
                sample_error_min = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_min'};
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
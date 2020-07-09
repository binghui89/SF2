dt_rtd = 5; % min
dt_rtpd = 15; % min

%% Load NDA data
T_nda = readtable('C:\Users\bxl180002\Downloads\Transformed_NDA_errors_BL_es.xlsx');
T_nda.TIME_START = datetime(T_nda.timestart_full, 'InputFormat', 'yyyy-MM-dd''T''HH:mmXXX', 'TimeZone', 'UTC');
T_nda.TIME = T_nda.TIME_START + duration(0, dt_rtpd, 0);
T_nda.HOUR_START = datetime(T_nda.TIME_START.Year, T_nda.TIME_START.Month, T_nda.TIME_START.Day, T_nda.TIME_START.Hour, 0, 0, 'TimeZone', 'UTC');
T_nda.DATE_START = datetime(T_nda.HOUR_START.Year, T_nda.HOUR_START.Month, T_nda.HOUR_START.Day, 'TimeZone', 'UTC');
T_nda.TIME_local = datetime(T_nda.TIME, 'TimeZone', 'America/Los_Angeles');
T_nda.TIME_START_local = datetime(T_nda.TIME_START, 'TimeZone', 'America/Los_Angeles');
T_nda.HOUR_START_local = datetime(T_nda.HOUR_START, 'TimeZone', 'America/Los_Angeles');
T_nda.DATE_START_local = datetime(T_nda.HOUR_START_local.Year, T_nda.HOUR_START_local.Month, T_nda.HOUR_START_local.Day, 'TimeZone', 'America/Los_Angeles');

T_nda.error_max = T_nda.transformedRampUp;
T_nda.error_min = -T_nda.transformedRampDown;

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

%% Two-dim classifier, multi-site average, RTPD

% Select month and k
% this_year = 2019;
% this_month = 10;
% karray = 5:5:60;
karray = 30;

% Result container
cell_baseline_rtpd = cell(numel(karray), 1); % k, site, classifier
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

for s = 1
        fprintf('Baseline\n');
    for k = karray
%         T_baseline_rtpd = array2table(unique(T_nda.HOUR_START((T_nda.HOUR_START_local.Month==this_month)&(T_nda.HOUR_START_local.Year==this_year))), 'VariableNames', {'HOUR_START'}); % Result container
%         T_baseline_rtpd.HOUR_START_local = datetime(T_baseline_rtpd.HOUR_START, 'TimeZone', 'America/Los_Angeles');
        T_baseline_rtpd = T_nda((T_nda.HOUR_START_local.Month==this_month)&(T_nda.HOUR_START_local.Year==this_year), :);
        [~, IA, IC ] = unique(T_baseline_rtpd.HOUR_START);
        T_baseline_rtpd = T_baseline_rtpd(IA, {'HOUR_START', 'HOUR_START_local', 'DATE_START', 'DATE_START_local'});
        for i = 1: size(T_baseline_rtpd, 1)
            this_date = T_baseline_rtpd.DATE_START(i);
            this_date_local = T_baseline_rtpd.DATE_START_local(i);
            this_hour_local = T_baseline_rtpd.HOUR_START_local.Hour(i);
            selected_days = ismember(T_nda.DATE_START_local, return_history_days(this_date_local, k))&(T_nda.HOUR_START_local.Hour==this_hour_local); % We use 30 previous days
            sample_error_max = T_nda{selected_days, 'error_max'};
            sample_error_min = T_nda{selected_days, 'error_min'};
            [f,x] = ecdf(sample_error_max(:));
            T_baseline_rtpd.FRU(i) = interp1(f, x, 0.975);
            [f,x] = ecdf(sample_error_min(:));
            T_baseline_rtpd.FRD(i) = interp1(f, x, 0.025);
        end

        % This is the actual need of FRP
        f_errormax_rtpd = T_nda{ismember(T_nda.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
        fru_need_rtpd = max(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4), [], 1)';
        f_errormin_rtpd = T_nda{ismember(T_nda.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min	
        frd_need_rtpd = min(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4), [], 1)';	

        % Calculate baseline FRP imbalance
        T_baseline_rtpd.FRU_error = T_baseline_rtpd.FRU - fru_need_rtpd;
        T_baseline_rtpd.FRD_error = T_baseline_rtpd.FRD - frd_need_rtpd;
        
        cell_baseline_rtpd{karray==k, s} = T_baseline_rtpd;

        fprintf('k = %g\n', k);
    end
    
    fprintf('kNN\n');
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
                if any(isnan(T_results_rtpd{i, {'classifier_1', 'classifier_2'}}), 2)
                    selected_days = ismember(T_nda.DATE_START_local, return_history_days(this_date_local, k))&(T_nda.HOUR_START_local.Hour==this_hour_local); % We use 30 previous days
                else
                    isweekday = T_results_rtpd.isweekday(i);
                    T_sample = T_pwr_hourly((T_pwr_hourly.DATE_START_local<this_date_local)&(T_pwr_hourly.HOUR_START_local.Hour==this_hour_local)&(T_pwr_hourly.isweekday==isweekday), :);
                    T_sample.dist = sqrt(sum((T_sample{:, {'classifier_1', 'classifier_2'}}-T_results_rtpd{i, {'classifier_1', 'classifier_2'}}).^2, 2)); % Euclidean distance
                    T_sample_sorted = sortrows(T_sample, 'dist');
                    selected_days = ismember(T_nda.DATE_START_local, T_sample_sorted.DATE_START_local(1:k));
                    T_results_rtpd.use_knn(i) = true;
                end

                sample_error_max = T_nda{selected_days&(T_nda.HOUR_START_local.Hour==this_hour_local), 'error_max'};
                sample_error_min = T_nda{selected_days&(T_nda.HOUR_START_local.Hour==this_hour_local), 'error_min'};
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

%% Save csv file
T_results_rtpd = cell_results_rtpd{karray==30, 1, 1}; % Use the first classifier
T_save = T_nda((T_nda.HOUR_START_local<=max(T_results_rtpd.HOUR_START_local))&(T_nda.HOUR_START_local>=min(T_results_rtpd.HOUR_START_local)), {'timestamp', 'hour', 'RTPDInterval', 'transformedRampUp', 'transformedRampDown', 'offset_negative', 'timestart_full'});
T_save.isweekday = reshape(repmat(T_results_rtpd.isweekday, 1, 4)', size(T_results_rtpd, 1)*4, []);
T_save.use_knn = reshape(repmat(T_results_rtpd.use_knn, 1, 4)', size(T_results_rtpd, 1)*4, []);
T_save.FRU_baseline = reshape(repmat(T_baseline_rtpd.FRU, 1, 4)', size(T_baseline_rtpd, 1)*4, []);
T_save.FRD_baseline = -reshape(repmat(T_baseline_rtpd.FRD, 1, 4)', size(T_baseline_rtpd, 1)*4, []);
T_save.FRU_knn = reshape(repmat(T_results_rtpd.FRU, 1, 4)', size(T_results_rtpd, 1)*4, []);
T_save.FRD_knn = -reshape(repmat(T_results_rtpd.FRD, 1, 4)', size(T_results_rtpd, 1)*4, []);
xlsxname = strcat('NDA_results_2020', num2str(this_month, '%02g'), '.xlsx');
writetable(T_save,xlsxname,'Sheet',1);

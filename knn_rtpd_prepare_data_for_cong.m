%% Load CAISO data RTD and RTPD
load_caiso_data;

%% Load IBM data, forecast 15-min
load_ibm_5sites;

%% Prepare and save data for Cong
ar_datetime = datetime(2019,08,01, 0, 0, 0, 'TimeZone', 'UTC-8'):duration(1, 0, 0): datetime(2020,04,29, 23, 0, 0, 'TimeZone', 'UTC-8');
ar_datetime = ar_datetime(:);
ar_datetime.TimeZone = 'UTC';

cell_pwr_hourly = cell(size(cell_pwr));
cell_pwr_hourly_forcong = cell(size(cell_pwr));
for s = 1: 5
    fprintf('s = %g\n', s);
    
    T_pwr = cell_pwr{s}; % The 5th is CA_Topaz site
    T_pwr.TIME_START  = T_pwr.TIME - duration(0, dt_rtpd, 0);
    T_pwr.HOUR_START  = datetime(T_pwr.TIME_START.Year, T_pwr.TIME_START.Month, T_pwr.TIME_START.Day, T_pwr.TIME_START.Hour, 0, 0, 'TimeZone', 'UTC');
    
    T_pwr_hourly = table(unique(T_pwr.HOUR_START), 'VariableNames', {'HOUR_START'});
    
    % % Classifier 1: k (50 percentile), mean
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'k_p050'}), {'HOUR_START'}, 'mean'); 
    T_pwr_hourly.mean_k = tmp.mean_k_p050;

    % % Classifier 2: k (50 percentile), std.
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'k_p050'}), {'HOUR_START'}, 'std'); 
    T_pwr_hourly.std_k = tmp.std_k_p050;

    % % Classifier 3: k (50 percentile), variability
    T_pwr.dk = [nan; diff(T_pwr.k_p050)]; % delta k
    T_pwr.dk_sq = [nan; diff(T_pwr.k_p050)].^2; % % (delta k)^2
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'dk_sq', 'k_width'}), {'HOUR_START'}, 'mean'); 
    tmp.v = sqrt(tmp.mean_dk_sq); % This is variability within each hour, following Inman et al. 2013, section 2.6.1
    T_pwr_hourly.v_k = tmp.v;

    % % Classifier 4: k_pv (50 percentile), mean
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'kpv_p050'}), {'HOUR_START'}, 'mean'); 
    T_pwr_hourly.mean_kpv = tmp.mean_kpv_p050;

    % % Classifier 5: k_pv (50 percentile), std
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'kpv_p050'}), {'HOUR_START'}, 'std'); 
    T_pwr_hourly.std_kpv = tmp.std_kpv_p050;

    % % Classifier 6: k (50 percentile), variability
    T_pwr.dkpv = [nan; diff(T_pwr.kpv_p050)]; % delta k
    T_pwr.dkpv_sq = [nan; diff(T_pwr.kpv_p050)].^2; % % (delta k)^2
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'dkpv_sq'}), {'HOUR_START'}, 'mean'); 
    tmp.vpv = sqrt(tmp.mean_dkpv_sq); % This is variability within each hour, following Inman et al. 2013, section 2.6.1
    T_pwr_hourly.v_kpv = tmp.vpv;

    % % Classifier 7: width of k (75 - 25 percentile), mean
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'k_width'}), {'HOUR_START'}, 'mean'); 
    T_pwr_hourly.mean_w = tmp.mean_k_width;

    % % Classifier 8: width of k (75 - 25 percentile), std
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'k_width'}), {'HOUR_START'}, 'std'); 
    T_pwr_hourly.std_w = tmp.std_k_width;

    % % Classifier 9: width of k (75 - 25 percentile), variability
    T_pwr.dw = [nan; diff(T_pwr.k_width)]; % delta k width
    T_pwr.dw_sq = [nan; diff(T_pwr.dw)].^2; % % (delta k)^2
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'dw_sq'}), {'HOUR_START'}, 'mean'); 
    tmp.vw = sqrt(tmp.mean_dw_sq); % This is variability within each hour, following Inman et al. 2013, section 2.6.1
    T_pwr_hourly.v_w = tmp.vw;

    % % Classifier 10: width of kpv (75 - 25 percentile), mean
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'kpv_width'}), {'HOUR_START'}, 'mean'); 
    T_pwr_hourly.mean_wpv = tmp.mean_kpv_width;

    % % Classifier 11: width of k (75 - 25 percentile), std
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'kpv_width'}), {'HOUR_START'}, 'std'); 
    T_pwr_hourly.std_wpv = tmp.std_kpv_width;

    % % Classifier 12: width of k (75 - 25 percentile), variability
    T_pwr.dwpv = [nan; diff(T_pwr.kpv_width)]; % delta k width
    T_pwr.dwpv_sq = [nan; diff(T_pwr.dwpv)].^2; % % (delta k)^2
    tmp = grpstats(T_pwr(:, {'HOUR_START', 'dwpv_sq'}), {'HOUR_START'}, 'mean'); 
    tmp.vwpv = sqrt(tmp.mean_dwpv_sq); % This is variability within each hour, following Inman et al. 2013, section 2.6.1
    T_pwr_hourly.v_wpv = tmp.vwpv;
    
    cell_pwr_hourly{s} = T_pwr_hourly;
    T_pwr_hourly_forcong = T_pwr_hourly((T_pwr_hourly{:, 'HOUR_START'}>=min(ar_datetime))&(T_pwr_hourly{:, 'HOUR_START'}<=max(ar_datetime)), :);
    cell_pwr_hourly_forcong{s} = T_pwr_hourly_forcong;

end

dirhome = pwd;
T_rtpd_forcong = T_rtpd((T_rtpd{:, 'HOUR_START'}>=min(ar_datetime))&(T_rtpd{:, 'HOUR_START'}<=max(ar_datetime)), {'HOUR_START', 'error_max', 'error_min'}); % Only look at the hours with IBM forecast
T_rtpd_forcong = grpstats(T_rtpd_forcong, 'HOUR_START', {'max', 'min'}, 'DataVars', {'error_max', 'error_min'});
T_rtpd_forcong = T_rtpd_forcong(:, {'HOUR_START', 'max_error_max', 'min_error_min'});
T_rtpd_forcong.Properties.VariableNames{'max_error_max'} = 'FRU';
T_rtpd_forcong.Properties.VariableNames{'min_error_min'} = 'FRD';

if write_flag
    cd('C:\Users\bxl180002\Downloads\RampSolar');
    writetable(T_rtpd_forcong, 'rtpd.csv');
    for s = 1:5
        writetable(cell_pwr_hourly_forcong{s}, strcat('T_site_', num2str(s),'.csv'));
    end
    cd(dirhome);

end

T_rtpd_forcong = T_rtpd((T_rtpd{:, 'HOUR_START'}>=min(ar_datetime))&(T_rtpd{:, 'HOUR_START'}<=max(ar_datetime)), {'HOUR_START', 'error_max', 'error_min'});
T_rtpd_forcong = T_rtpd_forcong(:, {'HOUR_START', 'error_max', 'error_min'});
T_rtpd_forcong = grpstats(T_rtpd_forcong, 'HOUR_START', {'max', 'min'}, 'DataVars', {'error_max', 'error_min'});
T_rtpd_forcong = T_rtpd_forcong(:, {'HOUR_START', 'max_error_max', 'min_error_min'});
T_rtpd_forcong.Properties.VariableNames{'max_error_max'} = 'FRU';
T_rtpd_forcong.Properties.VariableNames{'min_error_min'} = 'FRD';
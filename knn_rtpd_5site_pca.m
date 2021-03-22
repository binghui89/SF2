% Use PCA to reduce dimension

%% Load CAISO data RTD and RTPD
load_caiso_data;

%% Load IBM data, forecast 15-min
load_ibm_5sites;


%% Prepare and save data for Cong
% Fix the errors in the calculation of classifiers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    T_pwr.k_p025(isinf(T_pwr.k_p025)) = nan;
    T_pwr.k_p050(isinf(T_pwr.k_p050)) = nan;
    T_pwr.k_p075(isinf(T_pwr.k_p075)) = nan;
    T_pwr.kpv_p025(isnan(T_pwr.k_p025)) = nan;
    T_pwr.kpv_p050(isnan(T_pwr.k_p050)) = nan;
    T_pwr.kpv_p075(isnan(T_pwr.k_p075)) = nan;
    T_pwr.kpv_p025(isinf(T_pwr.kpv_p025)) = nan;
    T_pwr.kpv_p050(isinf(T_pwr.kpv_p050)) = nan;
    T_pwr.kpv_p075(isinf(T_pwr.kpv_p075)) = nan;
    
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
    % T_pwr.dw_sq = [nan; diff(T_pwr.dw)].^2; % % (delta k)^2, this is not correct
    T_pwr.dw_sq = T_pwr.dw.^2; % % (delta k)^2
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
    % T_pwr.dwpv_sq = [nan; diff(T_pwr.dwpv)].^2; % % (delta k)^2, this is not correct
    T_pwr.dwpv_sq = T_pwr.dwpv.^2; % % (delta k)^2
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

write_flag = false;
if write_flag
    cd('C:\Users\bxl180002\Downloads\RampSolar');
    writetable(T_rtpd_forcong, 'rtpd.csv');
    for s = 1:5
        writetable(cell_pwr_hourly_forcong{s}, strcat('T_site_', num2str(s),'.csv'));
    end
    cd(dirhome);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data cleansing
colnames = {'mean_k','std_k','v_k','mean_kpv','std_kpv','v_kpv','mean_w','std_w','v_w','mean_wpv','std_wpv','v_wpv'};
cell_pwr_hourly_norm = cell(size(cell_pwr_hourly_forcong));
for s = 1: length(cell_pwr_hourly_forcong)
    T_pwr_hourly_norm = cell_pwr_hourly_forcong{s};
    
    for ic = 1: length(colnames)
        T_pwr_hourly_norm{isinf(T_pwr_hourly_norm{:, colnames{ic}}), colnames{ic}} = nan;
    end

    T_pwr_hourly_norm.mean_k(isinf(T_pwr_hourly_norm.mean_k))= nan;
    T_pwr_hourly_norm.mean_k = (T_pwr_hourly_norm.mean_k - mean(T_pwr_hourly_norm.mean_k, 'omitnan'))./std(T_pwr_hourly_norm.mean_k, 'omitnan');
    T_pwr_hourly_norm.std_k = (T_pwr_hourly_norm.std_k - mean(T_pwr_hourly_norm.std_k, 'omitnan'))./std(T_pwr_hourly_norm.std_k, 'omitnan');
    T_pwr_hourly_norm.v_k(isinf(T_pwr_hourly_norm.v_k)) = nan;
    T_pwr_hourly_norm.v_k = (T_pwr_hourly_norm.v_k - mean(T_pwr_hourly_norm.v_k, 'omitnan'))./std(T_pwr_hourly_norm.v_k, 'omitnan');

    T_pwr_hourly_norm.mean_kpv(isinf(T_pwr_hourly_norm.mean_kpv))=nan;
    T_pwr_hourly_norm.mean_kpv = (T_pwr_hourly_norm.mean_kpv - mean(T_pwr_hourly_norm.mean_kpv, 'omitnan'))./std(T_pwr_hourly_norm.mean_kpv, 'omitnan');
    T_pwr_hourly_norm.std_kpv = (T_pwr_hourly_norm.std_kpv - mean(T_pwr_hourly_norm.std_kpv, 'omitnan'))./std(T_pwr_hourly_norm.std_kpv, 'omitnan');
    T_pwr_hourly_norm.v_kpv = (T_pwr_hourly_norm.v_kpv - mean(T_pwr_hourly_norm.v_kpv, 'omitnan'))./std(T_pwr_hourly_norm.v_kpv, 'omitnan');

    T_pwr_hourly_norm.mean_w(T_pwr_hourly_norm.mean_w<0) = nan;
    T_pwr_hourly_norm.mean_w = (T_pwr_hourly_norm.mean_w - mean(T_pwr_hourly_norm.mean_w, 'omitnan'))./std(T_pwr_hourly_norm.mean_w, 'omitnan');
    T_pwr_hourly_norm.std_w = (T_pwr_hourly_norm.std_w - mean(T_pwr_hourly_norm.std_w, 'omitnan'))./std(T_pwr_hourly_norm.std_w, 'omitnan');
    T_pwr_hourly_norm.v_w = (T_pwr_hourly_norm.v_w - mean(T_pwr_hourly_norm.v_w, 'omitnan'))./std(T_pwr_hourly_norm.v_w, 'omitnan');

    T_pwr_hourly_norm.mean_wpv(T_pwr_hourly_norm.mean_wpv<0) = nan;
    T_pwr_hourly_norm.mean_wpv = (T_pwr_hourly_norm.mean_wpv - mean(T_pwr_hourly_norm.mean_wpv, 'omitnan'))./std(T_pwr_hourly_norm.mean_wpv, 'omitnan');
    T_pwr_hourly_norm.std_wpv = (T_pwr_hourly_norm.std_wpv - mean(T_pwr_hourly_norm.std_wpv, 'omitnan'))./std(T_pwr_hourly_norm.std_wpv, 'omitnan');
    T_pwr_hourly_norm.v_wpv = (T_pwr_hourly_norm.v_wpv - mean(T_pwr_hourly_norm.v_wpv, 'omitnan'))./std(T_pwr_hourly_norm.v_wpv, 'omitnan');

    figure();
    subplot(4, 3, 1);
    histogram(T_pwr_hourly_norm.mean_k, 100); 
    subplot(4, 3, 2);
    histogram(T_pwr_hourly_norm.std_k, 100);
    subplot(4, 3, 3);
    histogram(T_pwr_hourly_norm.v_k, 100);

    subplot(4, 3, 4);
    histogram(T_pwr_hourly_norm.mean_kpv, 100);
    subplot(4, 3, 5);
    histogram(T_pwr_hourly_norm.std_kpv, 100);
    subplot(4, 3, 6);
    histogram(T_pwr_hourly_norm.v_kpv, 100);

    subplot(4, 3, 7);
    histogram(T_pwr_hourly_norm.mean_w, 100);
    subplot(4, 3, 8);
    histogram(T_pwr_hourly_norm.std_w, 100);
    subplot(4, 3, 9);
    histogram(T_pwr_hourly_norm.v_w, 100);

    subplot(4, 3, 10);
    histogram(T_pwr_hourly_norm.mean_wpv, 100);
    subplot(4, 3, 11);
    histogram(T_pwr_hourly_norm.std_wpv, 100);
    subplot(4, 3, 12);
    histogram(T_pwr_hourly_norm.v_wpv, 100);
    
    cell_pwr_hourly_norm{s} = T_pwr_hourly_norm;
end

%% PCA, Multi-site
% this_year = 2020;
% this_month = 2;
ar_all = [ cell_pwr_hourly_norm{1}{:, colnames} cell_pwr_hourly_norm{2}{:, colnames} cell_pwr_hourly_norm{3}{:, colnames} cell_pwr_hourly_norm{4}{:, colnames} cell_pwr_hourly_norm{5}{:, colnames}];
[coeff,score,latent,tsquared,explained,mu] = pca(ar_all);
T_score_5site = table(T_pwr_hourly_norm.HOUR_START, score, 'VariableNames', {'HOUR_START', 'score'});
figure(); pareto(explained);

karray = 5:5:60;

% Result container
cell_results_rtpd  = cell(numel(karray), 1, 4); % site, k, classifier

%% kNN, multi-site PCA
fprintf('kNN, 5-site PCA\n');
% Select classifier
T_socre = T_score_5site;
s = 1;
for classifier = 1: 4 % # of PCA components
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

    T_pwr_hourly = T_socre;
    T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
    fprintf('classifier = %g\n', classifier);
    w = ones(classifier, 1);
    for k = karray
        % Test one month knn
        T_results_rtpd = T_pwr_hourly(T_pwr_hourly.DATE.Month==this_month, :); % Result container
        for i = 1: size(T_results_rtpd, 1)
            this_date = datetime(T_results_rtpd.HOUR_START.Year(i), T_results_rtpd.HOUR_START.Month(i), T_results_rtpd.HOUR_START.Day(i), 'TimeZone', 'UTC');
            this_hour = T_results_rtpd.HOUR_START.Hour(i);
            T_sample = T_pwr_hourly((T_pwr_hourly.DATE<this_date)&(T_pwr_hourly.HOUR_START.Hour==this_hour), :);
            if any(isnan(T_results_rtpd.score(i, 1:classifier)), 2)
                T_sample_sorted = T_sample(T_sample.DATE>=this_date-days(k), :); % We use 30 previous days, i.e., baseline, if data is nan
                selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); % 30 the nearest days for baseline
            else
                T_sample.dist = abs(T_sample.score(:, 1:classifier)-T_results_rtpd.score(i, 1:classifier))*w; % Euclidean distance
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


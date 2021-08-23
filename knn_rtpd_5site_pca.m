% Use PCA to reduce dimension


%% Prepare data and write to disk
if ~exist('write_flag')
    write_flag = input('Write file (Yes: 1/No: 0)?');
end
if ~exist('this_year')
    this_year  = input('Year (2019 or 2020)?');
end
if ~exist('this_month')
    this_month = input('Month (1 to 12)?');
end
knn_rtpd_prepare_data_for_cong;

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

if write_flag
    T_score_forcong = T_score_5site(:, 'HOUR_START');
    T_score_forcong{:, '1'} = T_score_5site.score(:, 1);
    T_score_forcong{:, '2'} = T_score_5site.score(:, 2);
    T_score_forcong{:, '3'} = T_score_5site.score(:, 3);
    T_score_forcong{:, '4'} = T_score_5site.score(:, 4);
    
    cd('C:\Users\bxl180002\Downloads\RampSolar');
    writetable(T_score_forcong, 'pca_score.csv');
    cd(dirhome);
end

karray = 5:5:60;

% Result container
cell_results_rtpd  = cell(numel(karray), 1, 4); % site, k, classifier

%% kNN, multi-site PCA
fprintf('kNN, 5-site PCA\n');
% Select classifier
T_socre = T_score_5site;
s = 1;
for classifier = 1: 4 % # of PCA components
    if classifier == 1 % We only need to run baseline once
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
    end

    T_pwr_hourly = T_socre;
    T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
    fprintf('classifier = %g\n', classifier);
    w = ones(classifier, 1);
    for k = karray
        % Test one month knn
        T_results_rtpd = T_pwr_hourly((T_pwr_hourly.DATE.Month==this_month)&(T_pwr_hourly.DATE.Year==this_year), :); % Result container
        for i = 1: size(T_results_rtpd, 1)
            this_date = datetime(T_results_rtpd.HOUR_START.Year(i), T_results_rtpd.HOUR_START.Month(i), T_results_rtpd.HOUR_START.Day(i), 'TimeZone', 'UTC');
            this_hour = T_results_rtpd.HOUR_START.Hour(i);
%             T_sample = T_pwr_hourly((T_pwr_hourly.DATE<this_date)&(T_pwr_hourly.HOUR_START.Hour==this_hour), :);
            T_sample = T_pwr_hourly((T_pwr_hourly.DATE<this_date)&(T_pwr_hourly.HOUR_START.Hour==this_hour)&(T_pwr_hourly.HOUR_START>=this_date-days(180)), :);
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


%% Load CAISO data RTD and RTPD
load_caiso_data;

%% Load IBM data, forecast 15-min
load_ibm_5sites;

%% Explore different classifiers, one-dim, RTPD
% Select month and k
% this_year = 2019;
% this_month = 10;
karray = 5:5:60;

% Result container
cell_baseline_rtpd = cell(numel(karray), 5); % site, k, classifier
cell_results_rtpd  = cell(numel(karray), 5, 12); % site, k, classifier

for s = 1: 5
    fprintf('s = %g\n', s);
    
    T_pwr = cell_pwr{s}; % The 5th is CA_Topaz site
    T_pwr.TIME_START  = T_pwr.TIME - duration(0, dt_rtpd, 0);
    T_pwr.HOUR_START  = datetime(T_pwr.TIME_START.Year, T_pwr.TIME_START.Month, T_pwr.TIME_START.Day, T_pwr.TIME_START.Hour, 0, 0, 'TimeZone', 'UTC');

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
        switch classifier
            case 1
                % % Classifier 1: k (50 percentile), mean
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'k_p050'}), {'HOUR_START'}, 'mean'); 
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
                T_pwr_hourly.Properties.VariableNames{'mean_k_p050'} = 'classifier_1';
            case 2
                % % Classifier 2: k (50 percentile), std.
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'k_p050'}), {'HOUR_START'}, 'std'); 
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
                T_pwr_hourly.Properties.VariableNames{'std_k_p050'} = 'classifier_1';
            case 3
                % % Classifier 3: k (50 percentile), variability
                T_pwr.dk = [nan; diff(T_pwr.k_p050)]; % delta k
                T_pwr.dk_sq = [nan; diff(T_pwr.k_p050)].^2; % % (delta k)^2
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'dk_sq', 'k_width'}), {'HOUR_START'}, 'mean'); 
                T_pwr_hourly.v = sqrt(T_pwr_hourly.mean_dk_sq); % This is variability within each hour, following Inman et al. 2013, section 2.6.1
                T_pwr_hourly.Properties.VariableNames{'v'} = 'classifier_1';
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
            case 4
                % % Classifier 4: k_pv (50 percentile), mean
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'kpv_p050'}), {'HOUR_START'}, 'mean'); 
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
                T_pwr_hourly.Properties.VariableNames{'mean_kpv_p050'} = 'classifier_1';
            case 5
                % % Classifier 5: k_pv (50 percentile), std
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'kpv_p050'}), {'HOUR_START'}, 'std'); 
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
                T_pwr_hourly.Properties.VariableNames{'std_kpv_p050'} = 'classifier_1';
            case 6
                % % Classifier 6: k (50 percentile), variability
                T_pwr.dkpv = [nan; diff(T_pwr.kpv_p050)]; % delta k
                T_pwr.dkpv_sq = [nan; diff(T_pwr.kpv_p050)].^2; % % (delta k)^2
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'dkpv_sq'}), {'HOUR_START'}, 'mean'); 
                T_pwr_hourly.vpv = sqrt(T_pwr_hourly.mean_dkpv_sq); % This is variability within each hour, following Inman et al. 2013, section 2.6.1
                T_pwr_hourly.Properties.VariableNames{'vpv'} = 'classifier_1';
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
            case 7
                % % Classifier 7: width of k (75 - 25 percentile), mean
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'k_width'}), {'HOUR_START'}, 'mean'); 
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
                T_pwr_hourly.Properties.VariableNames{'mean_k_width'} = 'classifier_1';
            case 8
                % % Classifier 8: width of k (75 - 25 percentile), std
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'k_width'}), {'HOUR_START'}, 'std'); 
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
                T_pwr_hourly.Properties.VariableNames{'std_k_width'} = 'classifier_1';
            case 9
                % % Classifier 9: width of k (75 - 25 percentile), variability
                T_pwr.dw = [nan; diff(T_pwr.k_width)]; % delta k width
                T_pwr.dw_sq = [nan; diff(T_pwr.dw)].^2; % % (delta k)^2
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'dw_sq'}), {'HOUR_START'}, 'mean'); 
                T_pwr_hourly.vw = sqrt(T_pwr_hourly.mean_dw_sq); % This is variability within each hour, following Inman et al. 2013, section 2.6.1
                T_pwr_hourly.Properties.VariableNames{'vw'} = 'classifier_1';
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
            case 10
                % % Classifier 10: width of kpv (75 - 25 percentile), mean
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'kpv_width'}), {'HOUR_START'}, 'mean'); 
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
                T_pwr_hourly.Properties.VariableNames{'mean_kpv_width'} = 'classifier_1';
            case 11
                % % Classifier 11: width of k (75 - 25 percentile), std
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'kpv_width'}), {'HOUR_START'}, 'std'); 
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
                T_pwr_hourly.Properties.VariableNames{'std_kpv_width'} = 'classifier_1';
            case 12
                % % Classifier 12: width of k (75 - 25 percentile), variability
                T_pwr.dwpv = [nan; diff(T_pwr.kpv_width)]; % delta k width
                T_pwr.dwpv_sq = [nan; diff(T_pwr.dwpv)].^2; % % (delta k)^2
                T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'dwpv_sq'}), {'HOUR_START'}, 'mean'); 
                T_pwr_hourly.vwpv = sqrt(T_pwr_hourly.mean_dwpv_sq); % This is variability within each hour, following Inman et al. 2013, section 2.6.1
                T_pwr_hourly.Properties.VariableNames{'vwpv'} = 'classifier_1';
                T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');
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
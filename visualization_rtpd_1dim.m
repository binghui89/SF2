%% Visualization, RTPD

risk_factor = 1; % This is the coefficient for 1 MW of reserve shortage

% This is the actual need of FRP
f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
fru_need_rtpd = max(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4), [], 1)';
f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
frd_need_rtpd = min(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4), [], 1)';
T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});


for s = 2
    % Magnitude of imbalance, including over and under supply, hourly
    frp_rtpd_imbalance_baseline_hourly = zeros(numel(karray), 1, 24);
    fru_rtpd_over_baseline_hourly  = zeros(numel(karray), 1, 24);
    fru_rtpd_short_baseline_hourly = zeros(numel(karray), 1, 24);
    frd_rtpd_over_baseline_hourly  = zeros(numel(karray), 1, 24);
    frd_rtpd_short_baseline_hourly = zeros(numel(karray), 1, 24);
    
    frp_rtpd_imbalance_knn_hourly = zeros(numel(karray), numel(classifier), 24);
    fru_rtpd_over_knn_hourly  = zeros(numel(karray), numel(classifier), 24);
    fru_rtpd_short_knn_hourly = zeros(numel(karray), numel(classifier), 24);
    frd_rtpd_over_knn_hourly  = zeros(numel(karray), numel(classifier), 24);
    frd_rtpd_short_knn_hourly = zeros(numel(karray), numel(classifier), 24);
    
    % Frequency of shortage, use hourly data, hourly
    frp_rtpd_freqshort_baseline_hourly = zeros(numel(karray), 1, 24);
    fru_rtpd_freqshort_baseline_hourly = zeros(numel(karray), 1, 24);
    frd_rtpd_freqshort_baseline_hourly = zeros(numel(karray), 1, 24);
    
    frp_rtpd_freqshort_knn_hourly = zeros(numel(karray), numel(classifier), 24);
    fru_rtpd_freqshort_knn_hourly = zeros(numel(karray), numel(classifier), 24);
    frd_rtpd_freqshort_knn_hourly = zeros(numel(karray), numel(classifier), 24);
    
    % Frequency of shortage, use quarterly data, hourly
    frp_rtpd_freqshort_baseline_hourly_hd = zeros(numel(karray), 1, 24);
    fru_rtpd_freqshort_baseline_hourly_hd = zeros(numel(karray), 1, 24);
    frd_rtpd_freqshort_baseline_hourly_hd = zeros(numel(karray), 1, 24);
    
    frp_rtpd_freqshort_knn_hourly_hd = zeros(numel(karray), numel(classifier), 24);
    fru_rtpd_freqshort_knn_hourly_hd = zeros(numel(karray), numel(classifier), 24);
    frd_rtpd_freqshort_knn_hourly_hd = zeros(numel(karray), numel(classifier), 24);

    % Mean shortage, hourly, using 15-min data
    frp_rtpd_meanshort_baseline_hourly_hd = zeros(numel(karray), 1, 24);
    fru_rtpd_meanshort_baseline_hourly_hd = zeros(numel(karray), 1, 24);
    frd_rtpd_meanshort_baseline_hourly_hd = zeros(numel(karray), 1, 24);
    
    frp_rtpd_meanshort_knn_hourly_hd = zeros(numel(karray), numel(classifier), 24);
    fru_rtpd_meanshort_knn_hourly_hd = zeros(numel(karray), numel(classifier), 24);
    frd_rtpd_meanshort_knn_hourly_hd = zeros(numel(karray), numel(classifier), 24);

    % Mean over supply, hourly, using 15-min data
    frp_rtpd_meanover_baseline_hourly_hd = zeros(numel(karray), 1, 24);
    fru_rtpd_meanover_baseline_hourly_hd = zeros(numel(karray), 1, 24);
    frd_rtpd_meanover_baseline_hourly_hd = zeros(numel(karray), 1, 24);
    
    frp_rtpd_meanover_knn_hourly_hd = zeros(numel(karray), numel(classifier), 24);
    fru_rtpd_meanover_knn_hourly_hd = zeros(numel(karray), numel(classifier), 24);
    frd_rtpd_meanover_knn_hourly_hd = zeros(numel(karray), numel(classifier), 24);

    % Total imbalance, lumped sum of all time
    frp_rtpd_imbalance_baseline = zeros(numel(karray), 1);
    fru_rtpd_over_baseline  = zeros(numel(karray), 1);
    fru_rtpd_short_baseline = zeros(numel(karray), 1);
    frd_rtpd_over_baseline  = zeros(numel(karray), 1);
    frd_rtpd_short_baseline = zeros(numel(karray), 1);
    
    frp_rtpd_imbalance_knn = zeros(numel(karray), numel(classifier));
    fru_rtpd_over_knn  = zeros(numel(karray), numel(classifier));
    fru_rtpd_short_knn = zeros(numel(karray), numel(classifier));
    frd_rtpd_over_knn  = zeros(numel(karray), numel(classifier));
    frd_rtpd_short_knn = zeros(numel(karray), numel(classifier));
        
        % Total frequency of shortage, use quarterly data
        frp_rtpd_freqshort_baseline_hd = zeros(numel(karray), 1);
        fru_rtpd_freqshort_baseline_hd = zeros(numel(karray), 1);
        frd_rtpd_freqshort_baseline_hd = zeros(numel(karray), 1);
        
        frp_rtpd_freqshort_knn_hd = zeros(numel(karray), numel(classifier));
        fru_rtpd_freqshort_knn_hd = zeros(numel(karray), numel(classifier));
        frd_rtpd_freqshort_knn_hd = zeros(numel(karray), numel(classifier));
            
    for k = karray
        T_baseline_rtpd = cell_baseline_rtpd{karray==k, s};
        % Calculate baseline FRP imbalance
        T_baseline_rtpd.FRU_error = T_baseline_rtpd.FRU - fru_need_rtpd;
        T_baseline_rtpd.FRD_error = T_baseline_rtpd.FRD - frd_need_rtpd;

        T_baseline_rtpd.hour_start_local = datetime(T_baseline_rtpd.HOUR_START, 'TimeZone', '-08:00');
        T_fruerror_rtpd = array2table(T_baseline_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
        T_frderror_rtpd = array2table(T_baseline_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
        T_baseline_rtpd = [T_baseline_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
        T_baseline_rtpd.oclock_local = T_baseline_rtpd.hour_start_local.Hour;
        T_baseline_rtpd.FRU_short_freq_hd = sum(T_baseline_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
        T_baseline_rtpd.FRD_short_freq_hd = sum(T_baseline_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

        for ih = 1: 24
            h = ih - 1;
            fru_rtpd_over_baseline_hourly(karray==k, 1, ih)  = abs(sum(T_baseline_rtpd.FRU_error( (T_baseline_rtpd.FRU_error>=0) & (T_baseline_rtpd.oclock_local==h) )));
            fru_rtpd_short_baseline_hourly(karray==k, 1, ih) = abs(sum(T_baseline_rtpd.FRU_error( (T_baseline_rtpd.FRU_error<=0) & (T_baseline_rtpd.oclock_local==h) )));
            frd_rtpd_over_baseline_hourly(karray==k, 1, ih)  = abs(sum(T_baseline_rtpd.FRD_error( (T_baseline_rtpd.FRD_error<=0) & (T_baseline_rtpd.oclock_local==h) )));
            frd_rtpd_short_baseline_hourly(karray==k, 1, ih) = abs(sum(T_baseline_rtpd.FRD_error( (T_baseline_rtpd.FRD_error>=0) & (T_baseline_rtpd.oclock_local==h) )));

%             fru_rtpd_freqshort_baseline_hourly(karray==k, 1, ih) = sum((T_baseline_rtpd.FRU_error<=0) & (T_baseline_rtpd.oclock_local==h))/size(T_baseline_rtpd, 1);
%             frd_rtpd_freqshort_baseline_hourly(karray==k, 1, ih) = sum((T_baseline_rtpd.FRD_error>=0) & (T_baseline_rtpd.oclock_local==h))/size(T_baseline_rtpd, 1);
            fru_rtpd_freqshort_baseline_hourly(karray==k, 1, ih) = sum((T_baseline_rtpd.FRU_error<=0) & (T_baseline_rtpd.oclock_local==h))/sum(T_baseline_rtpd.oclock_local==h);
            frd_rtpd_freqshort_baseline_hourly(karray==k, 1, ih) = sum((T_baseline_rtpd.FRD_error>=0) & (T_baseline_rtpd.oclock_local==h))/sum(T_baseline_rtpd.oclock_local==h);
            
            fru_rtpd_freqshort_baseline_hourly_hd(karray==k, 1, ih) = mean(T_baseline_rtpd.FRU_short_freq_hd(T_baseline_rtpd.oclock_local==h));
            frd_rtpd_freqshort_baseline_hourly_hd(karray==k, 1, ih) = mean(T_baseline_rtpd.FRD_short_freq_hd(T_baseline_rtpd.oclock_local==h));
            
            tmp = T_baseline_rtpd{T_baseline_rtpd.oclock_local==h, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}};
            fru_rtpd_meanshort_baseline_hourly_hd(karray==k, 1, ih) = abs(sum(tmp(tmp<0))/sum(tmp(:)<0));
            tmp = T_baseline_rtpd{T_baseline_rtpd.oclock_local==h, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}};
            frd_rtpd_meanshort_baseline_hourly_hd(karray==k, 1, ih) = abs(sum(tmp(tmp>0))/sum(tmp(:)>0));
            
            tmp = T_baseline_rtpd{T_baseline_rtpd.oclock_local==h, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}};
            fru_rtpd_meanover_baseline_hourly_hd(karray==k, 1, ih) = abs(sum(tmp(tmp<0))/sum(tmp(:)<0));
            tmp = T_baseline_rtpd{T_baseline_rtpd.oclock_local==h, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}};
            frd_rtpd_meanover_baseline_hourly_hd(karray==k, 1, ih) = abs(sum(tmp(tmp>0))/sum(tmp(:)>0));
        end

        fru_rtpd_over_baseline(karray==k)  = abs(sum(T_baseline_rtpd.FRU_error(T_baseline_rtpd.FRU_error>=0)));
        fru_rtpd_short_baseline(karray==k) = abs(sum(T_baseline_rtpd.FRU_error(T_baseline_rtpd.FRU_error<=0)));
        frd_rtpd_over_baseline(karray==k)  = abs(sum(T_baseline_rtpd.FRD_error(T_baseline_rtpd.FRD_error<=0)));
        frd_rtpd_short_baseline(karray==k) = abs(sum(T_baseline_rtpd.FRD_error(T_baseline_rtpd.FRD_error>=0)));
        frp_rtpd_imbalance_baseline(karray==k) = fru_rtpd_over_baseline(karray==k) + risk_factor*fru_rtpd_short_baseline(karray==k) + frd_rtpd_over_baseline(karray==k) + risk_factor*frd_rtpd_short_baseline(karray==k);
        
        tmp = T_baseline_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
        fru_rtpd_freqshort_baseline_hd(karray==k) = sum(tmp(:))/numel(tmp);
        tmp = T_baseline_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
        frd_rtpd_freqshort_baseline_hd(karray==k) = sum(tmp(:))/numel(tmp);

        for classifier = 1: 12
            T_results_rtpd = cell_results_rtpd{karray==k, s, classifier};
            
            % FRU and FRD errors compared with actual needs
            T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
            T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

            T_results_rtpd.hour_start_local = datetime(T_results_rtpd.HOUR_START, 'TimeZone', '-08:00');
            T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
            T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
            T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
            T_results_rtpd.oclock_local = T_results_rtpd.hour_start_local.Hour;
            T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
            T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;
            
            for ih = 1: 24
                h = ih - 1;
                fru_rtpd_over_knn_hourly(karray==k, classifier, ih)  = abs(sum(T_results_rtpd.FRU_error( (T_results_rtpd.FRU_error>=0) & (T_results_rtpd.oclock_local==h) )));
                fru_rtpd_short_knn_hourly(karray==k, classifier, ih) = abs(sum(T_results_rtpd.FRU_error( (T_results_rtpd.FRU_error<=0) & (T_results_rtpd.oclock_local==h) )));
                frd_rtpd_over_knn_hourly(karray==k, classifier, ih)  = abs(sum(T_results_rtpd.FRD_error( (T_results_rtpd.FRD_error<=0) & (T_results_rtpd.oclock_local==h) )));
                frd_rtpd_short_knn_hourly(karray==k, classifier, ih) = abs(sum(T_results_rtpd.FRD_error( (T_results_rtpd.FRD_error>=0) & (T_results_rtpd.oclock_local==h) )));
                
%                 fru_rtpd_freqshort_knn_hourly(karray==k, classifier, ih) = sum((T_results_rtpd.FRU_error<=0) & (T_results_rtpd.oclock_local==h))/size(T_results_rtpd, 1);
%                 frd_rtpd_freqshort_knn_hourly(karray==k, classifier, ih) = sum((T_results_rtpd.FRD_error>=0) & (T_results_rtpd.oclock_local==h))/size(T_results_rtpd, 1);
                fru_rtpd_freqshort_knn_hourly(karray==k, classifier, ih) = sum((T_results_rtpd.FRU_error<=0) & (T_results_rtpd.oclock_local==h))/sum(T_results_rtpd.oclock_local==h);
                frd_rtpd_freqshort_knn_hourly(karray==k, classifier, ih) = sum((T_results_rtpd.FRD_error>=0) & (T_results_rtpd.oclock_local==h))/sum(T_results_rtpd.oclock_local==h);
                
                fru_rtpd_freqshort_knn_hourly_hd(karray==k, classifier, ih) = mean(T_results_rtpd.FRU_short_freq_hd(T_results_rtpd.oclock_local==h));
                frd_rtpd_freqshort_knn_hourly_hd(karray==k, classifier, ih) = mean(T_results_rtpd.FRD_short_freq_hd(T_results_rtpd.oclock_local==h));
                
                tmp = T_results_rtpd{T_results_rtpd.oclock_local==h, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}};
                fru_rtpd_meanshort_knn_hourly_hd(karray==k, classifier, ih) = abs(sum(tmp(tmp<0))/sum(tmp(:)<0));
                tmp = T_results_rtpd{T_results_rtpd.oclock_local==h, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}};
                frd_rtpd_meanshort_knn_hourly_hd(karray==k, classifier, ih) = abs(sum(tmp(tmp>0))/sum(tmp(:)>0));
                
                tmp = T_results_rtpd{T_results_rtpd.oclock_local==h, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}};
                fru_rtpd_meanover_knn_hourly_hd(karray==k, classifier, ih) = abs(sum(tmp(tmp<0))/sum(tmp(:)<0));
                tmp = T_results_rtpd{T_results_rtpd.oclock_local==h, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}};
                frd_rtpd_meanover_knn_hourly_hd(karray==k, classifier, ih) = abs(sum(tmp(tmp>0))/sum(tmp(:)>0));
            end

            fru_rtpd_over_knn(karray==k, classifier)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
            fru_rtpd_short_knn(karray==k, classifier) = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error<=0)));
            frd_rtpd_over_knn(karray==k, classifier)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
            frd_rtpd_short_knn(karray==k, classifier) = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error>=0)));
            frp_rtpd_imbalance_knn(karray==k, classifier) = fru_rtpd_over_knn(karray==k, classifier) + risk_factor*fru_rtpd_short_knn(karray==k, classifier) + frd_rtpd_over_knn(karray==k, classifier) + risk_factor*frd_rtpd_short_knn(karray==k, classifier);

            tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
            fru_rtpd_freqshort_knn_hd(karray==k, classifier) = sum(tmp(:))/numel(tmp);
            tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
            frd_rtpd_freqshort_knn_hd(karray==k, classifier) = sum(tmp(:))/numel(tmp);

        end
        frp_rtpd_imbalance_baseline_hourly = fru_rtpd_over_baseline_hourly + risk_factor.*fru_rtpd_short_baseline_hourly + frd_rtpd_over_baseline_hourly + risk_factor.*frd_rtpd_short_baseline_hourly;
        frp_rtpd_imbalance_knn_hourly = fru_rtpd_over_knn_hourly + risk_factor.*fru_rtpd_short_knn_hourly + frd_rtpd_over_knn_hourly + risk_factor.*frd_rtpd_short_knn_hourly;
        frp_rtpd_freqshort_baseline_hourly = fru_rtpd_freqshort_baseline_hourly + frd_rtpd_freqshort_baseline_hourly;
        frp_rtpd_freqshort_baseline_hourly_hd = fru_rtpd_freqshort_baseline_hourly_hd + frd_rtpd_freqshort_baseline_hourly_hd;
        frp_rtpd_freqshort_knn_hourly = fru_rtpd_freqshort_knn_hourly + frd_rtpd_freqshort_knn_hourly;
        frp_rtpd_freqshort_knn_hourly_hd = fru_rtpd_freqshort_knn_hourly_hd + frd_rtpd_freqshort_knn_hourly_hd;
        frp_rtpd_meanshort_baseline_hourly_hd = fru_rtpd_meanshort_baseline_hourly_hd + frd_rtpd_meanshort_baseline_hourly_hd;
        frp_rtpd_meanshort_knn_hourly_hd = fru_rtpd_meanshort_knn_hourly_hd + frd_rtpd_meanshort_knn_hourly_hd;
        frp_rtpd_meanover_baseline_hourly_hd = fru_rtpd_meanshort_baseline_hourly_hd + frd_rtpd_meanshort_baseline_hourly_hd;
        frp_rtpd_meanover_knn_hourly_hd = fru_rtpd_meanshort_knn_hourly_hd + frd_rtpd_meanshort_knn_hourly_hd;
        frp_rtpd_freqshort_baseline_hd = fru_rtpd_freqshort_baseline_hd + frd_rtpd_freqshort_baseline_hd;
        frp_rtpd_freqshort_knn_hd = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
    end

    %% Risk-adjusted imbalance
    figure();
    h = plot(karray, frp_rtpd_imbalance_knn, karray, frp_rtpd_imbalance_baseline, '-k.');
    set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
    set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
    set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
    set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
    set(h([2, 5, 8, 11]), 'LineStyle', '--');
    set(h([3, 6, 9, 12]), 'LineStyle', ':');
    set(h(1:12), 'Marker', '.');
    title(strcat('Risk-adjusted imbalance: ', uniquegensite(s)));
%     legend({'\mu(k)', '\sigma(k)', 'u(k)', '\mu(k_P)', '\sigma(k_P)', 'u(k_P)', '\mu(w)', '\sigma(w)', 'u(w)', '\mu(w_P)', '\sigma(w_P)', 'u(w_P)', 'Baseline'});
    ylabel('MW'); xlabel('k');

    %% Over supply, shortage, 4-in-1 plot
    figure();
    subplot(3, 2, 1);
    h = plot(karray, fru_rtpd_over_knn, karray, fru_rtpd_over_baseline, '-k.');
    set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
    set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
    set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
    set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
    set(h([2, 5, 8, 11]), 'LineStyle', '--');
    set(h([3, 6, 9, 12]), 'LineStyle', ':');
    set(h(1:12), 'Marker', '.');
    title('FRU over supply');
    ylabel('MW'); xlabel('k');

    subplot(3, 2, 2);
    h = plot(karray, frd_rtpd_over_knn, karray, frd_rtpd_over_baseline, '-k.'); 
    set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
    set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
    set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
    set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
    set(h([2, 5, 8, 11]), 'LineStyle', '--');
    set(h([3, 6, 9, 12]), 'LineStyle', ':');
    set(h(1:12), 'Marker', '.');
    title('FRD over supply');
    ylabel('MW'); xlabel('k');

    subplot(3, 2, 3);
    h = plot(karray, fru_rtpd_short_knn, karray, fru_rtpd_short_baseline, '-k.'); 
    set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
    set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
    set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
    set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
    set(h([2, 5, 8, 11]), 'LineStyle', '--');
    set(h([3, 6, 9, 12]), 'LineStyle', ':');
    set(h(1:12), 'Marker', '.');
    title('FRU shortage');
    ylabel('MW'); xlabel('k');

    subplot(3, 2, 4);
    h = plot(karray, frd_rtpd_short_knn, karray, frd_rtpd_short_baseline, '-k.'); 
    set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
    set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
    set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
    set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
    set(h([2, 5, 8, 11]), 'LineStyle', '--');
    set(h([3, 6, 9, 12]), 'LineStyle', ':');
    set(h(1:12), 'Marker', '.');
    title('FRD shortage');
    ylabel('MW'); xlabel('k');

    subplot(3, 2, 5);
    h = plot(karray, fru_rtpd_freqshort_knn_hd, karray, fru_rtpd_freqshort_baseline_hd, '-k.'); 
    set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
    set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
    set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
    set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
    set(h([2, 5, 8, 11]), 'LineStyle', '--');
    set(h([3, 6, 9, 12]), 'LineStyle', ':');
    set(h(1:12), 'Marker', '.');
    title('FRU shortage frequency');
    ylabel('MW'); xlabel('k');

    subplot(3, 2, 6);
    h = plot(karray, frd_rtpd_freqshort_knn_hd, karray, frd_rtpd_freqshort_baseline_hd, '-k.'); 
    set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
    set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
    set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
    set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
    set(h([2, 5, 8, 11]), 'LineStyle', '--');
    set(h([3, 6, 9, 12]), 'LineStyle', ':');
    set(h(1:12), 'Marker', '.');
    title('FRD shortage frequency');
    ylabel('MW'); xlabel('k');
    sgtitle(strcat('Risk-adjusted imbalance: ', uniquegensite(s)));
    
    %% Risk-adjusted imbalance by hour
    figure();
    h_plot = 9:16;
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        subplot(3, 3, ih);
        h = plot(karray, frp_rtpd_imbalance_knn_hourly(:, :, hr), karray, squeeze(frp_rtpd_imbalance_baseline_hourly(:, :, hr)), '-k.');
        set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
        set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
        set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
        set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
        set(h([2, 5, 8, 11]), 'LineStyle', '--');
        set(h([3, 6, 9, 12]), 'LineStyle', ':');
        set(h(1:12), 'Marker', '.');
        xlim([min(karray), max(karray)]);
        title(strcat('Hour ', num2str(hr)));
    end
    sgtitle(strcat('Risk adjusted imbalance, site:', uniquegensite(s)));
    
    %% Over and under supply of FRU and FRD by hour
    figure();
    h_plot = 9:16;
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        subplot(3, 3, ih);
        h = plot(karray, fru_rtpd_over_knn_hourly(:, :, hr), karray, squeeze(fru_rtpd_over_baseline_hourly(:, :, hr)), '-k.');
        set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
        set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
        set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
        set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
        set(h([2, 5, 8, 11]), 'LineStyle', '--');
        set(h([3, 6, 9, 12]), 'LineStyle', ':');
        set(h(1:12), 'Marker', '.');
        xlim([min(karray), max(karray)]);
        title(strcat('Hour ', num2str(hr)));
    end
    sgtitle(strcat('FRU over procurement, site:', uniquegensite(s)));
    
    figure();
    h_plot = 9:16;
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        subplot(3, 3, ih);
        h = plot(karray, fru_rtpd_short_knn_hourly(:, :, hr), karray, squeeze(fru_rtpd_short_baseline_hourly(:, :, hr)), '-k.');
        set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
        set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
        set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
        set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
        set(h([2, 5, 8, 11]), 'LineStyle', '--');
        set(h([3, 6, 9, 12]), 'LineStyle', ':');
        set(h(1:12), 'Marker', '.');
        xlim([min(karray), max(karray)]);
        title(strcat('Hour ', num2str(hr)));
    end
    sgtitle(strcat('FRU shortage, site:', uniquegensite(s)));
    
    figure();
    h_plot = 9:16;
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        subplot(3, 3, ih);
        h = plot(karray, frd_rtpd_over_knn_hourly(:, :, hr), karray, squeeze(frd_rtpd_over_baseline_hourly(:, :, hr)), '-k.');
        set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
        set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
        set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
        set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
        set(h([2, 5, 8, 11]), 'LineStyle', '--');
        set(h([3, 6, 9, 12]), 'LineStyle', ':');
        set(h(1:12), 'Marker', '.');
        xlim([min(karray), max(karray)]);
        title(strcat('Hour ', num2str(hr)));
    end
    sgtitle(strcat('FRD over procurement, site:', uniquegensite(s)));
    
    figure();
    h_plot = 9:16;
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        subplot(3, 3, ih);
        h = plot(karray, frd_rtpd_short_knn_hourly(:, :, hr), karray, squeeze(frd_rtpd_short_baseline_hourly(:, :, hr)), '-k.');
        set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
        set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
        set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
        set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
        set(h([2, 5, 8, 11]), 'LineStyle', '--');
        set(h([3, 6, 9, 12]), 'LineStyle', ':');
        set(h(1:12), 'Marker', '.');
        xlim([min(karray), max(karray)]);
        title(strcat('Hour ', num2str(hr)));
    end
    sgtitle(strcat('FRD shortage, site:', uniquegensite(s)));
    
    %% Frequency of shortage by hour
    figure();
    h_plot = 9:16;
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        subplot(3, 3, ih);
        h = plot(karray, frp_rtpd_freqshort_knn_hourly(:, :, hr), karray, squeeze(frp_rtpd_freqshort_baseline_hourly(:, :, hr)), '-k.');
        set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
        set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
        set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
        set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
        set(h([2, 5, 8, 11]), 'LineStyle', '--');
        set(h([3, 6, 9, 12]), 'LineStyle', ':');
        set(h(1:12), 'Marker', '.');
        xlim([min(karray), max(karray)]);
        title(strcat('Hour ', num2str(hr)));
    end
    sgtitle(strcat('FRP LOLP, site:', uniquegensite(s)));
    
    figure();
    h_plot = 9:16;
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        subplot(3, 3, ih);
        h = plot(karray, fru_rtpd_freqshort_knn_hourly(:, :, hr), karray, squeeze(fru_rtpd_freqshort_baseline_hourly(:, :, hr)), '-k.');
        set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
        set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
        set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
        set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
        set(h([2, 5, 8, 11]), 'LineStyle', '--');
        set(h([3, 6, 9, 12]), 'LineStyle', ':');
        set(h(1:12), 'Marker', '.');
        xlim([min(karray), max(karray)]);
        title(strcat('Hour ', num2str(hr)));
    end
    sgtitle(strcat('FRU LOLP, site:', uniquegensite(s)));
    
    figure();
    h_plot = 9:16;
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        subplot(3, 3, ih);
        h = plot(karray, frd_rtpd_freqshort_knn_hourly(:, :, hr), karray, squeeze(frd_rtpd_freqshort_baseline_hourly(:, :, hr)), '-k.');
        set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
        set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
        set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
        set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
        set(h([2, 5, 8, 11]), 'LineStyle', '--');
        set(h([3, 6, 9, 12]), 'LineStyle', ':');
        set(h(1:12), 'Marker', '.');
        xlim([min(karray), max(karray)]);
        title(strcat('Hour ', num2str(hr)));
    end
    sgtitle(strcat('FRD LOLP, site:', uniquegensite(s)));
    
    %% HD frequency of shortage by hour (quarterly data used)
    figure();
    h_plot = 7:18;
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        subplot(4, 3, ih);
        h = plot(karray, frp_rtpd_freqshort_knn_hourly_hd(:, :, hr), karray, squeeze(frp_rtpd_freqshort_baseline_hourly_hd(:, :, hr)), '-k.');
        set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
        set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
        set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
        set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
        set(h([2, 5, 8, 11]), 'LineStyle', '--');
        set(h([3, 6, 9, 12]), 'LineStyle', ':');
        set(h(1:12), 'Marker', '.');
        xlim([min(karray), max(karray)]);
        ylim([0, 0.3]);
        title(strcat('Hour ', num2str(hr)));
    end
    sgtitle(strcat('FRP LOLP HD, site:', uniquegensite(s)));
    
    figure();
    h_plot = 7:18;
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        subplot(4, 3, ih);
        h = plot(karray, fru_rtpd_freqshort_knn_hourly_hd(:, :, hr), karray, squeeze(fru_rtpd_freqshort_baseline_hourly_hd(:, :, hr)), '-k.');
        set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
        set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
        set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
        set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
        set(h([2, 5, 8, 11]), 'LineStyle', '--');
        set(h([3, 6, 9, 12]), 'LineStyle', ':');
        set(h(1:12), 'Marker', '.');
        xlim([min(karray), max(karray)]);
        ylim([0, 0.3]);
        title(strcat('Hour ', num2str(hr)));
    end
    sgtitle(strcat('FRU LOLP HD, site:', uniquegensite(s)));
    
    figure();
    h_plot = 7:18;
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        subplot(4, 3, ih);
        h = plot(karray, frd_rtpd_freqshort_knn_hourly_hd(:, :, hr), karray, squeeze(frd_rtpd_freqshort_baseline_hourly_hd(:, :, hr)), '-k.');
        set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
        set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
        set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
        set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
        set(h([2, 5, 8, 11]), 'LineStyle', '--');
        set(h([3, 6, 9, 12]), 'LineStyle', ':');
        set(h(1:12), 'Marker', '.');
        xlim([min(karray), max(karray)]);
        ylim([0, 0.3]);
        title(strcat('Hour ', num2str(hr)));
    end
    sgtitle(strcat('FRD LOLP HD, site:', uniquegensite(s)));
    %% Pareto frontier
    figure();
    frp_rtpd_over_baseline = fru_rtpd_over_baseline + frd_rtpd_over_baseline;
    frp_rtpd_over_knn = fru_rtpd_over_knn + frd_rtpd_over_knn;
    hold on;
    h = nan(13, 1);
    for classifier = 1: 12
        h(classifier) = scatter(frp_rtpd_freqshort_knn_hd(:, classifier), frp_rtpd_over_knn(:, classifier)./1E3, 60);
    end
    h(end) = scatter(frp_rtpd_freqshort_baseline_hd, frp_rtpd_over_baseline./1E3, 60, 'ko');
    set(h(end), 'MarkerFaceColor', 'k'); 

    set(h(1:3), 'MarkerEdgeColor', [0, 0.4470, 0.7410]); % Blue
    set(h(4:6), 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980]); % Blue
    set(h(7:9), 'MarkerEdgeColor', [0.9290, 0.6940, 0.1250]); % Red
    set(h(10:12), 'MarkerEdgeColor', [0.4940, 0.1840, 0.5560]); % Red
    set(h(1:3), 'MarkerFaceColor', [0, 0.4470, 0.7410]); % Blue
    set(h(4:6), 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]); % Blue
    set(h(7:9), 'MarkerFaceColor', [0.9290, 0.6940, 0.1250]); % Red
    set(h(10:12), 'MarkerFaceColor', [0.4940, 0.1840, 0.5560]); % Red
    set(h([1, 4, 7, 10]), 'Marker', 's');
    set(h([2, 5, 8, 11]), 'Marker', 'd');
    set(h([3, 6, 9, 12]), 'Marker', 'h');
    set(h(:), 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
%     title(strcat('Pareto frontier: ', uniquegensite(s)));
    ylabel('FRP over supply (GW)'); xlabel('Probability of FRP shortage');
    xline(frp_rtpd_freqshort_baseline_hd(karray==30), 'Color','red', 'LineStyle','--');
    yline(frp_rtpd_over_baseline(karray==30)./1E3, 'Color','red', 'LineStyle','--');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    box on;
    legend(h, {'\mu(k)', '\sigma(k)', 'u(k)', '\mu(k_P)', '\sigma(k_P)', 'v(k_P)', '\mu(w)', '\sigma(w)', 'v(w)', '\mu(w_P)', '\sigma(w_P)', 'v(w_P)', 'Baseline'},'Location', 'eastoutside');

    fig = figure();
    subplot(2, 1, 1);
    h = plot(karray, frp_rtpd_freqshort_knn_hd, karray, frp_rtpd_freqshort_baseline_hd, '-k.'); 
    set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
    set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
    set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
    set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
    set(h([2, 5, 8, 11]), 'LineStyle', '--');
    set(h([3, 6, 9, 12]), 'LineStyle', ':');
    set(h, 'LineWidth', 1);
    set(h(1:12), 'Marker', '.');
%     xlabel('k');
    ylabel('FRP shortage');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    legend({'\mu(k)', '\sigma(k)', 'u(k)', '\mu(k_P)', '\sigma(k_P)', 'v(k_P)', '\mu(w)', '\sigma(w)', 'v(w)', '\mu(w_P)', '\sigma(w_P)', 'v(w_P)', 'Baseline'}, 'Location', 'eastoutside');
    
    set(gca,'xtick',[]);
    
%     yyaxis right;
%     figure();
    subplot(2, 1, 2);
    h = plot(karray, frp_rtpd_over_knn./1E3, karray, frp_rtpd_over_baseline./1E3, '-k.'); 
    set(h(1:3), 'Color', [0, 0.4470, 0.7410]);
    set(h(4:6), 'Color', [0.8500, 0.3250, 0.0980]);
    set(h(7:9), 'Color', [0.9290, 0.6940, 0.1250]);
    set(h(10:12), 'Color', [0.4940, 0.1840, 0.5560]);
    set(h([2, 5, 8, 11]), 'LineStyle', '--');
    set(h([3, 6, 9, 12]), 'LineStyle', ':');
    set(h, 'LineWidth', 1);
    set(h(1:12), 'Marker', '.');
    xlabel('k');
    ylabel('Over supply (GW)');
%     title(uniquegensite(s));
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
%     legend({'\mu(k)', '\sigma(k)', 'u(k)', '\mu(k_P)', '\sigma(k_P)', 'v(k_P)', '\mu(w)', '\sigma(w)', 'v(w)', '\mu(w_P)', '\sigma(w_P)', 'v(w_P)', 'Baseline'}, 'Location', 'eastoutside');

%     ha=get(fig,'children');
%     box3 = ha(3).Position;
%     box1 = ha(1).Position;
%     set(ha(1), 'Position', [box1(1), box1(2), box3(3), box1(4)]);
    %% Paper figures, only use k related classifiers
    
    % FRD risk and over procurement
    h_plot = [7, 10, 12, 15, 18];
    
    % FRU risk and over procurement
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        figure();
        yyaxis left;
        h = plot(karray, fru_rtpd_freqshort_knn_hourly_hd(:, 1:3, hr), karray, squeeze(fru_rtpd_freqshort_baseline_hourly_hd(:, :, hr)));
        set(h(1), 'LineStyle', '--', 'LineWidth', 2);
        set(h(2), 'LineStyle', '-.', 'LineWidth', 2);
        set(h(3), 'LineStyle', ':', 'LineWidth', 2);
        set(h(4), 'LineStyle', '-', 'LineWidth', 2);
        xlim([min(karray), max(karray)]);
        ylim([0, 0.3]);
        xlabel('k');
        ylabel('Risk of shortage');
        
        yyaxis right;
        h = plot(karray, fru_rtpd_over_knn_hourly(:, 1:3, hr)./1E3, karray, squeeze(fru_rtpd_over_baseline_hourly(:, :, hr))./1E3);
        set(h(1), 'LineStyle', '--', 'LineWidth', 2);
        set(h(2), 'LineStyle', '-.', 'LineWidth', 2);
        set(h(3), 'LineStyle', ':', 'LineWidth', 2);
        set(h(4), 'LineStyle', '-', 'LineWidth', 2);
        ylabel('Over supply (GW)');
%         title(strcat('Hour ', num2str(hr)));
        set(findall(gcf,'-property','FontSize'),'FontSize',22);
        ylim([0, 40]);
        hold on;
        h = plot(nan, nan, '-k', nan, nan, '--k', nan, nan, '-.k', nan, nan, ':k');
        legend(h, {'Baseline', '\mu(k)', '\sigma(k)', 'v(k)'}, 'NumColumns', 2);
    end
    
    % FRD risk and over procurement
    for ih = 1:numel(h_plot)
        hr = h_plot(ih);
        figure();
        yyaxis left;
        h = plot(karray, frd_rtpd_freqshort_knn_hourly_hd(:, 1:3, hr), karray, squeeze(frd_rtpd_freqshort_baseline_hourly_hd(:, :, hr)));
        set(h(1), 'LineStyle', '--', 'LineWidth', 2);
        set(h(2), 'LineStyle', '-.', 'LineWidth', 2);
        set(h(3), 'LineStyle', ':', 'LineWidth', 2);
        set(h(4), 'LineStyle', '-', 'LineWidth', 2);
        xlim([min(karray), max(karray)]);
        ylim([0, 0.3]);
        xlabel('k');
        ylabel('Risk of shortage');
        
        yyaxis right;
        h = plot(karray, frd_rtpd_over_knn_hourly(:, 1:3, hr)./1E3, karray, squeeze(frd_rtpd_over_baseline_hourly(:, :, hr))./1E3);
        set(h(1), 'LineStyle', '--', 'LineWidth', 2);
        set(h(2), 'LineStyle', '-.', 'LineWidth', 2);
        set(h(3), 'LineStyle', ':', 'LineWidth', 2);
        set(h(4), 'LineStyle', '-', 'LineWidth', 2);
        ylabel('Over supply (GW)');
%         title(strcat('Hour ', num2str(hr)));
        set(findall(gcf,'-property','FontSize'),'FontSize',22);
        ylim([0, 40]);
        hold on;
        h = plot(nan, nan, '-k', nan, nan, '--k', nan, nan, '-.k', nan, nan, ':k');
        legend(h, {'Baseline', '\mu(k)', '\sigma(k)', 'v(k)'}, 'NumColumns', 2);
    end

    

end

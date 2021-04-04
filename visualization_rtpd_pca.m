% Figures for the invited paper (IEEE trans. sustainable energy)
% Based on the first submission (decision 3/14/21), add PCA

load('knn_rtpd_puresolar_2.mat', 'T_rtpd', 'T_baseline_rtpd', 'karray', 'cell_baseline_rtpd', 'cell_results_rtpd');

% This is the actual need of FRP
f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
fru_need_rtpd = max(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4), [], 1)';
f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
frd_need_rtpd = min(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4), [], 1)';
T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});

s = 2;

% Calculate metrics
for k = karray
    T_baseline_rtpd = cell_baseline_rtpd{karray==k, s};
    % Calculate baseline FRP imbalance
    T_baseline_rtpd.FRU_error = T_baseline_rtpd.FRU - fru_need_rtpd;
    T_baseline_rtpd.FRD_error = T_baseline_rtpd.FRD - frd_need_rtpd;

    T_fruerror_rtpd = array2table(T_baseline_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
    T_frderror_rtpd = array2table(T_baseline_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
    T_baseline_rtpd = [T_baseline_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
    T_baseline_rtpd.FRU_short_freq_hd = sum(T_baseline_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
    T_baseline_rtpd.FRD_short_freq_hd = sum(T_baseline_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

    fru_rtpd_over_baseline(karray==k)  = abs(sum(T_baseline_rtpd.FRU_error(T_baseline_rtpd.FRU_error>=0)));
    frd_rtpd_over_baseline(karray==k)  = abs(sum(T_baseline_rtpd.FRD_error(T_baseline_rtpd.FRD_error<=0)));

    tmp = T_baseline_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
    fru_rtpd_freqshort_baseline_hd(karray==k) = sum(tmp(:))/numel(tmp);
    tmp = T_baseline_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
    frd_rtpd_freqshort_baseline_hd(karray==k) = sum(tmp(:))/numel(tmp);

    for classifier = 1: size(cell_results_rtpd, 3)
        T_results_rtpd = cell_results_rtpd{karray==k, s, classifier};
        T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
        T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

        T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
        T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
        T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
        T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
        T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;
        tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
        fru_rtpd_freqshort_knn_hd(karray==k, classifier) = sum(tmp(:))/numel(tmp);
        tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
        frd_rtpd_freqshort_knn_hd(karray==k, classifier) = sum(tmp(:))/numel(tmp);

        fru_rtpd_over_knn(karray==k, classifier)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
        frd_rtpd_over_knn(karray==k, classifier)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
    end
end
frp_rtpd_freqshort_knn_hd = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
frp_rtpd_freqshort_baseline_hd = fru_rtpd_freqshort_baseline_hd + frd_rtpd_freqshort_baseline_hd;
frp_rtpd_over_knn = fru_rtpd_over_knn + frd_rtpd_over_knn;
frp_rtpd_over_baseline = fru_rtpd_over_baseline + frd_rtpd_over_baseline;

T_metric_over   = table(karray(:), frp_rtpd_over_baseline(:), frp_rtpd_over_knn, 'VariableNames', {'k', 'baseline', 'knn1d'});
T_metric_pshort = table(karray(:), frp_rtpd_freqshort_baseline_hd(:), frp_rtpd_freqshort_knn_hd, 'VariableNames', {'k', 'baseline', 'knn1d'});

%% Pareto frontier

figure();
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
ylabel('FRP over supply (GW)'); xlabel('Probability of FRP shortage');
xline(frp_rtpd_freqshort_baseline_hd(karray==30), 'Color','red', 'LineStyle','--');
yline(frp_rtpd_over_baseline(karray==30)./1E3, 'Color','red', 'LineStyle','--');
set(findall(gcf,'-property','FontSize'),'FontSize',14);
box on;
clearvars -except 'T_metric_over' 'T_metric_pshort' 'h';

%% Load single-site PCA+kNN results
load('knn_rtpd_puresolar_2.pca.mat', 'T_rtpd', 'T_baseline_rtpd', 'karray', 'cell_baseline_rtpd', 'cell_results_rtpd');

% This is the actual need of FRP
f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
fru_need_rtpd = max(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4), [], 1)';
f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
frd_need_rtpd = min(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4), [], 1)';
T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});

s = 2;

% Calculate metrics
for k = karray
    T_baseline_rtpd = cell_baseline_rtpd{karray==k, s};
    % Calculate baseline FRP imbalance
    T_baseline_rtpd.FRU_error = T_baseline_rtpd.FRU - fru_need_rtpd;
    T_baseline_rtpd.FRD_error = T_baseline_rtpd.FRD - frd_need_rtpd;

    T_fruerror_rtpd = array2table(T_baseline_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
    T_frderror_rtpd = array2table(T_baseline_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
    T_baseline_rtpd = [T_baseline_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
    T_baseline_rtpd.FRU_short_freq_hd = sum(T_baseline_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
    T_baseline_rtpd.FRD_short_freq_hd = sum(T_baseline_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

    fru_rtpd_over_baseline(karray==k)  = abs(sum(T_baseline_rtpd.FRU_error(T_baseline_rtpd.FRU_error>=0)));
    frd_rtpd_over_baseline(karray==k)  = abs(sum(T_baseline_rtpd.FRD_error(T_baseline_rtpd.FRD_error<=0)));

    tmp = T_baseline_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
    fru_rtpd_freqshort_baseline_hd(karray==k) = sum(tmp(:))/numel(tmp);
    tmp = T_baseline_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
    frd_rtpd_freqshort_baseline_hd(karray==k) = sum(tmp(:))/numel(tmp);

    for classifier = 1: size(cell_results_rtpd, 3)
        T_results_rtpd = cell_results_rtpd{karray==k, s, classifier};
        T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
        T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

        T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
        T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
        T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
        T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
        T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;
        tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
        fru_rtpd_freqshort_knn_hd(karray==k, classifier) = sum(tmp(:))/numel(tmp);
        tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
        frd_rtpd_freqshort_knn_hd(karray==k, classifier) = sum(tmp(:))/numel(tmp);

        fru_rtpd_over_knn(karray==k, classifier)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
        frd_rtpd_over_knn(karray==k, classifier)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
    end
end
frp_rtpd_freqshort_knn_hd = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
frp_rtpd_freqshort_baseline_hd = fru_rtpd_freqshort_baseline_hd + frd_rtpd_freqshort_baseline_hd;
frp_rtpd_over_knn = fru_rtpd_over_knn + frd_rtpd_over_knn;
frp_rtpd_over_baseline = fru_rtpd_over_baseline + frd_rtpd_over_baseline;
for classifier = 1: size(cell_results_rtpd, 3)
    h(numel(h) + 1) = plot(frp_rtpd_freqshort_knn_hd(:, classifier), frp_rtpd_over_knn(:, classifier)./1E3, '^');
end
set(h(numel(h)-2), 'Color', 'r', 'MarkerFaceColor', 'r', 'Marker', '^', 'MarkerSize', 6);
set(h(numel(h)-1), 'Color', 'g', 'MarkerFaceColor', 'g', 'Marker', '^', 'MarkerSize', 6);
set(h(numel(h)), 'Color', 'm', 'MarkerFaceColor', 'm', 'Marker', '^', 'MarkerSize', 6);

T_metric_over.pcaknn = frp_rtpd_over_knn;
T_metric_pshort.pcaknn = frp_rtpd_freqshort_knn_hd;

clearvars -except 'T_metric_over' 'T_metric_pshort' 'h';

%% Load 1-D kNN optimization results
load('knn_post_rtpd_puresolar_2_complete.mat', 'T_rtpd', 'T_baseline_rtpd', 'cell_optimalknn_rtpd', 'cell_results_rtpd', 'array_ndays');

s = 2;
fprintf('s=%g\n', s);
for i = 1 : size(cell_optimalknn_rtpd, 1)
    T_optimal_rtpd = cell_optimalknn_rtpd{i, s};
    f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_optimal_rtpd.HOUR_START), 'error_max'}; % 15-min
    f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_optimal_rtpd.HOUR_START), 'error_min'}; % 15-min
    T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
    T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});
    T_optimal_rtpd = [T_optimal_rtpd T_errormax_rtpd T_errormin_rtpd];
    T_optimal_rtpd.FRU_NEED = max(T_optimal_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, [], 2);
    T_optimal_rtpd.FRD_NEED = min(T_optimal_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, [], 2);

    % Calculate FRU/FRD errors for that hour and each 15-min interval
    T_fruerror_rtpd = array2table(T_optimal_rtpd.FRU-T_optimal_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
    T_frderror_rtpd = array2table(T_optimal_rtpd.FRD-T_optimal_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
    T_optimal_rtpd = [T_optimal_rtpd T_fruerror_rtpd T_frderror_rtpd];
    T_optimal_rtpd.FRU_error = T_optimal_rtpd.FRU - T_optimal_rtpd.FRU_NEED;
    T_optimal_rtpd.FRD_error = T_optimal_rtpd.FRD - T_optimal_rtpd.FRD_NEED;

    % Calculate evaluation metrics: Total oversupply
    fru_rtpd_over_knn = abs(sum(T_optimal_rtpd.FRU_error(T_optimal_rtpd.FRU_error>=0)));
    frd_rtpd_over_knn = abs(sum(T_optimal_rtpd.FRD_error(T_optimal_rtpd.FRD_error<=0)));
    frp_rtpd_over_knn(i) = fru_rtpd_over_knn + frd_rtpd_over_knn;

    tmp = T_optimal_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
    fru_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
    tmp = T_optimal_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
    frd_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
    frp_rtpd_freqshort_knn_hd(i) = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
    fprintf('n=%g, reliability=%.2f, over=%.2f\n', 5*i, frp_rtpd_freqshort_knn_hd(i), frp_rtpd_over_knn(i));
end
T_metric_over_opt   = table(array_ndays(:), frp_rtpd_over_knn(:), 'VariableNames', {'n', 'knnopt'});
T_metric_pshort_opt = table(array_ndays(:), frp_rtpd_freqshort_knn_hd(:), 'VariableNames', {'n', 'knnopt'});

h(numel(h)+1) = plot(frp_rtpd_freqshort_knn_hd, frp_rtpd_over_knn./1E3, '-^');
clearvars -except 'T_metric_over' 'T_metric_pshort' 'h' 'T_metric_over_opt' 'T_metric_pshort_opt';

%% Load multi-site n-D PCA kNN optimization results
load('knn_post_rtpd_puresolar_2_complete.pca.mat', 'T_rtpd', 'T_baseline_rtpd', 'cell_optimalknn_rtpd', 'cell_results_rtpd', 'array_ndays');

s = 1; % Note in multi-site run, s always = 1
fprintf('s=%g\n', s);
for i = 1 : size(cell_optimalknn_rtpd, 1)
    T_optimal_rtpd = cell_optimalknn_rtpd{i, s};
    f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_optimal_rtpd.HOUR_START), 'error_max'}; % 15-min
    f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_optimal_rtpd.HOUR_START), 'error_min'}; % 15-min
    T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
    T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});
    T_optimal_rtpd = [T_optimal_rtpd T_errormax_rtpd T_errormin_rtpd];
    T_optimal_rtpd.FRU_NEED = max(T_optimal_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, [], 2);
    T_optimal_rtpd.FRD_NEED = min(T_optimal_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, [], 2);

    % Calculate FRU/FRD errors for that hour and each 15-min interval
    T_fruerror_rtpd = array2table(T_optimal_rtpd.FRU-T_optimal_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
    T_frderror_rtpd = array2table(T_optimal_rtpd.FRD-T_optimal_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
    T_optimal_rtpd = [T_optimal_rtpd T_fruerror_rtpd T_frderror_rtpd];
    T_optimal_rtpd.FRU_error = T_optimal_rtpd.FRU - T_optimal_rtpd.FRU_NEED;
    T_optimal_rtpd.FRD_error = T_optimal_rtpd.FRD - T_optimal_rtpd.FRD_NEED;

    % Calculate evaluation metrics: Total oversupply
    fru_rtpd_over_knn = abs(sum(T_optimal_rtpd.FRU_error(T_optimal_rtpd.FRU_error>=0)));
    frd_rtpd_over_knn = abs(sum(T_optimal_rtpd.FRD_error(T_optimal_rtpd.FRD_error<=0)));
    frp_rtpd_over_knn(i) = fru_rtpd_over_knn + frd_rtpd_over_knn;

    tmp = T_optimal_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
    fru_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
    tmp = T_optimal_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
    frd_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
    frp_rtpd_freqshort_knn_hd(i) = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
    fprintf('n=%g, reliability=%.2f, over=%.2f\n', 5*i, frp_rtpd_freqshort_knn_hd(i), frp_rtpd_over_knn(i));
end
plot(frp_rtpd_freqshort_knn_hd, frp_rtpd_over_knn./1E3, '-^');
T_metric_over_opt.pcaknn = frp_rtpd_over_knn(:);
T_metric_pshort_opt.pcaknn = frp_rtpd_freqshort_knn_hd(:);
clearvars -except 'T_metric_over' 'T_metric_pshort' 'h' 'T_metric_over_opt' 'T_metric_pshort_opt';
%%
figure();
h = [];
h(1) = scatter(T_metric_pshort.baseline(:), T_metric_over.baseline(:)./1E3, 60, 'ko');
set(h(1), 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', 'k');
% set(h(1), 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', [.5, .5, .5], 'Alpha', .5, 'MarkerEdgeAlpha', .5, 'MarkerSize', 6);

hold on;
h(2) = scatter(T_metric_pshort.knn1d(:), T_metric_over.knn1d(:)./1E3, 30, '.');
set(h(2), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

h(numel(h) + 1) = plot(T_metric_pshort_opt.knnopt, T_metric_over_opt.knnopt./1E3, 'Marker', '^', 'Color', [.85, .33, .1]); % Dynamic kNN 
text(0.08, 340, 'N=30', 'FontSize', 18, 'Color', [.85, .33, .1]);
text(0.08, 340, 'N=5', 'FontSize', 18, 'Color', [.85, .33, .1]);

for classifier = 1: size(T_metric_pshort.pcaknn, 2)
    h(numel(h) + 1) = plot(T_metric_pshort.pcaknn(:, classifier), T_metric_over.pcaknn(:, classifier)./1E3, '--^');
end
set(h(numel(h)-2), 'Color', [0.68, 0.85, 1], 'MarkerFaceColor', [0.68, 0.85, 1], 'Marker', '^', 'MarkerSize', 6);
set(h(numel(h)-1), 'Color', [0.36, 0.36, 0.99], 'MarkerFaceColor', [0.36, 0.36, 0.99], 'Marker', '^', 'MarkerSize', 6);
set(h(numel(h)), 'Color', [0.44, 0.15, 0.77], 'MarkerFaceColor', [0.44, 0.15, 0.77], 'Marker', '^', 'MarkerSize', 6);
   
h(numel(h) + 1) = plot(T_metric_pshort_opt.pcaknn, T_metric_over_opt.pcaknn./1E3, 'Color', [0, 0.61, 0], 'Marker', '^'); % Dynamic PCA-kNN
text(0.08, 340, 'N=30', 'FontSize', 18, 'Color', [0, 0.61, 0]);
text(0.08, 340, 'N=5', 'FontSize', 18, 'Color', [0, 0.61, 0]);

ylabel('FRP over supply (GWh)'); xlabel('Probability of FRP shortage');
karray = 5:5:60;
xline(T_metric_pshort.baseline(karray==30), 'Color','red', 'LineStyle','--');
yline(T_metric_over.baseline(karray==30)./1E3, 'Color','red', 'LineStyle','--');
set(findall(gcf,'-property','FontSize'),'FontSize',22);
box on;
legend(h(:), {'Baseline', 'kNN 1-site', 's=2, Opt.', 'PCA 1-d', 'PCA 2-d' ,'PCA 3-d', 'Opt. PCA'}, 'FontSize', 14);
% text(0.08, 340, 'PCA 1-d', 'FontSize', 18, 'Color', [0.68, 0.85, 1]);
% text(0.08, 340, 'PCA 2-d', 'FontSize', 18, 'Color', [0.36, 0.36, 0.99]);
% text(0.08, 340, 'PCA 3-d', 'FontSize', 18, 'Color', [0.53, 0.36, 1]);
% text(0.08, 340, 'PCA dynamic', 'FontSize', 18, 'Color', [0, 0.61, 0]);
ylim([180, 360]);
xlim([0.05, 0.14])
grid on

clearvars -except 'T_metric_over' 'T_metric_pshort' 'h' 'T_metric_over_opt' 'T_metric_pshort_opt' 'karray';

%% Tiled plot
figure();
t = tiledlayout(1, 2);
ax1 = nexttile;
dimension = 1;

switch dimension
    case 1
        load knn_post_rtpd_puresolar_2_complete;
    case 2
        load('knn_post_rtpd_puresolar_2_complete.2dim.mat');
end
%
% figure();
hold on;

hh = nan(5, 1);
for s = 1:5
    fprintf('s=%g\n', s);
    for i = 1 : 6
        T_optimal_rtpd = cell_optimalknn_rtpd{i, s};
        f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_optimal_rtpd.HOUR_START), 'error_max'}; % 15-min
        f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_optimal_rtpd.HOUR_START), 'error_min'}; % 15-min
        T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
        T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});
        T_optimal_rtpd = [T_optimal_rtpd T_errormax_rtpd T_errormin_rtpd];
        T_optimal_rtpd.FRU_NEED = max(T_optimal_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, [], 2);
        T_optimal_rtpd.FRD_NEED = min(T_optimal_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, [], 2);
        
        % Calculate FRU/FRD errors for that hour and each 15-min interval
        T_fruerror_rtpd = array2table(T_optimal_rtpd.FRU-T_optimal_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
        T_frderror_rtpd = array2table(T_optimal_rtpd.FRD-T_optimal_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
        T_optimal_rtpd = [T_optimal_rtpd T_fruerror_rtpd T_frderror_rtpd];
        T_optimal_rtpd.FRU_error = T_optimal_rtpd.FRU - T_optimal_rtpd.FRU_NEED;
        T_optimal_rtpd.FRD_error = T_optimal_rtpd.FRD - T_optimal_rtpd.FRD_NEED;
        
        % Calculate evaluation metrics: Total oversupply
        fru_rtpd_over_knn = abs(sum(T_optimal_rtpd.FRU_error(T_optimal_rtpd.FRU_error>=0)));
        frd_rtpd_over_knn = abs(sum(T_optimal_rtpd.FRD_error(T_optimal_rtpd.FRD_error<=0)));
        frp_rtpd_over_dknn(i) = fru_rtpd_over_knn + frd_rtpd_over_knn;

        tmp = T_optimal_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
        fru_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
        tmp = T_optimal_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
        frd_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
        frp_rtpd_freqshort_dknn_hd(i) = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
        fprintf('n=%g, reliability=%.2f, over=%.2f\n', 5*i, frp_rtpd_freqshort_dknn_hd(i), frp_rtpd_over_dknn(i));
    end
    hh(s) = plot(frp_rtpd_freqshort_dknn_hd, frp_rtpd_over_dknn./1E3, '-^');
end
%
for i = 1 : 12
    T_baseline_rtpd = cell_baseline_rtpd{i, 1};
    T_baseline_rtpd = T_baseline_rtpd(T_baseline_rtpd.HOUR_START.Month==2, :);
    f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
    f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
%     T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
%     T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});
%     T_baseline_rtpd = [T_baseline_rtpd T_errormax_rtpd T_errormin_rtpd];
    T_baseline_rtpd.FRU_NEED = max(T_baseline_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, [], 2);
    T_baseline_rtpd.FRD_NEED = min(T_baseline_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, [], 2);
    
    % Calculate FRU/FRD errors for that hour and each 15-min interval
    T_fruerror_rtpd = array2table(T_baseline_rtpd.FRU-T_baseline_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
    T_frderror_rtpd = array2table(T_baseline_rtpd.FRD-T_baseline_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
    T_baseline_rtpd = [T_baseline_rtpd T_fruerror_rtpd T_frderror_rtpd];
    T_baseline_rtpd.FRU_error = T_baseline_rtpd.FRU - T_baseline_rtpd.FRU_NEED;
    T_baseline_rtpd.FRD_error = T_baseline_rtpd.FRD - T_baseline_rtpd.FRD_NEED;
    
    % Calculate evaluation metrics: Total oversupply
    fru_rtpd_over_knn = abs(sum(T_baseline_rtpd.FRU_error(T_baseline_rtpd.FRU_error>=0)));
    frd_rtpd_over_knn = abs(sum(T_baseline_rtpd.FRD_error(T_baseline_rtpd.FRD_error<=0)));
    frp_rtpd_over_baseline(i) = fru_rtpd_over_knn + frd_rtpd_over_knn;

    tmp = T_baseline_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
    fru_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
    tmp = T_baseline_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
    frd_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
    frp_rtpd_freqshort_baseline_hd(i) = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
    fprintf('n=%g, reliability=%.2f, over=%.2f\n', 5*i, frp_rtpd_freqshort_baseline_hd(i), frp_rtpd_over_baseline(i));
end
h_baseline = scatter(frp_rtpd_freqshort_baseline_hd, frp_rtpd_over_baseline./1E3, 60, 'ko');
set(h_baseline, 'MarkerFaceColor', 'k');
set(h_baseline, 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
xline(frp_rtpd_freqshort_baseline_hd(6), 'Color','red', 'LineStyle','--');
yline(frp_rtpd_over_baseline(6)./1E3, 'Color','red', 'LineStyle','--');

%
switch dimension
    case 1
        load('knn_post_rtpd_puresolar_2.5site.mat', 'cell_optimalknn_rtpd');
    case 2
        load('knn_post_rtpd_puresolar_2.5site.2dim.mat', 'cell_optimalknn_rtpd');
end
s = 1;
for i = 1 : 6
    T_optimal_rtpd = cell_optimalknn_rtpd{i, s};
    f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_optimal_rtpd.HOUR_START), 'error_max'}; % 15-min
    f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_optimal_rtpd.HOUR_START), 'error_min'}; % 15-min
    T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
    T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});
    T_optimal_rtpd = [T_optimal_rtpd T_errormax_rtpd T_errormin_rtpd];
    T_optimal_rtpd.FRU_NEED = max(T_optimal_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, [], 2);
    T_optimal_rtpd.FRD_NEED = min(T_optimal_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, [], 2);
    % Calculate FRU/FRD errors for that hour and each 15-min interval
    T_fruerror_rtpd = array2table(T_optimal_rtpd.FRU-T_optimal_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
    T_frderror_rtpd = array2table(T_optimal_rtpd.FRD-T_optimal_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
    T_optimal_rtpd = [T_optimal_rtpd T_fruerror_rtpd T_frderror_rtpd];
    T_optimal_rtpd.FRU_error = T_optimal_rtpd.FRU - T_optimal_rtpd.FRU_NEED;
    T_optimal_rtpd.FRD_error = T_optimal_rtpd.FRD - T_optimal_rtpd.FRD_NEED;
    % Calculate evaluation metrics: Total oversupply
    fru_rtpd_over_knn = abs(sum(T_optimal_rtpd.FRU_error(T_optimal_rtpd.FRU_error>=0)));
    frd_rtpd_over_knn = abs(sum(T_optimal_rtpd.FRD_error(T_optimal_rtpd.FRD_error<=0)));
    frp_rtpd_over_dknn(i) = fru_rtpd_over_knn + frd_rtpd_over_knn;
    tmp = T_optimal_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
    fru_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
    tmp = T_optimal_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
    frd_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
    frp_rtpd_freqshort_dknn_hd(i) = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
    fprintf('n=%g, reliability=%.2f, over=%.2f\n', 5*i, frp_rtpd_freqshort_dknn_hd(i), frp_rtpd_over_dknn(i));
end
% h_mean = plot(frp_rtpd_freqshort_dknn_hd, frp_rtpd_over_dknn./1E3, '-k^');

% xlim([0.05, 0.16]);
% ylim([180, 360]);
ylim([180, 360]);
xlim([0.05, 0.14])

box on;
grid on;
% set(findall(gcf,'-property','FontSize'),'FontSize',22);
text(frp_rtpd_freqshort_dknn_hd(1), frp_rtpd_over_dknn(1)./1E3, 'N=5', 'FontSize', 18, 'Color', 'k');
text(frp_rtpd_freqshort_dknn_hd(end), frp_rtpd_over_dknn(end)./1E3, 'N=30', 'FontSize', 18, 'Color', 'k');
text(frp_rtpd_freqshort_dknn_hd(1), frp_rtpd_over_dknn(1)./1E3, '(a)', 'FontSize', 18, 'Color', 'k');
xlabel('Frequency of reserve shortage', 'FontSize', 22);
ylabel('Oversupply (GWh)', 'FontSize', 22);
legend([hh; h_baseline], {'Site 1', 'Site 2', 'Site 3', 'Site 4', 'Site 5', 'Baseline'}, 'FontSize', 14);
legend boxon;

ax2 = nexttile;
h = [];
h(1) = scatter(T_metric_pshort.baseline(:), T_metric_over.baseline(:)./1E3, 60, 'ko');
set(h(1), 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', 'k');
% set(h(1), 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', [.5, .5, .5], 'Alpha', .5, 'MarkerEdgeAlpha', .5, 'MarkerSize', 6);

hold on;
% h(2) = scatter(T_metric_pshort.knn1d(:), T_metric_over.knn1d(:)./1E3, 30, '.');
% set(h(2), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

h(numel(h) + 1) = plot(T_metric_pshort_opt.knnopt, T_metric_over_opt.knnopt./1E3, 'Marker', '^', 'Color', [.85, .33, .1], 'MarkerFaceColor', [.85, .33, .1]); % Dynamic kNN 
text(0.08, 340, 'N=30', 'FontSize', 18, 'Color', [.85, .33, .1]);
text(0.08, 340, 'N=5', 'FontSize', 18, 'Color', [.85, .33, .1]);

for classifier = 1: size(T_metric_pshort.pcaknn, 2)
    h(numel(h) + 1) = plot(T_metric_pshort.pcaknn(:, classifier), T_metric_over.pcaknn(:, classifier)./1E3, '--^');
end
set(h(numel(h)-2), 'Color', [0.68, 0.85, 1], 'MarkerFaceColor', [0.68, 0.85, 1], 'Marker', '^', 'MarkerSize', 6);
set(h(numel(h)-1), 'Color', [0.36, 0.36, 0.99], 'MarkerFaceColor', [0.36, 0.36, 0.99], 'Marker', '^', 'MarkerSize', 6);
set(h(numel(h)), 'Color', [0.44, 0.15, 0.77], 'MarkerFaceColor', [0.44, 0.15, 0.77], 'Marker', '^', 'MarkerSize', 6);
   
h(numel(h) + 1) = plot(T_metric_pshort_opt.pcaknn, T_metric_over_opt.pcaknn./1E3, 'Color', [0, 0.61, 0],'MarkerFaceColor', [0, 0.61, 0], 'Marker', '^'); % Dynamic PCA-kNN
text(0.08, 340, 'N=30', 'FontSize', 18, 'Color', [0, 0.61, 0]);
text(0.08, 340, 'N=5', 'FontSize', 18, 'Color', [0, 0.61, 0]);
text(0.08, 340, '(b)', 'FontSize', 18, 'Color', 'k');


% ylabel('FRP over supply (GWh)'); 
xlabel('Probability of FRP shortage');
karray = 5:5:60;
xline(T_metric_pshort.baseline(karray==30), 'Color','red', 'LineStyle','--');
yline(T_metric_over.baseline(karray==30)./1E3, 'Color','red', 'LineStyle','--');
box on;
legend(h(:), {'Baseline', 's=2, Opt.', 'PCA 1-d', 'PCA 2-d' ,'PCA 3-d', 'Opt. PCA'}, 'FontSize', 14);
% text(0.08, 340, 'PCA 1-d', 'FontSize', 18, 'Color', [0.68, 0.85, 1]);
% text(0.08, 340, 'PCA 2-d', 'FontSize', 18, 'Color', [0.36, 0.36, 0.99]);
% text(0.08, 340, 'PCA 3-d', 'FontSize', 18, 'Color', [0.53, 0.36, 1]);
% text(0.08, 340, 'PCA dynamic', 'FontSize', 18, 'Color', [0, 0.61, 0]);
ylim([180, 360]);
xlim([0.05, 0.14])
grid on

set(findall(gcf,'-property','FontSize'),'FontSize',14);
yticklabels(ax2,{});
t.TileSpacing = 'compact';
t.Padding = 'none';
clearvars -except 'T_metric_over' 'T_metric_pshort' 'h' 'T_metric_over_opt' 'T_metric_pshort_opt' 'karray';

%% Load Cong's data
load('knn_rtpd_puresolar_2.mat', 'cell_baseline_rtpd')
T_baseline_rtpd = cell_baseline_rtpd{karray == 30, 1};
dirhome = pwd;

% ML files
dirwork = '.\Cong_results\ML2';
cd(dirwork);
files_cong = dir('Predictions-1Feature*Site_2.csv');
files_abspath = cell(size(files_cong, 1), 1);
files_dim = nan(size(files_cong, 1), 1);
files_isup = nan(size(files_cong, 1), 1);
for i = 1: size(files_cong, 1)
    fname = files_cong(i).name;
    fdir  = files_cong(i).folder;
    files_abspath{i} = fullfile(fdir, fname);
    if contains(fname, 'FRU')
        files_isup(i) = 1;
    else
        files_isup(i) = 0;
    end
    if contains(fname, '1Feature')
        files_dim(i) = 1;
    elseif contains(fname, '2Feature')
        files_dim(i) = 2;
    elseif contains(fname, '3Feature')
        files_dim(i) = 3;
    elseif contains(fname, '4Feature')
        files_dim(i) = 4;
    end
end
T_dlresults = array2table(files_abspath, 'VariableNames', {'ABSPATH'});
T_dlresults{:, 'ISUP'} = files_isup;
T_dlresults{:, 'DIM'} = files_dim;
cd(dirhome);

allfeatures = {'mean_k', 'std_k', 'v_k', 'mean_w', 'std_w', 'v_w', 'mean_kpv', 'std_kpv', 'v_kpv', 'mean_wpv', 'std_wpv', 'v_wpv'};

% CNN files
dirwork = '.\Cong_results\CNN2';
cd(dirwork);
files_cnn = dir('Predictions_1Feature*');
for i = 1:size(files_cnn, 1)
    fname = files_cnn(i).name;
    fdir  = files_cnn(i).folder;
    files_abspath{i} = fullfile(fdir, fname);
    if contains(fname, 'FRU')
        files_isup(i) = 1;
    else
        files_isup(i) = 0;
    end
    if contains(fname, '1Feature')
        files_dim(i) = 1;
    elseif contains(fname, '2Feature')
        files_dim(i) = 2;
    elseif contains(fname, '3Feature')
        files_dim(i) = 3;
    elseif contains(fname, '4Feature')
        files_dim(i) = 4;
    end
end
T_cnnresults = array2table(files_abspath, 'VariableNames', {'ABSPATH'});
T_cnnresults{:, 'ISUP'} = files_isup;
T_cnnresults{:, 'DIM'} = files_dim;
cd(dirhome);

%% Load results from Cong: ML
cell_dl5site_rtpd = cell(height(T_dlresults), 1);
for i = 1: height(T_dlresults)
    tmp = readtable(T_dlresults{i, 'ABSPATH'}{1});
    tmp.df_merge_df_tidx_clean_which_test_.TimeZone = 'UTC';
    tmp.Properties.VariableNames{1} = 'HOUR_START';
    cell_dl5site_rtpd{i} = tmp;
    fprintf('Size: %2g, %2g, File: %s\n', size(tmp, 1), size(tmp, 2), T_dlresults{i, 'ABSPATH'}{1});
end

% Load results from Cong: CNN
cell_cnn5site_rtpd = cell(height(T_cnnresults), 1); % Time axis is missing
for i = 1: height(T_cnnresults)
    tmp = readtable(T_cnnresults{i, 'ABSPATH'}{1});
    tmp.HOUR_START = datetime(tmp.HOUR_START, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    tmp.HOUR_START = tmp.HOUR_START - duration(7, 0, 0);
    cell_cnn5site_rtpd{i} = tmp;
    fprintf('Size: %2g, %2g, File: %s\n', size(tmp, 1), size(tmp, 2), T_cnnresults{i, 'ABSPATH'}{1});
end


allmethods = {'svm_r', 'svm_l', 'svm_l2', 'svm_p', 'svm_p2', 'ann1', 'ann2', 'ann4', 'ann5', 'ann6', 'gbm1', 'gbm2', 'gbm3', 'gbm4', 'rf1', 'rf2'};
% Post-process: fill Cong's missing values with baseline values: ML
cell_dl5site_rtpd1 = cell(size(cell_dl5site_rtpd));
for i = 1: numel(cell_dl5site_rtpd)
    if ~isempty(cell_dl5site_rtpd{i})
        if T_dlresults{i, 'ISUP'} == 1
            column_baseline = 'FRU'; % Use CAISO baseline if DL has no results
        else
            column_baseline = 'FRD'; % Use CAISO baseline if DL has not results
        end
        tmp = T_baseline_rtpd(:, {'HOUR_START', column_baseline});

        for j = 1: length(allmethods)
            tmp{:, allmethods{j}} = tmp{:, column_baseline};
            [Lia, Locb] = ismember(cell_dl5site_rtpd{i}.HOUR_START, tmp.HOUR_START);
            tmp{Locb(Lia), allmethods{j}} = cell_dl5site_rtpd{i}{Lia, allmethods{j}};
            cell_dl5site_rtpd1{i} = tmp;
        end
    end
end

% Post-process: fill Cong's missing values with baseline values: CNN
cell_cnn5site_rtpd1 = cell(size(cell_cnn5site_rtpd));
for i = 1: numel(cell_cnn5site_rtpd)
    if ~isempty(cell_cnn5site_rtpd{i})
        if T_cnnresults{i, 'ISUP'} == 1
            column_baseline = 'FRU'; % Use CAISO baseline if DL has no results
        else
            column_baseline = 'FRD'; % Use CAISO baseline if DL has not results
        end
        tmp = T_baseline_rtpd(:, {'HOUR_START', column_baseline});

        tmp{:, 'CNN'} = tmp{:, column_baseline};
        [Lia, Locb] = ismember(cell_cnn5site_rtpd{i}.HOUR_START, tmp.HOUR_START);
        tmp{Locb(Lia), 'CNN'} = cell_cnn5site_rtpd{i}{Lia, 'Predicted'};
        cell_cnn5site_rtpd1{i} = tmp;
    end
end

%% Calculate evaluating metrics for Cong: 1D
load('knn_rtpd_puresolar_2.mat', 'T_rtpd', 'T_baseline_rtpd', 'karray', 'cell_baseline_rtpd', 'cell_results_rtpd');

% This is the actual need of FRP
f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
fru_need_rtpd = max(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4), [], 1)';
f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
frd_need_rtpd = min(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4), [], 1)';
T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});

frd_rtpd_freqshort_knn_hd = [];
fru_rtpd_freqshort_knn_hd = [];
fru_rtpd_over_knn = [];
frd_rtpd_over_knn = [];
% Calculate metrics: 1D feature
for i = 1: numel(allfeatures)
    % ML methods
    T_cong_d = cell_dl5site_rtpd1{contains(T_dlresults.ABSPATH, strcat('Predictions-1Feature_', allfeatures{i}, '_FRD_Site_2.csv'))};
    T_cong_u = cell_dl5site_rtpd1{contains(T_dlresults.ABSPATH, strcat('Predictions-1Feature_', allfeatures{i}, '_FRU_Site_2.csv'))};
    for j = 1: numel(allmethods)
        T_results_rtpd = table(T_cong_u{:, allmethods{j}}, T_cong_d{:, allmethods{j}}, 'VariableNames', {'FRU', 'FRD'});

        T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
        T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

        T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
        T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
        T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
        T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
        T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

        tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
        fru_rtpd_freqshort_knn_hd(i, j) = sum(tmp(:))/numel(tmp);
        tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
        frd_rtpd_freqshort_knn_hd(i, j) = sum(tmp(:))/numel(tmp);

        fru_rtpd_over_knn(i, j)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
        frd_rtpd_over_knn(i, j)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
    end
    
    % CNN
    T_cong_d = cell_cnn5site_rtpd1{contains(T_cnnresults.ABSPATH, strcat('Predictions_1Feature-', allfeatures{i}, '-FRD.csv'))};
    T_cong_u = cell_cnn5site_rtpd1{contains(T_cnnresults.ABSPATH, strcat('Predictions_1Feature-', allfeatures{i}, '-FRU.csv'))};
    T_results_rtpd = table(T_cong_u{:, 'CNN'}, T_cong_d{:, 'CNN'}, 'VariableNames', {'FRU', 'FRD'});

    T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
    T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

    T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
    T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
    T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
    T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
    T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

    tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
    fru_rtpd_freqshort_knn_hd(i, numel(allmethods)+1) = sum(tmp(:))/numel(tmp);
    tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
    frd_rtpd_freqshort_knn_hd(i, numel(allmethods)+1) = sum(tmp(:))/numel(tmp);

    fru_rtpd_over_knn(i, numel(allmethods)+1)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
    frd_rtpd_over_knn(i, numel(allmethods)+1)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));

end
frp_rtpd_freqshort_knn_hd = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
frp_rtpd_over_knn = fru_rtpd_over_knn + frd_rtpd_over_knn;

T_frp_rtpd_freqshort_knn_hd = [array2table(allfeatures(:), 'VariableNames', {'Features'}) array2table(frp_rtpd_freqshort_knn_hd, 'VariableNames', [allmethods, 'CNN'])];
T_frp_rtpd_over_knn = [array2table(allfeatures(:), 'VariableNames', {'Features'}) array2table(frp_rtpd_over_knn, 'VariableNames', [allmethods, 'CNN'])];

%% Load files :2D
dirwork = '.\Cong_results\ML2';
cd(dirwork);
files_cong = dir('Predictions-2Feature*MultiSite_.csv');
files_abspath = cell(size(files_cong, 1), 1);
files_dim = nan(size(files_cong, 1), 1);
files_isup = nan(size(files_cong, 1), 1);
for i = 1: size(files_cong, 1)
    fname = files_cong(i).name;
    fdir  = files_cong(i).folder;
    files_abspath{i} = fullfile(fdir, fname);
    if contains(fname, 'FRU')
        files_isup(i) = 1;
    else
        files_isup(i) = 0;
    end
    if contains(fname, '1Feature')
        files_dim(i) = 1;
    elseif contains(fname, '2Feature')
        files_dim(i) = 2;
    elseif contains(fname, '3Feature')
        files_dim(i) = 3;
    elseif contains(fname, '4Feature')
        files_dim(i) = 4;
    end
end
T_dlresults_2d = array2table(files_abspath, 'VariableNames', {'ABSPATH'});
T_dlresults_2d{:, 'ISUP'} = files_isup;
T_dlresults_2d{:, 'DIM'} = files_dim;
cd(dirhome);

allfeatures_2d = {'mean_k_mean_w', 'mean_k_v_k', 'mean_kpv_mean_wpv', 'mean_wpv_v_kpv'};

% CNN files
dirwork = '.\Cong_results\CNN2';
cd(dirwork);
files_cnn = dir('Predictions_2Feature*');
for i = 1:size(files_cnn, 1)
    fname = files_cnn(i).name;
    fdir  = files_cnn(i).folder;
    files_abspath{i} = fullfile(fdir, fname);
    if contains(fname, 'FRU')
        files_isup(i) = 1;
    else
        files_isup(i) = 0;
    end
    if contains(fname, '1Feature')
        files_dim(i) = 1;
    elseif contains(fname, '2Feature')
        files_dim(i) = 2;
    elseif contains(fname, '3Feature')
        files_dim(i) = 3;
    elseif contains(fname, '4Feature')
        files_dim(i) = 4;
    end
end
T_cnnresults_2d = array2table(files_abspath, 'VariableNames', {'ABSPATH'});
T_cnnresults_2d{:, 'ISUP'} = files_isup;
T_cnnresults_2d{:, 'DIM'} = files_dim;
cd(dirhome);

% Load results from Cong: ML 2d
cell_dl5site_rtpd_2d = cell(height(T_dlresults_2d), 1);
for i = 1: height(T_dlresults_2d)
    tmp = readtable(T_dlresults_2d{i, 'ABSPATH'}{1});
    tmp.df_merge_df_tidx_clean_which_test_.TimeZone = 'UTC';
    tmp.Properties.VariableNames{1} = 'HOUR_START';
    cell_dl5site_rtpd_2d{i} = tmp;
    fprintf('Size: %2g, %2g, File: %s\n', size(tmp, 1), size(tmp, 2), T_dlresults_2d{i, 'ABSPATH'}{1});
end

% Load results from Cong: CNN 2d
cell_cnn5site_rtpd_2d = cell(height(T_cnnresults_2d), 1); % Time axis is missing
for i = 1: height(T_cnnresults_2d)
    tmp = readtable(T_cnnresults_2d{i, 'ABSPATH'}{1});
    tmp.HOUR_START = datetime(tmp.HOUR_START, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    tmp.HOUR_START = tmp.HOUR_START - duration(7, 0, 0);
    cell_cnn5site_rtpd_2d{i} = tmp;
    fprintf('Size: %2g, %2g, File: %s\n', size(tmp, 1), size(tmp, 2), T_cnnresults_2d{i, 'ABSPATH'}{1});
end

% Post-process: fill Cong's missing values with baseline values: ML 2d
cell_dl5site_rtpd1_2d = cell(size(cell_dl5site_rtpd_2d));
for i = 1: numel(cell_dl5site_rtpd_2d)
    if ~isempty(cell_dl5site_rtpd_2d{i})
        if T_dlresults_2d{i, 'ISUP'} == 1
            column_baseline = 'FRU'; % Use CAISO baseline if DL has no results
        else
            column_baseline = 'FRD'; % Use CAISO baseline if DL has not results
        end
        tmp = T_baseline_rtpd(:, {'HOUR_START', column_baseline});

        for j = 1: length(allmethods)
            tmp{:, allmethods{j}} = tmp{:, column_baseline};
            [Lia, Locb] = ismember(cell_dl5site_rtpd_2d{i}.HOUR_START, tmp.HOUR_START);
            tmp{Locb(Lia), allmethods{j}} = cell_dl5site_rtpd_2d{i}{Lia, allmethods{j}};
            cell_dl5site_rtpd1_2d{i} = tmp;
        end
    end
end

% Post-process: fill Cong's missing values with baseline values: CNN 2d
cell_cnn5site_rtpd1_2d = cell(size(cell_cnn5site_rtpd_2d));
for i = 1: numel(cell_cnn5site_rtpd_2d)
    if ~isempty(cell_cnn5site_rtpd_2d{i})
        if T_cnnresults_2d{i, 'ISUP'} == 1
            column_baseline = 'FRU'; % Use CAISO baseline if DL has no results
        else
            column_baseline = 'FRD'; % Use CAISO baseline if DL has not results
        end
        tmp = T_baseline_rtpd(:, {'HOUR_START', column_baseline});

        tmp{:, 'CNN'} = tmp{:, column_baseline};
        [Lia, Locb] = ismember(cell_cnn5site_rtpd_2d{i}.HOUR_START, tmp.HOUR_START);
        tmp{Locb(Lia), 'CNN'} = cell_cnn5site_rtpd_2d{i}{Lia, 'Predicted'};
        cell_cnn5site_rtpd1_2d{i} = tmp;
    end
end

%% Calculate evaluating metrics for Cong: 2D
load('knn_rtpd_puresolar_2.mat', 'T_rtpd', 'T_baseline_rtpd', 'karray', 'cell_baseline_rtpd', 'cell_results_rtpd');

% This is the actual need of FRP
f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
fru_need_rtpd = max(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4), [], 1)';
f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
frd_need_rtpd = min(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4), [], 1)';
T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});

frd_rtpd_freqshort_knn_hd = [];
fru_rtpd_freqshort_knn_hd = [];
fru_rtpd_over_knn = [];
frd_rtpd_over_knn = [];
% Calculate metrics: 1D feature
for i = 1: numel(allfeatures_2d)
    % ML methods
    T_cong_d = cell_dl5site_rtpd1_2d{contains(T_dlresults_2d.ABSPATH, strcat('Predictions-2Feature_', allfeatures_2d{i}, '_FRD_MultiSite_.csv'))};
    T_cong_u = cell_dl5site_rtpd1_2d{contains(T_dlresults_2d.ABSPATH, strcat('Predictions-2Feature_', allfeatures_2d{i}, '_FRU_MultiSite_.csv'))};
    for j = 1: numel(allmethods)
        T_results_rtpd = table(T_cong_u{:, allmethods{j}}, T_cong_d{:, allmethods{j}}, 'VariableNames', {'FRU', 'FRD'});

        T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
        T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

        T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
        T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
        T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
        T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
        T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

        tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
        fru_rtpd_freqshort_knn_hd(i, j) = sum(tmp(:))/numel(tmp);
        tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
        frd_rtpd_freqshort_knn_hd(i, j) = sum(tmp(:))/numel(tmp);

        fru_rtpd_over_knn(i, j)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
        frd_rtpd_over_knn(i, j)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
    end
    
    % CNN
    T_cong_d = cell_cnn5site_rtpd1_2d{contains(T_cnnresults_2d.ABSPATH, strcat('Predictions_2Feature-', allfeatures_2d{i}, '-FRD.csv'))};
    T_cong_u = cell_cnn5site_rtpd1_2d{contains(T_cnnresults_2d.ABSPATH, strcat('Predictions_2Feature-', allfeatures_2d{i}, '-FRU.csv'))};
    T_results_rtpd = table(T_cong_u{:, 'CNN'}, T_cong_d{:, 'CNN'}, 'VariableNames', {'FRU', 'FRD'});

    T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
    T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

    T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
    T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
    T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
    T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
    T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

    tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
    fru_rtpd_freqshort_knn_hd(i, numel(allmethods)+1) = sum(tmp(:))/numel(tmp);
    tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
    frd_rtpd_freqshort_knn_hd(i, numel(allmethods)+1) = sum(tmp(:))/numel(tmp);

    fru_rtpd_over_knn(i, numel(allmethods)+1)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
    frd_rtpd_over_knn(i, numel(allmethods)+1)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));

end
frp_rtpd_freqshort_knn_hd = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
frp_rtpd_over_knn = fru_rtpd_over_knn + frd_rtpd_over_knn;

T_frp_rtpd_freqshort_knn_hd_2d = [array2table(allfeatures_2d(:), 'VariableNames', {'Features'}) array2table(frp_rtpd_freqshort_knn_hd, 'VariableNames', [allmethods, 'CNN'])];
T_frp_rtpd_over_knn_2d = [array2table(allfeatures_2d(:), 'VariableNames', {'Features'}) array2table(frp_rtpd_over_knn, 'VariableNames', [allmethods, 'CNN'])];

%% Frontier of the ML/DL-based methods
figure();
h = [];
h(1) = scatter(T_metric_pshort.baseline(:), T_metric_over.baseline(:)./1E3, 60, 'ko');
set(h(1), 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', 'k');
% set(h(1), 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', [.5, .5, .5], 'Alpha', .5, 'MarkerEdgeAlpha', .5, 'MarkerSize', 6);

hold on;
% for i = 
h(2) = scatter(T_frp_rtpd_freqshort_knn_hd{:, [allmethods, 'CNN']}(:), T_frp_rtpd_over_knn{:, [allmethods, 'CNN']}(:)./1E3, 30, '.');
h(3) = scatter(T_frp_rtpd_freqshort_knn_hd_2d{:, [allmethods, 'CNN']}(:), T_frp_rtpd_over_knn_2d{:, [allmethods, 'CNN']}(:)./1E3, 30, '.');
% set(h(2), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

xline(T_metric_pshort.baseline(karray==30), 'Color','red', 'LineStyle','--');
yline(T_metric_over.baseline(karray==30)./1E3, 'Color','red', 'LineStyle','--');
%% What if we scale them up, 1-d

alpha = 1:0.1:3;
frd_rtpd_freqshort_knn_hd = [];
fru_rtpd_freqshort_knn_hd = [];
fru_rtpd_over_knn = [];
frd_rtpd_over_knn = [];
% Calculate metrics: 1D feature
for m = 1: numel(alpha)
    for i = 1: numel(allfeatures)
        % ML methods
        T_cong_d = cell_dl5site_rtpd1{contains(T_dlresults.ABSPATH, strcat('Predictions-1Feature_', allfeatures{i}, '_FRD_Site_2.csv'))};
        T_cong_u = cell_dl5site_rtpd1{contains(T_dlresults.ABSPATH, strcat('Predictions-1Feature_', allfeatures{i}, '_FRU_Site_2.csv'))};
        for j = 1: numel(allmethods)
            T_results_rtpd = table(T_cong_u{:, allmethods{j}}.*alpha(m), T_cong_d{:, allmethods{j}}.*alpha(m), 'VariableNames', {'FRU', 'FRD'});

            T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
            T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

            T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
            T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
            T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
            T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
            T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

            tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
            fru_rtpd_freqshort_knn_hd(i, j, m) = sum(tmp(:))/numel(tmp);
            tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
            frd_rtpd_freqshort_knn_hd(i, j, m) = sum(tmp(:))/numel(tmp);

            fru_rtpd_over_knn(i, j, m)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
            frd_rtpd_over_knn(i, j, m)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
        end

        % CNN
        T_cong_d = cell_cnn5site_rtpd1{contains(T_cnnresults.ABSPATH, strcat('Predictions_1Feature-', allfeatures{i}, '-FRD.csv'))};
        T_cong_u = cell_cnn5site_rtpd1{contains(T_cnnresults.ABSPATH, strcat('Predictions_1Feature-', allfeatures{i}, '-FRU.csv'))};
        T_results_rtpd = table(T_cong_u{:, 'CNN'}.*alpha(m), T_cong_d{:, 'CNN'}.*alpha(m), 'VariableNames', {'FRU', 'FRD'});

        T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
        T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

        T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
        T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
        T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
        T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
        T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

        tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
        fru_rtpd_freqshort_knn_hd(i, numel(allmethods)+1, m) = sum(tmp(:))/numel(tmp);
        tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
        frd_rtpd_freqshort_knn_hd(i, numel(allmethods)+1, m) = sum(tmp(:))/numel(tmp);

        fru_rtpd_over_knn(i, numel(allmethods)+1, m)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
        frd_rtpd_over_knn(i, numel(allmethods)+1, m)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
    end
end

frp_rtpd_freqshort_knn_hd = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
frp_rtpd_over_knn = fru_rtpd_over_knn + frd_rtpd_over_knn;

ar_frp_rtpd_freqshort_knn_hd = frp_rtpd_freqshort_knn_hd;
ar_frp_rtpd_over_knn = frp_rtpd_over_knn;

%% What if we scale them up, 2-D
frd_rtpd_freqshort_knn_hd = [];
fru_rtpd_freqshort_knn_hd = [];
fru_rtpd_over_knn = [];
frd_rtpd_over_knn = [];
% Calculate metrics: 1D feature
for m = 1: numel(alpha)

    for i = 1: numel(allfeatures_2d)
        % ML methods
        T_cong_d = cell_dl5site_rtpd1_2d{contains(T_dlresults_2d.ABSPATH, strcat('Predictions-2Feature_', allfeatures_2d{i}, '_FRD_MultiSite_.csv'))};
        T_cong_u = cell_dl5site_rtpd1_2d{contains(T_dlresults_2d.ABSPATH, strcat('Predictions-2Feature_', allfeatures_2d{i}, '_FRU_MultiSite_.csv'))};
        for j = 1: numel(allmethods)
            T_results_rtpd = table(T_cong_u{:, allmethods{j}}.*alpha(m), T_cong_d{:, allmethods{j}}.*alpha(m), 'VariableNames', {'FRU', 'FRD'});

            T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
            T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

            T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
            T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
            T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
            T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
            T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

            tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
            fru_rtpd_freqshort_knn_hd(i, j, m) = sum(tmp(:))/numel(tmp);
            tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
            frd_rtpd_freqshort_knn_hd(i, j, m) = sum(tmp(:))/numel(tmp);

            fru_rtpd_over_knn(i, j, m)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
            frd_rtpd_over_knn(i, j, m)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
        end

        % CNN
        T_cong_d = cell_cnn5site_rtpd1_2d{contains(T_cnnresults_2d.ABSPATH, strcat('Predictions_2Feature-', allfeatures_2d{i}, '-FRD.csv'))};
        T_cong_u = cell_cnn5site_rtpd1_2d{contains(T_cnnresults_2d.ABSPATH, strcat('Predictions_2Feature-', allfeatures_2d{i}, '-FRU.csv'))};
        T_results_rtpd = table(T_cong_u{:, 'CNN'}.*alpha(m), T_cong_d{:, 'CNN'}.*alpha(m), 'VariableNames', {'FRU', 'FRD'});

        T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
        T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

        T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
        T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
        T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
        T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
        T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

        tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
        fru_rtpd_freqshort_knn_hd(i, numel(allmethods)+1, m) = sum(tmp(:))/numel(tmp);
        tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
        frd_rtpd_freqshort_knn_hd(i, numel(allmethods)+1, m) = sum(tmp(:))/numel(tmp);

        fru_rtpd_over_knn(i, numel(allmethods)+1, m)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
        frd_rtpd_over_knn(i, numel(allmethods)+1, m)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));

    end
end
frp_rtpd_freqshort_knn_hd = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
frp_rtpd_over_knn = fru_rtpd_over_knn + frd_rtpd_over_knn;

ar_frp_rtpd_freqshort_knn_hd_2d = frp_rtpd_freqshort_knn_hd;
ar_frp_rtpd_over_knn_2d = frp_rtpd_over_knn;

%% Plot 
figure();
h = [];
hold on;
for i = 1: size(ar_frp_rtpd_freqshort_knn_hd, 1)
    for j = 1: size(ar_frp_rtpd_freqshort_knn_hd, 2)
        line(squeeze(ar_frp_rtpd_freqshort_knn_hd(i, j, :)), squeeze(ar_frp_rtpd_over_knn(i, j, :))./1E3, 'Color',[0.5, 0.5, 0.5]);
    end
end

for i = 1: size(ar_frp_rtpd_freqshort_knn_hd_2d, 1)
    for j = 1: size(ar_frp_rtpd_freqshort_knn_hd_2d, 2)
        h = line(squeeze(ar_frp_rtpd_freqshort_knn_hd_2d(i, j, :)), squeeze(ar_frp_rtpd_over_knn_2d(i, j, :))./1E3, 'Color',[0.5, 0.5, 0.5]);
    end
end


h(numel(h) + 1) = scatter(T_metric_pshort.baseline(:), T_metric_over.baseline(:)./1E3, 60, 'ko');
set(h(end), 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', 'k');
% set(h(1), 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', [.5, .5, .5], 'Alpha', .5, 'MarkerEdgeAlpha', .5, 'MarkerSize', 6);

hold on;
% h(2) = scatter(T_metric_pshort.knn1d(:), T_metric_over.knn1d(:)./1E3, 30, '.');
% set(h(2), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

h(numel(h) + 1) = plot(T_metric_pshort_opt.knnopt, T_metric_over_opt.knnopt./1E3, 'Marker', '^', 'Color', [.85, .33, .1], 'MarkerFaceColor', [.85, .33, .1]); % Dynamic kNN 
% text(0.08, 340, 'N=30', 'FontSize', 18, 'Color', [.85, .33, .1]);
% text(0.08, 340, 'N=5', 'FontSize', 18, 'Color', [.85, .33, .1]);

% for classifier = 1: size(T_metric_pshort.pcaknn, 2)
%     h(numel(h) + 1) = plot(T_metric_pshort.pcaknn(:, classifier), T_metric_over.pcaknn(:, classifier)./1E3, '--^');
% end
% set(h(numel(h)-2), 'Color', [0.68, 0.85, 1], 'MarkerFaceColor', [0.68, 0.85, 1], 'Marker', '^', 'MarkerSize', 6);
% set(h(numel(h)-1), 'Color', [0.36, 0.36, 0.99], 'MarkerFaceColor', [0.36, 0.36, 0.99], 'Marker', '^', 'MarkerSize', 6);
% set(h(numel(h)), 'Color', [0.44, 0.15, 0.77], 'MarkerFaceColor', [0.44, 0.15, 0.77], 'Marker', '^', 'MarkerSize', 6);
   
h(numel(h) + 1) = plot(T_metric_pshort_opt.pcaknn, T_metric_over_opt.pcaknn./1E3, 'Color', [0, 0.61, 0], 'Marker', '^', 'MarkerFaceColor', [0, 0.61, 0]); % Dynamic PCA-kNN
% text(0.08, 340, 'N=30', 'FontSize', 18, 'Color', [0, 0.61, 0]);
% text(0.08, 340, 'N=5', 'FontSize', 18, 'Color', [0, 0.61, 0]);

ylabel('FRP over supply (GWh)'); xlabel('Probability of FRP shortage');
karray = 5:5:60;
xline(T_metric_pshort.baseline(karray==30), 'Color','red', 'LineStyle','--');
yline(T_metric_over.baseline(karray==30)./1E3, 'Color','red', 'LineStyle','--');
% legend(h(:), {'Baseline', 'kNN 1-site', 's=2, Opt.', 'PCA 1-d', 'PCA 2-d' ,'PCA 3-d', 'Opt. PCA'}, 'FontSize', 14);
% text(0.08, 340, 'PCA 1-d', 'FontSize', 18, 'Color', [0.68, 0.85, 1]);
% text(0.08, 340, 'PCA 2-d', 'FontSize', 18, 'Color', [0.36, 0.36, 0.99]);
% text(0.08, 340, 'PCA 3-d', 'FontSize', 18, 'Color', [0.53, 0.36, 1]);
% text(0.08, 340, 'PCA dynamic', 'FontSize', 18, 'Color', [0, 0.61, 0]);

scatter(reshape(ar_frp_rtpd_freqshort_knn_hd(:, :, 1), numel(ar_frp_rtpd_freqshort_knn_hd(:, :, 1)), 1), reshape(ar_frp_rtpd_over_knn(:, :, 1), numel(ar_frp_rtpd_over_knn(:, :, 1)), 1)./1E3, 30, 'k.');
h(numel(h) + 1) = scatter(reshape(ar_frp_rtpd_freqshort_knn_hd_2d(:, :, 1), numel(ar_frp_rtpd_freqshort_knn_hd_2d(:, :, 1)), 1), reshape(ar_frp_rtpd_over_knn_2d(:, :, 1), numel(ar_frp_rtpd_over_knn_2d(:, :, 1)), 1)./1E3, 30, 'k.');
% set(h(2), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');


grid on
text(0.08, 340, '$\beta=1$', 'FontSize', 22, 'Color', 'k', 'Interpreter', 'latex');
text(0.08, 340, '$\beta=3$', 'FontSize', 22, 'Color', 'k', 'Interpreter', 'latex');
set(findall(gcf,'-property','FontSize'),'FontSize',22);
box on;

legend(h, {'scaled ML', 'Baseline', 'kNN', 'PCA-kNN', 'ML'})
rectangle('Position',[0.05, 180, 0.09, 180], 'EdgeColor','red');


% ylim([180, 360]);
% xlim([0.05, 0.14])

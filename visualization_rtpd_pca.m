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

h(numel(h) + 1) = plot(T_metric_pshort_opt.knnopt, T_metric_over_opt.knnopt./1E3, '-k^'); % Dynamic kNN 

for classifier = 1: size(T_metric_pshort.pcaknn, 2)
    h(numel(h) + 1) = plot(T_metric_pshort.pcaknn(:, classifier), T_metric_over.pcaknn(:, classifier)./1E3, '--^');
end
set(h(numel(h)-2), 'Color', [0.68, 0.85, 1], 'MarkerFaceColor', [0.68, 0.85, 1], 'Marker', '^', 'MarkerSize', 6);
set(h(numel(h)-1), 'Color', [0.36, 0.36, 0.99], 'MarkerFaceColor', [0.36, 0.36, 0.99], 'Marker', '^', 'MarkerSize', 6);
set(h(numel(h)), 'Color', [0.53, 0.36, 1], 'MarkerFaceColor', [0.53, 0.36, 1], 'Marker', '^', 'MarkerSize', 6);
   
h(numel(h) + 1) = plot(T_metric_pshort_opt.pcaknn, T_metric_over_opt.pcaknn./1E3, '-b^'); % Dynamic PCA-kNN

ylabel('FRP over supply (GWh)'); xlabel('Probability of FRP shortage');
xline(T_metric_pshort.baseline(karray==30), 'Color','red', 'LineStyle','--');
yline(T_metric_over.baseline(karray==30)./1E3, 'Color','red', 'LineStyle','--');
set(findall(gcf,'-property','FontSize'),'FontSize',14);
box on;
legend(h(1:3), {'Baseline', 'kNN 1-site', 'Dynamic 1-site'});
text(0.08, 340, 'PCA 1-d', 'FontSize', 14, 'Color', [0.68, 0.85, 1]);
text(0.08, 340, 'PCA 2-d', 'FontSize', 14, 'Color', [0.36, 0.36, 0.99]);
text(0.08, 340, 'PCA 3-d', 'FontSize', 14, 'Color', [0.53, 0.36, 1]);
text(0.08, 340, 'PCA dynamic', 'FontSize', 14, 'Color', 'b');

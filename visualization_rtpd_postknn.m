dimension = 1;

switch dimension
    case 1
        load knn_post_rtpd_puresolar_2_complete;
    case 2
        load('knn_post_rtpd_puresolar_2_complete.2dim.mat');
end
%%
figure();
hold on;

h = nan(5, 1);
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
    h(s) = plot(frp_rtpd_freqshort_dknn_hd, frp_rtpd_over_dknn./1E3, '-^');
end
%%
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

%%
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
h_mean = plot(frp_rtpd_freqshort_dknn_hd, frp_rtpd_over_dknn./1E3, '-k^');

xlim([0.05, 0.12]);
ylim([220, 360]);
box on;
grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',22);
text(frp_rtpd_freqshort_dknn_hd(1), frp_rtpd_over_dknn(1)./1E3, 'N=5', 'FontSize', 18, 'Color', 'k');
text(frp_rtpd_freqshort_dknn_hd(end), frp_rtpd_over_dknn(end)./1E3, 'N=30', 'FontSize', 18, 'Color', 'k');
xlabel('Frequency of reserve shortage', 'FontSize', 22);
ylabel('Oversupply (GWh)', 'FontSize', 22);
legend([h; h_mean; h_baseline], {'Site 1', 'Site 2', 'Site 3', 'Site 4', 'Site 5', 'Mean', 'Baseline'}, 'FontSize', 14);

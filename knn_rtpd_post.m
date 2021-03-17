% Post-processing kNN training data
% First, concatenation results from multiple months
% knn_1 = load('knn_rtpd_puresolar_1.2dim.mat');
% knn_2 = load('knn_rtpd_puresolar_2.2dim.mat');

% These tables should be the same from different mat files
T_rtpd = knn_1.T_rtpd;
T_pwr = knn_1.T_pwr;

T_baseline_rtpd = [knn_1.cell_baseline_rtpd{1, 1, 1}; knn_2.cell_baseline_rtpd{1, 1, 1}];

cell_baseline_rtpd = cell(size(knn_1.cell_baseline_rtpd));
cell_results_rtpd  = cell(size(knn_1.cell_results_rtpd));
f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});

for i1 = 1:size(cell_baseline_rtpd, 1)
    for i2 = 1:size(cell_baseline_rtpd, 2)
        for i3 = 1: size(cell_baseline_rtpd, 3)
            T_baseline_rtpd = [knn_1.cell_baseline_rtpd{i1, i2, i3}; knn_2.cell_baseline_rtpd{i1, i2, i3}];
            cell_baseline_rtpd{i1, i2, i3} = [T_baseline_rtpd T_errormax_rtpd T_errormin_rtpd];
        end
    end
end

for i1 = 1:size(cell_results_rtpd, 1)
    for i2 = 1:size(cell_results_rtpd, 2)
        for i3 = 1: size(cell_results_rtpd, 3)
            T_results_rtpd = [knn_1.cell_results_rtpd{i1, i2, i3}; knn_2.cell_results_rtpd{i1, i2, i3}];
            T_results_rtpd = [T_results_rtpd T_errormax_rtpd T_errormin_rtpd];

            T_results_rtpd.FRU_NEED = max(T_results_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, [], 2);
            T_results_rtpd.FRD_NEED = min(T_results_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, [], 2);
            T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_results_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
            T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_results_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
            T_results_rtpd = [T_results_rtpd T_fruerror_rtpd T_frderror_rtpd];
            T_results_rtpd.FRU_error = T_results_rtpd.FRU - T_results_rtpd.FRU_NEED;
            T_results_rtpd.FRD_error = T_results_rtpd.FRD - T_results_rtpd.FRD_NEED;

            T_results_rtpd.fru_over = abs(max(0, T_results_rtpd.FRU_error));
            T_results_rtpd.frd_over = abs(min(0, T_results_rtpd.FRD_error));
            T_results_rtpd.fru_short = abs(min(0, T_results_rtpd.FRU_error));
            T_results_rtpd.frd_short = abs(max(0, T_results_rtpd.FRD_error));
            T_results_rtpd.fru_imbalance = T_results_rtpd.fru_over + T_results_rtpd.fru_short;
            T_results_rtpd.frd_imbalance = T_results_rtpd.frd_over + T_results_rtpd.frd_short;
            T_results_rtpd.fru_freqshort_hd = mean(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2);
            T_results_rtpd.frd_freqshort_hd = mean(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2);

            cell_results_rtpd{i1, i2, i3} = T_results_rtpd;
        end
    end
end

clear knn_1 knn_2;
%%
% Calculate the evaluation metrics for the optimal selection of kNN parameters
ns = size(cell_results_rtpd, 2); % Number of sites
array_ndays = 5*[1, 2, 3, 4, 5, 6];
cell_optimalknn_rtpd = cell(numel(array_ndays), ns);
for indays = 1: numel(array_ndays)
    ndays = array_ndays(indays);
    tic;
    for i1 = 1: size(cell_results_rtpd, 1)
        for i2 = 1: size(cell_results_rtpd, 2)
            for i3 = 1: size(cell_results_rtpd, 3)
                T_results_rtpd = cell_results_rtpd{i1, i2, i3};
%                 T_results_rtpd.FRU_NEED = max(T_results_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, [], 2);
%                 T_results_rtpd.FRD_NEED = min(T_results_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, [], 2);
%                 T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_results_rtpd{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
%                 T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_results_rtpd{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
%                 T_results_rtpd = [T_results_rtpd T_fruerror_rtpd T_frderror_rtpd];
%                 T_results_rtpd.FRU_error = T_results_rtpd.FRU - T_results_rtpd.FRU_NEED;
%                 T_results_rtpd.FRD_error = T_results_rtpd.FRD - T_results_rtpd.FRD_NEED;
% 
%                 T_results_rtpd.fru_over = abs(max(0, T_results_rtpd.FRU_error));
%                 T_results_rtpd.frd_over = abs(min(0, T_results_rtpd.FRD_error));
%                 T_results_rtpd.fru_short = abs(min(0, T_results_rtpd.FRU_error));
%                 T_results_rtpd.frd_short = abs(max(0, T_results_rtpd.FRD_error));
%                 T_results_rtpd.fru_imbalance = T_results_rtpd.fru_over + T_results_rtpd.fru_short;
%                 T_results_rtpd.frd_imbalance = T_results_rtpd.frd_over + T_results_rtpd.frd_short;
%                 T_results_rtpd.fru_freqshort_hd = mean(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2);
%                 T_results_rtpd.frd_freqshort_hd = mean(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2);

                % Evaluation metrics of last 30 days, place holders
                T_results_rtpd.fru_freqshort_hd_30day = nan(size(T_results_rtpd, 1), 1); % Place holder, frequency of FRU shortage
                T_results_rtpd.frd_freqshort_hd_30day = nan(size(T_results_rtpd, 1), 1); % Place holder, frequency of FRD shortage
                T_results_rtpd.fru_short_30day = nan(size(T_results_rtpd, 1), 1); % Place holder, FRU shortage
                T_results_rtpd.frd_short_30day = nan(size(T_results_rtpd, 1), 1); % Place holder, FRD shortage
                T_results_rtpd.fru_over_30day = nan(size(T_results_rtpd, 1), 1); % Place holder, FRU oversupply
                T_results_rtpd.frd_over_30day = nan(size(T_results_rtpd, 1), 1); % Place holder, FRD oversupply
                T_results_rtpd.fru_imbalance_30day = nan(size(T_results_rtpd, 1), 1); % Place holder, total imbalance
                T_results_rtpd.frd_imbalance_30day = nan(size(T_results_rtpd, 1), 1); % Place holder, total imbalance

                % History training data
                i_start = find(T_results_rtpd.HOUR_START.Month==2, true, 'first');
                for irow = i_start:size(T_results_rtpd, 1)
                    this_hour = T_results_rtpd.HOUR_START(irow);
                    selected = (T_results_rtpd.HOUR_START<this_hour)&(T_results_rtpd.HOUR_START>=this_hour-days(ndays))&(T_results_rtpd.HOUR_START.Hour==this_hour.Hour);

                    T_results_rtpd.fru_freqshort_hd_30day(irow) = mean(T_results_rtpd{selected, 'fru_freqshort_hd'});
                    T_results_rtpd.frd_freqshort_hd_30day(irow) = mean(T_results_rtpd{selected, 'frd_freqshort_hd'});
                    T_results_rtpd.fru_short_30day(irow) = sum(T_results_rtpd{selected, 'fru_short'});
                    T_results_rtpd.frd_short_30day(irow) = sum(T_results_rtpd{selected, 'frd_short'});
                    T_results_rtpd.fru_over_30day(irow) = sum(T_results_rtpd{selected, 'fru_over'});
                    T_results_rtpd.frd_over_30day(irow) = sum(T_results_rtpd{selected, 'frd_over'});
                    T_results_rtpd.fru_imbalance_30day(irow) = sum(T_results_rtpd{selected, 'fru_imbalance'});
                    T_results_rtpd.frd_imbalance_30day(irow) = sum(T_results_rtpd{selected, 'frd_imbalance'});
                end
                cell_results_rtpd{i1, i2, i3} = T_results_rtpd;
                fprintf('i1: %g, i2: %g, i3: %g, time: %.2f s\n', i1, i2, i3, toc);
            end
        end
    end

    %% Determine optimal FRU and FRD
    nK = size(cell_results_rtpd, 1);
    nS = size(cell_results_rtpd, 2);
    nC = size(cell_results_rtpd, 3);
    T_optimal_rtpd = T_baseline_rtpd(T_baseline_rtpd.HOUR_START.Month==2, 'HOUR_START');
    T_optimal_rtpd.FRU = nan(size(T_optimal_rtpd, 1), 1);
    T_optimal_rtpd.FRD = nan(size(T_optimal_rtpd, 1), 1);
    T_optimal_rtpd.k_optimal_fru = nan(size(T_optimal_rtpd, 1), 1);
    T_optimal_rtpd.c_optimal_fru = nan(size(T_optimal_rtpd, 1), 1);
    T_optimal_rtpd.k_optimal_frd = nan(size(T_optimal_rtpd, 1), 1);
    T_optimal_rtpd.c_optimal_frd = nan(size(T_optimal_rtpd, 1), 1);

    which_criteria = 1;

    for s = 1: ns
        for irow = 1: size(T_optimal_rtpd, 1)
            this_hour = T_optimal_rtpd.HOUR_START(irow);
            criteria_fru = nan(size(cell_results_rtpd, 1)*size(cell_results_rtpd, 3), 5);
            criteria_frd = nan(size(cell_results_rtpd, 1)*size(cell_results_rtpd, 3), 5);
            criteria_fru_30day = nan(size(cell_results_rtpd, 1)*size(cell_results_rtpd, 3), 5);
            criteria_frd_30day = nan(size(cell_results_rtpd, 1)*size(cell_results_rtpd, 3), 5);
            for ik = 1:nK
                for ic = 1: nC
                    T_results_rtpd = cell_results_rtpd{ik, s, ic};
                    criteria_fru(ic+(ik-1)*nC, 1) = ik;
                    criteria_fru(ic+(ik-1)*nC, 2) = ic;
                    criteria_fru(ic+(ik-1)*nC, 3) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'fru_freqshort_hd'};
                    criteria_fru(ic+(ik-1)*nC, 4) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'fru_short'};
                    criteria_fru(ic+(ik-1)*nC, 5) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'fru_over'};
                    criteria_fru(ic+(ik-1)*nC, 6) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'fru_imbalance'};

                    criteria_frd(ic+(ik-1)*nC, 1) = ik;
                    criteria_frd(ic+(ik-1)*nC, 2) = ic;
                    criteria_frd(ic+(ik-1)*nC, 3) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'frd_freqshort_hd'};
                    criteria_frd(ic+(ik-1)*nC, 4) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'frd_short'};
                    criteria_frd(ic+(ik-1)*nC, 5) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'frd_over'};
                    criteria_frd(ic+(ik-1)*nC, 6) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'frd_imbalance'};

                    criteria_fru_30day(ic+(ik-1)*nC, 1) = ik;
                    criteria_fru_30day(ic+(ik-1)*nC, 2) = ic;
                    criteria_fru_30day(ic+(ik-1)*nC, 3) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'fru_freqshort_hd_30day'};
                    criteria_fru_30day(ic+(ik-1)*nC, 4) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'fru_short_30day'};
                    criteria_fru_30day(ic+(ik-1)*nC, 5) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'fru_over_30day'};
                    criteria_fru_30day(ic+(ik-1)*nC, 6) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'fru_imbalance_30day'};

                    criteria_frd_30day(ic+(ik-1)*nC, 1) = ik;
                    criteria_frd_30day(ic+(ik-1)*nC, 2) = ic;
                    criteria_frd_30day(ic+(ik-1)*nC, 3) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'frd_freqshort_hd_30day'};
                    criteria_frd_30day(ic+(ik-1)*nC, 4) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'frd_short_30day'};
                    criteria_frd_30day(ic+(ik-1)*nC, 5) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'frd_over_30day'};
                    criteria_frd_30day(ic+(ik-1)*nC, 6) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'frd_imbalance_30day'};
                end
            end
            switch which_criteria
                case 1 % Use 30 days data, Sort by: (reliability, oversupply, imbalance, shortage, number of days)
                    criteria_fru_sorted = sortrows(criteria_fru_30day, [3, 5, 6, 5, 1]); 
                    criteria_frd_sorted = sortrows(criteria_frd_30day, [3, 5, 6, 4, 1]);
                case 2 % Use 30 days data, Sort by: (imbalance, reliability, total shortage, total oversupply, number of days)
                    criteria_fru_sorted = sortrows(criteria_fru_30day, [6, 3, 4, 5, 1]);
                    criteria_frd_sorted = sortrows(criteria_frd_30day, [6, 3, 4, 5, 1]);
                case 3 % Use real data, Sort by: (reliability, oversupply, imbalance, shortage, number of days)
                    criteria_fru_sorted = sortrows(criteria_fru, [3, 5, 6, 4, 1]); 
                    criteria_frd_sorted = sortrows(criteria_frd, [3, 5, 6, 4, 1]);
                case 4 % Use real data, Sort by: (imbalance, reliability, total shortage, total oversupply, number of days)
                    criteria_fru_sorted = sortrows(criteria_fru, [6, 3, 4, 5, 1]); 
                    criteria_frd_sorted = sortrows(criteria_frd, [6, 3, 4, 5, 1]);
            end

            k_optimal = criteria_fru_sorted(1, 1);
            c_optimal = criteria_fru_sorted(1, 2);
            T_results_rtpd = cell_results_rtpd{k_optimal, s, c_optimal};
            T_optimal_rtpd.FRU(irow) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'FRU'};
            T_optimal_rtpd.k_optimal_fru(irow) = k_optimal;
            T_optimal_rtpd.c_optimal_fru(irow) = c_optimal;

            k_optimal = criteria_frd_sorted(1, 1);
            c_optimal = criteria_frd_sorted(1, 2);
            T_results_rtpd = cell_results_rtpd{k_optimal, s, c_optimal};
            T_optimal_rtpd.FRD(irow) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'FRD'};
            T_optimal_rtpd.k_optimal_frd(irow) = k_optimal;
            T_optimal_rtpd.c_optimal_frd(irow) = c_optimal;
        end
        cell_optimalknn_rtpd{indays, s} = T_optimal_rtpd;
    end
end

%% Plot
figure();
hold on;
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
    plot(frp_rtpd_freqshort_dknn_hd, frp_rtpd_over_dknn./1E3, '-^');
end

%% Calculate evaluation metric for the optimal kNN
% Move actual FRU/FRD needs to the resulting table
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
frp_rtpd_over_knn = fru_rtpd_over_knn + frd_rtpd_over_knn;

% 
tmp = T_optimal_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
fru_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
tmp = T_optimal_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
frd_rtpd_freqshort_knn_hd = sum(tmp(:))/numel(tmp);
frp_rtpd_freqshort_knn_hd = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
% Post-processing kNN training data
% First, concatenation results from multiple months
knn_1 = load('knn_rtpd_puresolar_1.2dim.mat');
knn_2 = load('knn_rtpd_puresolar_2.2dim.mat');

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
            cell_results_rtpd{i1, i2, i3} = [T_results_rtpd T_errormax_rtpd T_errormin_rtpd];
        end
    end
end
%%
% Calculate the evaluation metrics for the optimal selection of kNN
% parameters
for i1 = 1:size(cell_results_rtpd, 1)
    for i2 = 1:size(cell_results_rtpd, 2)
        for i3 = 1: size(cell_results_rtpd, 3)
            T_results_rtpd = cell_results_rtpd{i1, i2, i3};
            T_results_rtpd.fru_freqshort_30day_hd = nan(size(T_results_rtpd, 1), 1); % Place holder, frequency of FRU shortage
            T_results_rtpd.frd_freqshort_30day_hd = nan(size(T_results_rtpd, 1), 1); % Place holder, frequency of FRD shortage
            T_results_rtpd.fru_short_30day_hd = nan(size(T_results_rtpd, 1), 1); % Place holder, FRU shortage
            T_results_rtpd.frd_short_30day_hd = nan(size(T_results_rtpd, 1), 1); % Place holder, FRD shortage
            T_results_rtpd.fru_over_30day = nan(size(T_results_rtpd, 1), 1); % Place holder, FRU oversupply
            T_results_rtpd.frd_over_30day = nan(size(T_results_rtpd, 1), 1); % Place holder, FRD oversupply
            i_start = find(T_results_rtpd.HOUR_START.Month==2, true, 'first');
            for irow = i_start:size(T_results_rtpd, 1)
                this_hour = T_results_rtpd.HOUR_START(irow);
                selected = (T_results_rtpd.HOUR_START<this_hour)&(T_results_rtpd.HOUR_START>=this_hour-days(30))&(T_results_rtpd.HOUR_START.Hour==this_hour.Hour);
                diff_fru = T_results_rtpd{selected, 'FRU'} - T_results_rtpd{selected, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}};
                diff_frd = T_results_rtpd{selected, 'FRD'} - T_results_rtpd{selected, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}};
                
                % Frequency of shortage, HD
                logical_fru_short = diff_fru<0;
                logical_frd_short = diff_frd>0;
                T_results_rtpd.fru_freqshort_30day_hd(irow) = sum(logical_fru_short(:))/numel(logical_fru_short);
                T_results_rtpd.frd_freqshort_30day_hd(irow) = sum(logical_frd_short(:))/numel(logical_frd_short);
                
                % Total shortage, HD
                T_results_rtpd.fru_short_30day_hd(irow) = abs(sum(diff_fru(diff_fru<0)));
                T_results_rtpd.frd_short_30day_hd(irow) = abs(sum(diff_frd(diff_frd>0)));
                
                % Total oversupply
                fru_error = max(diff_fru, [], 2);
                T_results_rtpd.fru_over_30day(irow) = abs(sum(fru_error(fru_error>0)));
                frd_error = min(diff_frd, [], 2);
                T_results_rtpd.frd_over_30day(irow) = abs(sum(frd_error(frd_error<0)));
            end
            cell_results_rtpd{i1, i2, i3} = T_results_rtpd;
        end
    end
end

%% Determine optimal FRU and FRD
s = 2;
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

for irow = 1: size(T_optimal_rtpd, 1)
    this_hour = T_optimal_rtpd.HOUR_START(irow);
    criteria_fru = nan(size(cell_results_rtpd, 1)*size(cell_results_rtpd, 3), 5);
    criteria_frd = nan(size(cell_results_rtpd, 1)*size(cell_results_rtpd, 3), 5);
    for ik = 1:nK
        for ic = 1: nC
            T_results_rtpd = cell_results_rtpd{ik, s, ic};
            criteria_fru(ic+(ik-1)*nC, 1) = ik;
            criteria_fru(ic+(ik-1)*nC, 2) = ic;
            criteria_fru(ic+(ik-1)*nC, 3) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'fru_freqshort_30day_hd'};
            criteria_fru(ic+(ik-1)*nC, 4) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'fru_short_30day_hd'};
            criteria_fru(ic+(ik-1)*nC, 5) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'fru_over_30day'};

            criteria_frd(ic+(ik-1)*nC, 1) = ik;
            criteria_frd(ic+(ik-1)*nC, 2) = ic;
            criteria_frd(ic+(ik-1)*nC, 3) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'frd_freqshort_30day_hd'};
            criteria_frd(ic+(ik-1)*nC, 4) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'frd_short_30day_hd'};
            criteria_frd(ic+(ik-1)*nC, 5) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'frd_over_30day'};
        end
    end
    criteria_fru_sorted = sortrows(criteria_fru, [3, 4, 5, 1]); % Sort by: First reliability, next total shortage, next total oversupply, next number of days
    k_optimal = criteria_fru_sorted(1, 1);
    c_optimal = criteria_fru_sorted(1, 2);
    T_results_rtpd = cell_results_rtpd{k_optimal, s, c_optimal};
    T_optimal_rtpd.FRU(irow) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'FRU'};
    T_optimal_rtpd.k_optimal_fru(irow) = k_optimal;
    T_optimal_rtpd.c_optimal_fru(irow) = c_optimal;

    criteria_frd_sorted = sortrows(criteria_frd, [3, 4, 5, 1]); % Sort by: First reliability, next total shortage, next total oversupply, next number of days
    k_optimal = criteria_frd_sorted(1, 1);
    c_optimal = criteria_frd_sorted(1, 2);
    T_results_rtpd = cell_results_rtpd{k_optimal, s, c_optimal};
    T_optimal_rtpd.FRD(irow) = T_results_rtpd{T_results_rtpd.HOUR_START==T_optimal_rtpd.HOUR_START(irow), 'FRD'};
    T_optimal_rtpd.k_optimal_frd(irow) = k_optimal;
    T_optimal_rtpd.c_optimal_frd(irow) = c_optimal;
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
fru_rtpd_freqshort_baseline_hd = sum(tmp(:))/numel(tmp);
tmp = T_optimal_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
frd_rtpd_freqshort_baseline_hd = sum(tmp(:))/numel(tmp);
frp_rtpd_freqshort_baseline_hd = fru_rtpd_freqshort_baseline_hd + frd_rtpd_freqshort_baseline_hd;
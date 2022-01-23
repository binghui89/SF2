%% Load Cong's data
if ~exist('this_month')
    this_month = 2; 
end

load(strcat('knn_rtpd_puresolar_', int2str(this_month), '.mat'), 'cell_baseline_rtpd', 'karray');
% load('knn_rtpd_puresolar_2.mat', 'T_rtpd', 'T_baseline_rtpd', 'karray', 'cell_baseline_rtpd', 'cell_results_rtpd');

T_baseline_rtpd = cell_baseline_rtpd{karray == 30, 1};
dirhome = pwd;

% ML files
dirwork = '.\Cong_results\ML_Longdata';
cd(dirwork);
files_cong = dir('Predictions-1Feature*.csv');
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
dirwork = '.\Cong_results\CNN_Longdata';
cd(dirwork);
files_cnn = dir('Predictions_1Feature*');
files_abspath = cell(size(files_cnn, 1), 1);
files_dim = nan(size(files_cnn, 1), 1);
files_isup = nan(size(files_cnn, 1), 1);
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

%% Add PCA results, there is only 1 so we just add it manually
T_cnnresults = [T_cnnresults; {'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\Cong_results\CNN_Longdata\Predictions_CNN_3PC_FRD.csv', 0, 3}];
T_cnnresults = [T_cnnresults; {'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\Cong_results\CNN_Longdata\Predictions_CNN_3PC_FRU.csv', 1, 3}];

T_dlresults  = [T_dlresults; {'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\Cong_results\ML_Longdata\Predictions-3PC_PC1_PC2_PC3_FRD.csv', 0, 3}];
T_dlresults  = [T_dlresults; {'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\Cong_results\ML_Longdata\Predictions-3PC_PC1_PC2_PC3_FRU.csv', 1, 3}];

% LSTM and GRU files
T_lstmresults = table(...
    {'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\Cong_results\LSTM_Longdata\Predictions_GRU_3PC_FRD.csv'; ...
    'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\Cong_results\LSTM_Longdata\Predictions_GRU_3PC_FRU.csv'; ...
    'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\Cong_results\LSTM_Longdata\Predictions_LSTM_3PC_FRD.csv'; ...
    'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\Cong_results\LSTM_Longdata\Predictions_LSTM_3PC_FRU.csv'}, ...
    [0; 1; 0; 1], ...
    [3; 3; 3; 3], ...
    'VariableNames', ...
    {'ABSPATH', 'ISUP', 'DIM'} ...
    );

%% Load results from Cong: 1-site, 1-d classifier
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

% Load results from Cong: LSTM
cell_lstm5site_rtpd = cell(height(T_lstmresults), 1); 
for i = 1: height(T_lstmresults)
    tmp = readtable(T_lstmresults{i, 'ABSPATH'}{1});
    tmp.HOUR_START = datetime(tmp.HOUR_START, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    tmp.HOUR_START = tmp.HOUR_START - duration(7, 0, 0);
    cell_lstm5site_rtpd{i} = tmp;
    fprintf('Size: %2g, %2g, File: %s\n', size(tmp, 1), size(tmp, 2), T_lstmresults{i, 'ABSPATH'}{1});
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

% Post-process: fill Cong's missing values with baseline values: LSTM
cell_lstm5site_rtpd1 = cell(size(cell_lstm5site_rtpd));
for i = 1: numel(cell_lstm5site_rtpd)
    if ~isempty(cell_lstm5site_rtpd{i})
        if T_lstmresults{i, 'ISUP'} == 1
            column_baseline = 'FRU'; % Use CAISO baseline if DL has no results
        else
            column_baseline = 'FRD'; % Use CAISO baseline if DL has not results
        end
        tmp = T_baseline_rtpd(:, {'HOUR_START', column_baseline});

        tmp{:, 'LSTM'} = tmp{:, column_baseline};
        [Lia, Locb] = ismember(cell_lstm5site_rtpd{i}.HOUR_START, tmp.HOUR_START);
        tmp{Locb(Lia), 'LSTM'} = cell_lstm5site_rtpd{i}{Lia, 'Predicted'};
        cell_lstm5site_rtpd1{i} = tmp;
    end
end


%% Calculate evaluating metrics for Cong: 1D
load(strcat('knn_rtpd_puresolar_', int2str(this_month), '.mat'), 'T_rtpd', 'T_baseline_rtpd', 'karray', 'cell_baseline_rtpd', 'cell_results_rtpd');

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
allfeature_sites = {};
for s = 1: 5
    for i = 1: numel(allfeatures)
        allfeature_sites{i+(s-1)*numel(allfeatures)} = strcat(allfeatures{i}, '-s', int2str(s));
        % ML methods
        T_cong_d = cell_dl5site_rtpd1{contains(T_dlresults.ABSPATH, strcat('Predictions-1Feature_', allfeatures{i}, '_FRD_Site_', int2str(s), '.csv'))};
        T_cong_u = cell_dl5site_rtpd1{contains(T_dlresults.ABSPATH, strcat('Predictions-1Feature_', allfeatures{i}, '_FRU_Site_', int2str(s), '.csv'))};
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
            fru_rtpd_freqshort_knn_hd(i+(s-1)*numel(allfeatures), j) = sum(tmp(:))/numel(tmp);
            tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
            frd_rtpd_freqshort_knn_hd(i+(s-1)*numel(allfeatures), j) = sum(tmp(:))/numel(tmp);

            fru_rtpd_over_knn(i+(s-1)*numel(allfeatures), j)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
            frd_rtpd_over_knn(i+(s-1)*numel(allfeatures), j)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
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
        fru_rtpd_freqshort_knn_hd(i+(s-1)*numel(allfeatures), numel(allmethods)+1) = sum(tmp(:))/numel(tmp);
        tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
        frd_rtpd_freqshort_knn_hd(i+(s-1)*numel(allfeatures), numel(allmethods)+1) = sum(tmp(:))/numel(tmp);

        fru_rtpd_over_knn(i+(s-1)*numel(allfeatures), numel(allmethods)+1)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
        frd_rtpd_over_knn(i+(s-1)*numel(allfeatures), numel(allmethods)+1)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));

    end
end

% Now calculate evaluating metrics for Cong: PCA
allfeature_sites{numel(allfeature_sites)+1} = 'PCA';
% ML methods
T_cong_d = cell_dl5site_rtpd1{contains(T_dlresults.ABSPATH, 'Predictions-3PC_PC1_PC2_PC3_FRD.csv')};
T_cong_u = cell_dl5site_rtpd1{contains(T_dlresults.ABSPATH, 'Predictions-3PC_PC1_PC2_PC3_FRU.csv')};
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
    fru_rtpd_freqshort_knn_hd(numel(allfeature_sites), j) = sum(tmp(:))/numel(tmp);
    tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
    frd_rtpd_freqshort_knn_hd(numel(allfeature_sites), j) = sum(tmp(:))/numel(tmp);

    fru_rtpd_over_knn(numel(allfeature_sites), j)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
    frd_rtpd_over_knn(numel(allfeature_sites), j)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
end

% Now calculate evaluating metrics for Cong: PCA

% CNN
T_cong_d = cell_cnn5site_rtpd1{contains(T_cnnresults.ABSPATH, 'Predictions_CNN_3PC_FRD.csv')};
T_cong_u = cell_cnn5site_rtpd1{contains(T_cnnresults.ABSPATH, 'Predictions_CNN_3PC_FRU.csv')};
T_results_rtpd = table(T_cong_u{:, 'CNN'}, T_cong_d{:, 'CNN'}, 'VariableNames', {'FRU', 'FRD'});

T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
fru_rtpd_freqshort_knn_hd(numel(allfeature_sites), numel(allmethods)+1) = sum(tmp(:))/numel(tmp);
tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
frd_rtpd_freqshort_knn_hd(numel(allfeature_sites), numel(allmethods)+1) = sum(tmp(:))/numel(tmp);

fru_rtpd_over_knn(numel(allfeature_sites), numel(allmethods)+1)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
frd_rtpd_over_knn(numel(allfeature_sites), numel(allmethods)+1)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));

% LSTM
fru_rtpd_freqshort_knn_hd(:, numel(allmethods)+2) = nan;
frd_rtpd_freqshort_knn_hd(:, numel(allmethods)+2) = nan;
fru_rtpd_over_knn(:, numel(allmethods)+2) = nan;
frd_rtpd_over_knn(:, numel(allmethods)+2) = nan;

T_cong_d = cell_lstm5site_rtpd1{contains(T_lstmresults.ABSPATH, 'Predictions_LSTM_3PC_FRD.csv')};
T_cong_u = cell_lstm5site_rtpd1{contains(T_lstmresults.ABSPATH, 'Predictions_LSTM_3PC_FRU.csv')};
T_results_rtpd = table(T_cong_u{:, 'LSTM'}, T_cong_d{:, 'LSTM'}, 'VariableNames', {'FRU', 'FRD'});

T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
fru_rtpd_freqshort_knn_hd(numel(allfeature_sites), numel(allmethods)+2) = sum(tmp(:))/numel(tmp);
tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
frd_rtpd_freqshort_knn_hd(numel(allfeature_sites), numel(allmethods)+2) = sum(tmp(:))/numel(tmp);

fru_rtpd_over_knn(numel(allfeature_sites), numel(allmethods)+2)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
frd_rtpd_over_knn(numel(allfeature_sites), numel(allmethods)+2)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));

% GRU
fru_rtpd_freqshort_knn_hd(:, numel(allmethods)+3) = nan;
frd_rtpd_freqshort_knn_hd(:, numel(allmethods)+3) = nan;
fru_rtpd_over_knn(:, numel(allmethods)+3) = nan;
frd_rtpd_over_knn(:, numel(allmethods)+3) = nan;

T_cong_d = cell_lstm5site_rtpd1{contains(T_lstmresults.ABSPATH, 'Predictions_GRU_3PC_FRD.csv')};
T_cong_u = cell_lstm5site_rtpd1{contains(T_lstmresults.ABSPATH, 'Predictions_GRU_3PC_FRU.csv')};
T_results_rtpd = table(T_cong_u{:, 'LSTM'}, T_cong_d{:, 'LSTM'}, 'VariableNames', {'FRU', 'FRD'});

T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;

T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;

tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
fru_rtpd_freqshort_knn_hd(numel(allfeature_sites), numel(allmethods)+3) = sum(tmp(:))/numel(tmp);
tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
frd_rtpd_freqshort_knn_hd(numel(allfeature_sites), numel(allmethods)+3) = sum(tmp(:))/numel(tmp);

fru_rtpd_over_knn(numel(allfeature_sites), numel(allmethods)+3)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
frd_rtpd_over_knn(numel(allfeature_sites), numel(allmethods)+3)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));

%
frp_rtpd_freqshort_knn_hd = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
frp_rtpd_over_knn = fru_rtpd_over_knn + frd_rtpd_over_knn;

T_frp_rtpd_freqshort_knn_hd = [array2table(allfeature_sites(:), 'VariableNames', {'Features'}) array2table(frp_rtpd_freqshort_knn_hd, 'VariableNames', [allmethods, 'CNN', 'LSTM', 'GRU'])];
T_frp_rtpd_over_knn = [array2table(allfeature_sites(:), 'VariableNames', {'Features'}) array2table(frp_rtpd_over_knn, 'VariableNames', [allmethods, 'CNN', 'LSTM', 'GRU'])];



% %% Load files: 5-site PCA
% dirwork = '.\Cong_results\\ML_Longdata';
% cd(dirwork);
% files_cong = dir('Predictions-3PC_*.csv');
% files_abspath = cell(size(files_cong, 1), 1);
% files_dim = nan(size(files_cong, 1), 1);
% files_isup = nan(size(files_cong, 1), 1);
% for i = 1: size(files_cong, 1)
%     fname = files_cong(i).name;
%     fdir  = files_cong(i).folder;
%     files_abspath{i} = fullfile(fdir, fname);
%     if contains(fname, 'FRU')
%         files_isup(i) = 1;
%     else
%         files_isup(i) = 0;
%     end
%     if contains(fname, '1PC')
%         files_dim(i) = 1;
%     elseif contains(fname, '2PC')
%         files_dim(i) = 2;
%     elseif contains(fname, '3PC')
%         files_dim(i) = 3;
% %     elseif contains(fname, '4Feature')
% %         files_dim(i) = 4;
%     end
% end
% T_dlresults_2d = array2table(files_abspath, 'VariableNames', {'ABSPATH'});
% T_dlresults_2d{:, 'ISUP'} = files_isup;
% T_dlresults_2d{:, 'DIM'} = files_dim;
% cd(dirhome);
% 
% allfeatures_2d = {'mean_k_mean_w', 'mean_k_v_k', 'mean_kpv_mean_wpv', 'mean_wpv_v_kpv'};
% 
% % CNN files
% dirwork = '.\Cong_results\CNN2';
% cd(dirwork);
% files_cnn = dir('Predictions_2Feature*');
% for i = 1:size(files_cnn, 1)
%     fname = files_cnn(i).name;
%     fdir  = files_cnn(i).folder;
%     files_abspath{i} = fullfile(fdir, fname);
%     if contains(fname, 'FRU')
%         files_isup(i) = 1;
%     else
%         files_isup(i) = 0;
%     end
%     if contains(fname, '1Feature')
%         files_dim(i) = 1;
%     elseif contains(fname, '2Feature')
%         files_dim(i) = 2;
%     elseif contains(fname, '3Feature')
%         files_dim(i) = 3;
%     elseif contains(fname, '4Feature')
%         files_dim(i) = 4;
%     end
% end
% T_cnnresults_2d = array2table(files_abspath, 'VariableNames', {'ABSPATH'});
% T_cnnresults_2d{:, 'ISUP'} = files_isup;
% T_cnnresults_2d{:, 'DIM'} = files_dim;
% cd(dirhome);
% 
% % Load results from Cong: ML 2d
% cell_dl5site_rtpd_2d = cell(height(T_dlresults_2d), 1);
% for i = 1: height(T_dlresults_2d)
%     tmp = readtable(T_dlresults_2d{i, 'ABSPATH'}{1});
%     tmp.df_merge_df_tidx_clean_which_test_.TimeZone = 'UTC';
%     tmp.Properties.VariableNames{1} = 'HOUR_START';
%     cell_dl5site_rtpd_2d{i} = tmp;
%     fprintf('Size: %2g, %2g, File: %s\n', size(tmp, 1), size(tmp, 2), T_dlresults_2d{i, 'ABSPATH'}{1});
% end
% 
% % Load results from Cong: CNN 2d
% cell_cnn5site_rtpd_2d = cell(height(T_cnnresults_2d), 1); % Time axis is missing
% for i = 1: height(T_cnnresults_2d)
%     tmp = readtable(T_cnnresults_2d{i, 'ABSPATH'}{1});
%     tmp.HOUR_START = datetime(tmp.HOUR_START, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
%     tmp.HOUR_START = tmp.HOUR_START - duration(7, 0, 0);
%     cell_cnn5site_rtpd_2d{i} = tmp;
%     fprintf('Size: %2g, %2g, File: %s\n', size(tmp, 1), size(tmp, 2), T_cnnresults_2d{i, 'ABSPATH'}{1});
% end
% 
% % Post-process: fill Cong's missing values with baseline values: ML 2d
% cell_dl5site_rtpd1_2d = cell(size(cell_dl5site_rtpd_2d));
% for i = 1: numel(cell_dl5site_rtpd_2d)
%     if ~isempty(cell_dl5site_rtpd_2d{i})
%         if T_dlresults_2d{i, 'ISUP'} == 1
%             column_baseline = 'FRU'; % Use CAISO baseline if DL has no results
%         else
%             column_baseline = 'FRD'; % Use CAISO baseline if DL has not results
%         end
%         tmp = T_baseline_rtpd(:, {'HOUR_START', column_baseline});
% 
%         for j = 1: length(allmethods)
%             tmp{:, allmethods{j}} = tmp{:, column_baseline};
%             [Lia, Locb] = ismember(cell_dl5site_rtpd_2d{i}.HOUR_START, tmp.HOUR_START);
%             tmp{Locb(Lia), allmethods{j}} = cell_dl5site_rtpd_2d{i}{Lia, allmethods{j}};
%             cell_dl5site_rtpd1_2d{i} = tmp;
%         end
%     end
% end
% 
% % Post-process: fill Cong's missing values with baseline values: CNN 2d
% cell_cnn5site_rtpd1_2d = cell(size(cell_cnn5site_rtpd_2d));
% for i = 1: numel(cell_cnn5site_rtpd_2d)
%     if ~isempty(cell_cnn5site_rtpd_2d{i})
%         if T_cnnresults_2d{i, 'ISUP'} == 1
%             column_baseline = 'FRU'; % Use CAISO baseline if DL has no results
%         else
%             column_baseline = 'FRD'; % Use CAISO baseline if DL has not results
%         end
%         tmp = T_baseline_rtpd(:, {'HOUR_START', column_baseline});
% 
%         tmp{:, 'CNN'} = tmp{:, column_baseline};
%         [Lia, Locb] = ismember(cell_cnn5site_rtpd_2d{i}.HOUR_START, tmp.HOUR_START);
%         tmp{Locb(Lia), 'CNN'} = cell_cnn5site_rtpd_2d{i}{Lia, 'Predicted'};
%         cell_cnn5site_rtpd1_2d{i} = tmp;
%     end
% end
% 
% %% Calculate evaluating metrics for Cong: 2D
% load('knn_rtpd_puresolar_2.mat', 'T_rtpd', 'T_baseline_rtpd', 'karray', 'cell_baseline_rtpd', 'cell_results_rtpd');
% 
% % This is the actual need of FRP
% f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
% fru_need_rtpd = max(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4), [], 1)';
% f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
% frd_need_rtpd = min(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4), [], 1)';
% T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
% T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});
% 
% frd_rtpd_freqshort_knn_hd = [];
% fru_rtpd_freqshort_knn_hd = [];
% fru_rtpd_over_knn = [];
% frd_rtpd_over_knn = [];
% % Calculate metrics: 1D feature
% for i = 1: numel(allfeatures_2d)
%     % ML methods
%     T_cong_d = cell_dl5site_rtpd1_2d{contains(T_dlresults_2d.ABSPATH, strcat('Predictions-2Feature_', allfeatures_2d{i}, '_FRD_MultiSite_.csv'))};
%     T_cong_u = cell_dl5site_rtpd1_2d{contains(T_dlresults_2d.ABSPATH, strcat('Predictions-2Feature_', allfeatures_2d{i}, '_FRU_MultiSite_.csv'))};
%     for j = 1: numel(allmethods)
%         T_results_rtpd = table(T_cong_u{:, allmethods{j}}, T_cong_d{:, allmethods{j}}, 'VariableNames', {'FRU', 'FRD'});
% 
%         T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
%         T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;
% 
%         T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
%         T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
%         T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
%         T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
%         T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;
% 
%         tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
%         fru_rtpd_freqshort_knn_hd(i, j) = sum(tmp(:))/numel(tmp);
%         tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
%         frd_rtpd_freqshort_knn_hd(i, j) = sum(tmp(:))/numel(tmp);
% 
%         fru_rtpd_over_knn(i, j)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
%         frd_rtpd_over_knn(i, j)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
%     end
%     
%     % CNN
%     T_cong_d = cell_cnn5site_rtpd1_2d{contains(T_cnnresults_2d.ABSPATH, strcat('Predictions_2Feature-', allfeatures_2d{i}, '-FRD.csv'))};
%     T_cong_u = cell_cnn5site_rtpd1_2d{contains(T_cnnresults_2d.ABSPATH, strcat('Predictions_2Feature-', allfeatures_2d{i}, '-FRU.csv'))};
%     T_results_rtpd = table(T_cong_u{:, 'CNN'}, T_cong_d{:, 'CNN'}, 'VariableNames', {'FRU', 'FRD'});
% 
%     T_results_rtpd.FRU_error = T_results_rtpd.FRU - fru_need_rtpd;
%     T_results_rtpd.FRD_error = T_results_rtpd.FRD - frd_need_rtpd;
% 
%     T_fruerror_rtpd = array2table(T_results_rtpd.FRU-T_errormax_rtpd{:, :}, 'VariableNames', {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'});
%     T_frderror_rtpd = array2table(T_results_rtpd.FRD-T_errormin_rtpd{:, :}, 'VariableNames', {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'});
%     T_results_rtpd = [T_results_rtpd, T_errormax_rtpd, T_errormin_rtpd, T_fruerror_rtpd, T_frderror_rtpd];
%     T_results_rtpd.FRU_short_freq_hd = sum(T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0, 2)./4;
%     T_results_rtpd.FRD_short_freq_hd = sum(T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0, 2)./4;
% 
%     tmp = T_results_rtpd{:, {'FRU_error_1', 'FRU_error_2', 'FRU_error_3', 'FRU_error_4'}}<0;
%     fru_rtpd_freqshort_knn_hd(i, numel(allmethods)+1) = sum(tmp(:))/numel(tmp);
%     tmp = T_results_rtpd{:, {'FRD_error_1', 'FRD_error_2', 'FRD_error_3', 'FRD_error_4'}}>0;
%     frd_rtpd_freqshort_knn_hd(i, numel(allmethods)+1) = sum(tmp(:))/numel(tmp);
% 
%     fru_rtpd_over_knn(i, numel(allmethods)+1)  = abs(sum(T_results_rtpd.FRU_error(T_results_rtpd.FRU_error>=0)));
%     frd_rtpd_over_knn(i, numel(allmethods)+1)  = abs(sum(T_results_rtpd.FRD_error(T_results_rtpd.FRD_error<=0)));
% 
% end
% frp_rtpd_freqshort_knn_hd = fru_rtpd_freqshort_knn_hd + frd_rtpd_freqshort_knn_hd;
% frp_rtpd_over_knn = fru_rtpd_over_knn + frd_rtpd_over_knn;
% 
% T_frp_rtpd_freqshort_knn_hd_2d = [array2table(allfeatures_2d(:), 'VariableNames', {'Features'}) array2table(frp_rtpd_freqshort_knn_hd, 'VariableNames', [allmethods, 'CNN'])];
% T_frp_rtpd_over_knn_2d = [array2table(allfeatures_2d(:), 'VariableNames', {'Features'}) array2table(frp_rtpd_over_knn, 'VariableNames', [allmethods, 'CNN'])];
% 
% %% Frontier of the ML/DL-based methods
% figure();
% h = [];
% h(1) = scatter(T_metric_pshort.baseline(:), T_metric_over.baseline(:)./1E3, 60, 'ko');
% set(h(1), 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', 'k');
% % set(h(1), 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', [.5, .5, .5], 'Alpha', .5, 'MarkerEdgeAlpha', .5, 'MarkerSize', 6);
% 
% hold on;
% % for i = 
% h(2) = scatter(T_frp_rtpd_freqshort_knn_hd{:, [allmethods, 'CNN']}(:), T_frp_rtpd_over_knn{:, [allmethods, 'CNN']}(:)./1E3, 30, '.');
% h(3) = scatter(T_frp_rtpd_freqshort_knn_hd_2d{:, [allmethods, 'CNN']}(:), T_frp_rtpd_over_knn_2d{:, [allmethods, 'CNN']}(:)./1E3, 30, '.');
% % set(h(2), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% 
% xline(T_metric_pshort.baseline(karray==30), 'Color','red', 'LineStyle','--');
% yline(T_metric_over.baseline(karray==30)./1E3, 'Color','red', 'LineStyle','--');
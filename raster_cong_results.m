%% Load Cong's results
dirhome = pwd;
dirwork = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\Cong_results\ML';
cd(dirwork);
tmp = readtable('TimeIndex.csv');
ar_datetime_hour_cong93 = tmp{:, 'x'};
ar_datetime_hour_cong93.TimeZone = 'UTC'; % 93 time slices
tmp = readtable('TimeIndex_2.csv');
ar_datetime_hour_cong77 = tmp{:, 'x'};
ar_datetime_hour_cong77.TimeZone = 'UTC'; % 77 time slices

files_cong = dir('Predictions*.csv');
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

cell_dlkcs_rtpd = cell(height(T_dlresults), 1);
for i = 1: height(T_dlresults)
    tmp = readtable(T_dlresults{i, 'ABSPATH'}{1});
    if size(tmp, 1) == numel(ar_datetime_hour_cong93)
        tmp{:, 'HOUR_START'} = ar_datetime_hour_cong93;
    elseif size(tmp, 1) == numel(ar_datetime_hour_cong77)
        tmp{:, 'HOUR_START'} = ar_datetime_hour_cong77;
    else
        fprintf('Row numbers do not match: %g, File: %s\n', size(tmp, 1), T_dlresults{i, 'ABSPATH'}{1});
        continue;
    end
    cell_dlkcs_rtpd{i} = tmp;
    fprintf('Size: %2g, %2g, File: %s\n', size(tmp, 1), size(tmp, 2), T_dlresults{i, 'ABSPATH'}{1});
end

% CNN results
dirwork = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\Cong_results\CNN';
cd(dirwork);
files_cnn = dir('Predictions*.csv');
cd(dirhome);
files_abspath = cell(size(files_cnn, 1), 1);
files_name = cell(size(files_cnn, 1), 1);
for i = 1: size(files_cnn, 1)
    fname = files_cnn(i).name;
    fdir  = files_cnn(i).folder;
    files_abspath{i} = fullfile(fdir, fname);
    files_name{i} = fname;
end
T_cnnresults = array2table([files_abspath, files_name], 'VariableNames', {'ABSPATH', 'FNAME'});

for i = 1: height(T_cnnresults)
    if contains(T_dlresults{i, 'ABSPATH'}{1}, T_cnnresults{i, 'FNAME'}{1})
        tmp = readtable(T_cnnresults{i, 'ABSPATH'}{1});
        cell_dlkcs_rtpd{i}{:, 'CNN'} = tmp{:, 'Predicted'};
    else
        fprintf('Filename mismatch!\n');
    end
end

%% Show evaluating metrics of Cong's results
T_baseline_rtpd{:, 'istrain'} = 0;
T_baseline_rtpd{(T_baseline_rtpd.HOUR_START<=trainhour_end)&(T_baseline_rtpd.HOUR_START>=trainhour_start), 'istrain'} = 1;

allmethods = {'svm_r', 'svm_l', 'svm_l2', 'svm_p', 'svm_p2', 'ann1', 'ann2', 'ann4', 'ann5', 'ann6', 'gbm1', 'gbm2', 'gbm3', 'gbm4', 'rf1', 'rf2', 'CNN'};
% Post-process: fill Cong's missing values with baseline values
cell_dlkcs_rtpd1 = cell(size(cell_dlkcs_rtpd));
for i = 1: numel(cell_dlkcs_rtpd)
    if ~isempty(cell_dlkcs_rtpd{i})
        if T_dlresults{i, 'ISUP'} == 1
            column_baseline = 'FRU_k30'; % Use CAISO baseline if DL has no results
        else
            column_baseline = 'FRD_k30'; % Use CAISO baseline if DL has not results
        end
        tmp = T_baseline_rtpd(T_baseline_rtpd.istrain==0, {'HOUR_START', column_baseline});

        for j = 1: length(allmethods)
            tmp{:, allmethods{j}} = tmp{:, column_baseline};
            [Lia, Locb] = ismember(cell_dlkcs_rtpd{i}.HOUR_START, tmp.HOUR_START);
            tmp{Locb, allmethods{j}} = cell_dlkcs_rtpd{i}{:, allmethods{j}};
            cell_dlkcs_rtpd1{i} = tmp;
        end
    end
end

% Calculate oversupply and Pshortage
tmp_abspath = {};
tmp_isup = [];
tmp_dim = [];
tmp_method = {};
tmp_o = [];
tmp_pr = [];

for i = 1: numel(cell_dlkcs_rtpd)
    tmp = cell_dlkcs_rtpd1{i};
    for j = 1: length(allmethods)
        if T_dlresults{i, 'ISUP'} == 1
            [o, pr] = calculate_oversupply_and_pshort(ar_actual_error_fru((T_actual_error.istrain==0)), tmp{:, allmethods{j}});
        else
            [o, pr] = calculate_oversupply_and_pshort(-ar_actual_error_frd((T_actual_error.istrain==0)), -tmp{:, allmethods{j}});
        end
        tmp_abspath = [tmp_abspath; T_dlresults{i, 'ABSPATH'}];
        tmp_isup = [tmp_isup; T_dlresults{i, 'ISUP'}];
        tmp_dim = [tmp_dim; T_dlresults{i, 'DIM'}];
        tmp_method = [tmp_method; allmethods{j}];
        tmp_o = [tmp_o; o];
        tmp_pr = [tmp_pr; pr];
    end
end

T_dlstats = array2table(tmp_abspath, 'VariableNames', {'ABSPATH'});
T_dlstats{:, 'ISUP'} = tmp_isup;
T_dlstats{:, 'DIM'} = tmp_dim;
T_dlstats{:, 'METHOD'} = tmp_method;
T_dlstats{:, 'O'} = tmp_o;
T_dlstats{:, 'PR'} = tmp_pr;


%% Plot
for k = karray
    [otest_baseline_fru(karray==k), ptest_baseline_fru(karray==k)] = calculate_oversupply_and_pshort(ar_actual_error_fru((T_actual_error.istrain==0), :), T_baseline_rtpd{T_baseline_rtpd.istrain==0, strcat('FRU_k', num2str(k))});
    [otest_baseline_frd(karray==k), ptest_baseline_frd(karray==k)] = calculate_oversupply_and_pshort(-ar_actual_error_frd((T_actual_error.istrain==0), :), -T_baseline_rtpd{T_baseline_rtpd.istrain==0, strcat('FRD_k', num2str(k))});
end


figure(); % UP
hold on;
for i = 1:4
    h(i) = scatter(T_dlstats{(T_dlstats.ISUP==1)&(T_dlstats.DIM==i), 'PR'}, T_dlstats{(T_dlstats.ISUP==1)&(T_dlstats.DIM==i), 'O'}./1E3, '.');
end
scatter(ptest_baseline_fru, otest_baseline_fru./1E3, 60, 'ko');
set(gca, 'YScale', 'log');
xline(ptest_baseline_fru(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(otest_baseline_fru(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
title('FRU');
legend(h, '1d' ,'2d', '3d', '4d');

figure(); % Down
hold on;
for i = 1:4
    h(i) = scatter(T_dlstats{(T_dlstats.ISUP==0)&(T_dlstats.DIM==i), 'PR'}, T_dlstats{(T_dlstats.ISUP==0)&(T_dlstats.DIM==i), 'O'}./1E3, '.');
end
scatter(ptest_baseline_frd, otest_baseline_frd./1E3, 60, 'ko');
set(gca, 'YScale', 'log');
xline(ptest_baseline_frd(karray==30), 'Color', 'r', 'LineStyle', '--');
yline(otest_baseline_frd(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
title('FRD');
legend(h, '1d' ,'2d', '3d', '4d');

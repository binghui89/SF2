% For each months, each direction (FRU or FRD)
% Baseline: 12 K values: in total 12 nMAE
% 1-d 1-site: 5 sites, 12 classifiers, 12 K values: in total 720 nMAEs
% Dynamic 1-site: 5 sites, 6 N values: in total 30 nMAE
% Dynamic PCA kNN: 6 N values: in total 6 nMAEs
% ML: 12 classifiers, 5 sites, 16 methods: in total 960 nMAEs
% PCA-ML: 16 methods: in total 16 nMAEs
% CNN: 12 classifiers:  in total 12 nMAEs (it's using 5-site mean? why no site)
% PCA-CNN: in total 1 nMAE

function visualization_rtpd_nMAE(this_month)
if ~exist('this_month')
    this_month = 2; 
end


% Baseline
load(strcat('knn_rtpd_puresolar_', int2str(this_month), '.mat'), 'T_rtpd', 'T_baseline_rtpd', 'karray', 'cell_baseline_rtpd', 'cell_results_rtpd');

% This is the actual need of FRP
f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
fru_need_rtpd = max(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4), [], 1)';
f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
frd_need_rtpd = min(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4), [], 1)';
T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});

s = 1;
for k = karray
    T_baseline_rtpd = cell_baseline_rtpd{karray==k, s};
    nmae_fru_baseline(karray==k) = nmae(fru_need_rtpd, T_baseline_rtpd.FRU, 1);
    nmae_frd_baseline(karray==k) = nmae(frd_need_rtpd, T_baseline_rtpd.FRD, 0);
end

% 1-d 1 site and dynamic (opt) 1 site
load(strcat('knn_post_rtpd_puresolar_', int2str(this_month), '_complete.mat'), 'cell_optimalknn_rtpd', 'cell_results_rtpd', 'array_ndays');
nmae_fru_1d1site = nan(size(cell_results_rtpd));
nmae_frd_1d1site = nan(size(cell_results_rtpd));
nmae_fru_opt1site = nan(size(cell_optimalknn_rtpd));
nmae_frd_opt1site = nan(size(cell_optimalknn_rtpd));

for i1 = 1: size(cell_results_rtpd, 1) 
    for i2 = 1: size(cell_results_rtpd, 2)
        for i3 = 1: size(cell_results_rtpd, 3)
            T_results_rtpd = cell_results_rtpd{i1, i2, i3};
            nmae_fru_1d1site(i1, i2, i3) = nmae(fru_need_rtpd, T_results_rtpd.FRU, 1);
            nmae_frd_1d1site(i1, i2, i3) = nmae(frd_need_rtpd, T_results_rtpd.FRD, 0);
        end
    end
end

for i1 = 1: size(cell_optimalknn_rtpd, 1) 
    for i2 = 1: size(cell_optimalknn_rtpd, 2)
        for i3 = 1: size(cell_optimalknn_rtpd, 3)
            T_optimal_rtpd = cell_optimalknn_rtpd{i1, i2, i3};
            nmae_fru_opt1site(i1, i2, i3) = nmae(fru_need_rtpd, T_optimal_rtpd.FRU, 1);
            nmae_frd_opt1site(i1, i2, i3) = nmae(frd_need_rtpd, T_optimal_rtpd.FRD, 0);
        end
    end
end

% Dynamic PCA kNN 
load(strcat('knn_post_rtpd_puresolar_', int2str(this_month), '_complete.pca.mat'), 'cell_optimalknn_rtpd', 'cell_results_rtpd', 'array_ndays');
nmae_fru_optpca = nan(size(cell_optimalknn_rtpd));
nmae_frd_optpca = nan(size(cell_optimalknn_rtpd));
for i1 = 1: size(cell_optimalknn_rtpd, 1) 
    for i2 = 1: size(cell_optimalknn_rtpd, 2)
        T_optimal_rtpd = cell_optimalknn_rtpd{i1, i2};
        nmae_fru_optpca(i1, i2) = nmae(fru_need_rtpd, T_optimal_rtpd.FRU, 1);
        nmae_frd_optpca(i1, i2) = nmae(frd_need_rtpd, T_optimal_rtpd.FRD, 0);
    end
end

load_cong_results;
% ML: 12 classifiers, 5 sites, 16 methods: in total 960 nMAEs
nmae_ml = nan(120, 16);
for i = 1: 120
    T_ml = cell_dl5site_rtpd1{i};
    for j = 1: 16
        isup = T_dlresults.ISUP(i);
        if isup
            nmae_ml(i, j) = nmae(fru_need_rtpd, T_ml{:, allmethods{j}}, isup);
        else
            nmae_ml(i, j) = nmae(frd_need_rtpd, T_ml{:, allmethods{j}}, isup);
        end
    end
end
nmae_fru_ml = nmae_ml(T_dlresults.ISUP(1:120)==1, :);
nmae_frd_ml = nmae_ml(T_dlresults.ISUP(1:120)==0, :);

% PCA-ML: 16 methods: in total 16 nMAEs
nmae_fru_mlpca = nan(1, 16);
nmae_frd_mlpca = nan(1, 16);
T_mlpca = cell_dl5site_rtpd1{121, 1};
for j = 1: 16
    nmae_frd_mlpca(1, j) = nmae(frd_need_rtpd, T_mlpca{:, allmethods{j}}, 0);
end
T_mlpca = cell_dl5site_rtpd1{122, 1};
for j = 1: 16
    nmae_fru_mlpca(1, j) = nmae(fru_need_rtpd, T_mlpca{:, allmethods{j}}, 1);
end

% CNN: 12 classifiers:  in total 12 nMAEs (it's using 5-site mean? why no site)
nmae_cnn = nan(24, 1);
for i = 1: 24
    T_cnn = cell_cnn5site_rtpd1{i};
    isup = T_cnnresults.ISUP(i);
    if isup
        nmae_cnn(i) = nmae(fru_need_rtpd, T_cnn{:, 'CNN'}, isup);
    else
        nmae_cnn(i) = nmae(frd_need_rtpd, T_cnn{:, 'CNN'}, isup);
    end
end
nmae_fru_cnn = nmae_cnn(T_cnnresults.ISUP(1:24)==1);
nmae_frd_cnn = nmae_cnn(T_cnnresults.ISUP(1:24)==0);

% PCA-CNN: in total 1 nMAE
T_cnnpca = cell_cnn5site_rtpd1{25};
nmae_frd_cnnpca = nmae(frd_need_rtpd, T_cnnpca{:, 'CNN'}, 0);
T_cnnpca = cell_cnn5site_rtpd1{26};
nmae_fru_cnnpca = nmae(fru_need_rtpd, T_cnnpca{:, 'CNN'}, 1);
end


function result = nmae(actual, predict, isup)
if isup
    actual(actual<0) = 0;
    predict(predict<0) = 0;
else
    actual(actual>0) = 0;
    predict(predict>0) = 0;
end
isnotzero = (actual~=0);
% result = mean(abs(predict(isnotzero) - actual(isnotzero))./abs(actual(isnotzero)));
result = mean(abs(predict(isnotzero) - actual(isnotzero)));
end
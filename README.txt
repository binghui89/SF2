Process CAISO data
1. Collect_CAISO_xxxx.ipynb -> Download raw CAISO OASIS files (directory CAISO_OASIS)
2. Process_CAISO_xxxx.ipynb -> Process raw CAISO data and generate for_flexiramp_summary_RTD/RTPD.csv, only format processing and no calculation
3. Peep_CAISO_xxxx.ipynb -> Read for_flexiramp_summary_RTD or RTPD.csv and generate ./tmp_excels/df_rtd.csv or df_rfpd.csv, calculate advisory renewable forecasts and net load errors using persistence forecast and show snapshots

Process of the revision in August 2021.
1. Use lab_ghi2power.m to convert raw IBM GHI into power.
2. load_ibm_5sites_carlo.m (previously: load_ibm_5sites.m)
3. load_caiso_data.m (previously the same, but it is updated. Need ./tmp_excels/df_rtd.csv and df_rtpd.csv)

Process of all kNN runs in August 2021 revision
1. 1-dim classifiers (example in lab_revision_1.m)
1.1. Run single site with all 1-dim classifiers.

    % Sample MATLAB code
    for this_month = 3:4
        this_year = 2020; 
        knn_rtpd_1site_1dim;
        close all;
        save(strcat('knn_rtpd_puresolar_', int2str(this_month), '.mat'));
        clearvars -except 'this_month';
    end

1.2. Run dynamic selection processing

    % Sample MATLAB code
    clear;
    for this_month = 4:12 % Target month for dynamic selection
        tic;
        knn_1 = load(strcat('knn_rtpd_puresolar_', int2str(this_month-1), '.mat')); 
        knn_2 = load(strcat('knn_rtpd_puresolar_', int2str(this_month), '.mat'));
        % Sometimes the elements of cell_baseline_rtpd are empty except the
        % first column, so this is a preprocessing process.
        knn_1.cell_baseline_rtpd = knn_1.cell_baseline_rtpd(:, 1);
        knn_2.cell_baseline_rtpd = knn_2.cell_baseline_rtpd(:, 1);
        knn_rtpd_post;
        save(strcat('knn_post_rtpd_puresolar_', int2str(this_month), '_complete.mat'));
        fprintf('%s written to file.\n', strcat('knn_post_rtpd_puresolar_', int2str(this_month), '_complete.mat'));
        toc;
        clearvars -except 'this_month';
        close all;
    end

2. PCA classifiers (example in lab_revision.m)
2.1. Run all PCA classifiers (1 principal component to 3 principal components)

    % Sample MATLAB code
    for this_month = 6:12
        tic;
        this_year = 2020; 
        write_flag = 0;
        knn_rtpd_5site_pca;
        close all;
        save(strcat('knn_rtpd_puresolar_', int2str(this_month), '.5site.pca.mat'));
        toc;
        clearvars -except 'this_month';
    end

2.2. Run dynamic selection

    clear;
    for this_month = 4:12 % Target month for dynamic selection
        tic;
        knn_1 = load(strcat('knn_rtpd_puresolar_', int2str(this_month-1), '.5site.pca.mat')); 
        knn_2 = load(strcat('knn_rtpd_puresolar_', int2str(this_month), '.5site.pca.mat'));
        knn_rtpd_post;
        save(strcat('knn_post_rtpd_puresolar_', int2str(this_month), '_complete.pca.mat'));
        fprintf('%s written to file.\n', strcat('knn_post_rtpd_puresolar_', int2str(this_month), '_complete.pca.mat'));
        toc;
        clearvars -except 'this_month';
        close all;
    end

Visualization
For PCA results, which include almost everything in Fig.9 of the revision:
>> this_month = 8; visualization_rtpd_pca_other_month;
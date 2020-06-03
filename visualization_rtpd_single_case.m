figure();
T_results_rtpd.TIME   = T_results_rtpd.HOUR_START + duration(1, 0, 0); % We always use the end of an interval as time stamp
T_baseline_rtpd.TIME = T_baseline_rtpd.HOUR_START + duration(1, 0, 0); % We always use the end of an interval as time stamp

T_baseline_rtpd.TIME = T_baseline_rtpd.HOUR_START + duration(1, 0, 0); % We always use the end of an interval as time stamp
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), T_rtpd.UP_RTPD, '-b');
hold on;
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), -1.*T_rtpd.DOWN_RTPD, '-b');
h(1) = stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), T_rtpd.error_max, '-r');
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), T_rtpd.error_min, '-r');
h(2) = stairs(T_results_rtpd.TIME-duration(1, 0, 0), T_results_rtpd.FRU, '-g');
stairs(T_results_rtpd.TIME-duration(1, 0, 0), T_results_rtpd.FRD, '-g');
h(3) = stairs(T_baseline_rtpd.TIME-duration(1, 0, 0), T_baseline_rtpd.FRU, '-k');
stairs(T_baseline_rtpd.TIME-duration(1, 0, 0), T_baseline_rtpd.FRD, '-k');

% Use box plot to show changes of FRP, RTD
fru_compare_rtpd = [T_baseline_rtpd{(T_baseline_rtpd.HOUR_START.Hour>=16)&(T_baseline_rtpd.HOUR_START.Hour<=24), {'FRU'}} T_results_rtpd{(T_results_rtpd.HOUR_START.Hour>=16)&(T_results_rtpd.HOUR_START.Hour<=24), {'FRU'}}];
delta_fru_percent = (fru_compare_rtpd(:, 2)-fru_compare_rtpd(:, 1))./fru_compare_rtpd(:, 1);

figure();
boxplot(reshape(delta_fru_percent, 8, size(delta_fru_percent, 1)/8)', 'Label', {'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
xlabel('Time (PST)');
ylabel('Relative change to baseline');
title('FRU');

frd_compare_rtpd = [T_baseline_rtpd{(T_baseline_rtpd.HOUR_START.Hour>=16)&(T_baseline_rtpd.HOUR_START.Hour<=24), {'FRD'}} T_results_rtpd{(T_results_rtpd.HOUR_START.Hour>=16)&(T_results_rtpd.HOUR_START.Hour<=24), {'FRD'}}];
delta_frd_percent = -(frd_compare_rtpd(:, 2)-frd_compare_rtpd(:, 1))./frd_compare_rtpd(:, 1);

figure();
% boxplot(reshape(delta_frd_percent.*100, 8, size(delta_frd_percent, 1)/8)', 'Label', {'16', '17', '18', '19', '20', '21', '22', '23'});
boxplot(reshape(delta_frd_percent, 8, size(delta_frd_percent, 1)/8)', 'Label', {'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
xlabel('Time (PST)');
ylabel('Relative change to baseline');
title('FRD');

% Frequency of reserve shortage
fru_need_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START)&(T_rtpd.HOUR_START.Hour>=16)&(T_rtpd.HOUR_START.Hour<=24), 'error_max'};
fru_need_rtpd = reshape(fru_need_rtpd, 4, size(fru_need_rtpd, 1)/4)';

fru_shortage_baseline_rtpd = fru_need_rtpd - fru_compare_rtpd(:, 1);
fru_shortage_results_rtpd  = fru_need_rtpd - fru_compare_rtpd(:, 2);
fru_n_shortage_hourly_results_rtpd = sum(reshape(sum(fru_shortage_results_rtpd>0, 2), 8, size(fru_shortage_results_rtpd, 1)/8), 2);
fru_n_not_nan_hourly_results_rtpd  = sum(reshape(~isnan(fru_shortage_results_rtpd), 8, numel(fru_shortage_results_rtpd)/8), 2);
fru_n_shortage_hourly_baseline_rtpd = sum(reshape(sum(fru_shortage_baseline_rtpd>0, 2), 8, size(fru_shortage_baseline_rtpd, 1)/8), 2);
fru_n_not_nan_hourly_baseline_rtpd  = sum(reshape(~isnan(fru_shortage_baseline_rtpd), 8, numel(fru_shortage_baseline_rtpd)/8), 2);
figure(); 
bar([fru_n_shortage_hourly_results_rtpd./fru_n_not_nan_hourly_results_rtpd fru_n_shortage_hourly_baseline_rtpd./fru_n_not_nan_hourly_baseline_rtpd]); % Percentage of violations each hour
legend({'New', 'Baseline'});
set(gca,'xticklabel',{'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
ylabel('Time of FRP shortage');
title('FRU');

frd_need_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START)&(T_rtpd.HOUR_START.Hour>=16)&(T_rtpd.HOUR_START.Hour<=24), 'error_min'};
frd_need_rtpd = reshape(frd_need_rtpd, 4, size(frd_need_rtpd, 1)/4)';

frd_shortage_baseline_rtpd = frd_need_rtpd - frd_compare_rtpd(:, 1);
frd_shortage_results_rtpd  = frd_need_rtpd - frd_compare_rtpd(:, 2);
frd_n_shortage_hourly_results_rtpd = sum(reshape(sum(frd_shortage_results_rtpd<0, 2), 8, size(frd_shortage_results_rtpd, 1)/8), 2);
frd_n_not_nan_hourly_results_rtpd  = sum(reshape(~isnan(frd_shortage_results_rtpd), 8, numel(frd_shortage_results_rtpd)/8), 2);
frd_n_shortage_hourly_baseline_rtpd = sum(reshape(sum(frd_shortage_baseline_rtpd<0, 2), 8, size(frd_shortage_baseline_rtpd, 1)/8), 2);
frd_n_not_nan_hourly_baseline_rtpd  = sum(reshape(~isnan(frd_shortage_baseline_rtpd), 8, numel(frd_shortage_baseline_rtpd)/8), 2);
figure(); 
bar([frd_n_shortage_hourly_results_rtpd./frd_n_not_nan_hourly_results_rtpd frd_n_shortage_hourly_baseline_rtpd./frd_n_not_nan_hourly_baseline_rtpd]); % Percentage of violations each hour
legend({'New', 'Baseline'});
set(gca,'xticklabel',{'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
ylabel('Time of FRP shortage');
title('FRD');

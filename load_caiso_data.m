dt_rtd = 5; % min
dt_rtpd = 15; % min

%% Load CAISO data RTD
if ispc
    T_rtd = readtable('.\tmp_excels\df_rtd.csv', 'Delimiter',',');
elseif isunix
    T_rtd = readtable('/home/bxl180002/git/SF2/tmp_excels/df_rtd.csv', 'Delimiter',',');
end
T_rtd.TIME = datetime(T_rtd.Var1, 'InputFormat', 'yyyy-MM-dd'' ''HH:mm:ssXXX', 'TimeZone', 'UTC'); % This is the time of the end of an interval
T_rtd.TIME_START = T_rtd.TIME - duration(0, dt_rtd, 0); % This is the time of the start of an interval
T_rtd.HOUR_START = datetime(T_rtd.TIME_START.Year, T_rtd.TIME_START.Month, T_rtd.TIME_START.Day, T_rtd.TIME_START.Hour, 0, 0, 'TimeZone', 'UTC');
% T_rtd.local_time = datetime(T_rtd.TIME, 'TimeZone', 'America/Los_Angeles');
T_rtd.local_time = datetime(T_rtd.TIME, 'TimeZone', '-08:00'); % We use PST as local time

T_rtd.FORECAST_ERROR_Brtd_Artd_solar = (-T_rtd.Solar_NP15_B_RTD-T_rtd.Solar_ZP26_B_RTD-T_rtd.Solar_SP15_B_RTD) - (-T_rtd.Solar_NP15_A_RTD-T_rtd.Solar_ZP26_A_RTD-T_rtd.Solar_SP15_A_RTD); % Pure solar forecasting errors

% Demonstrate FRP requirements and binding interval forecasts, RTD
% ax1 = subplot(2, 1, 1);
% plot(T_rtd.TIME, T_rtd.UP_RTD, '-b', T_rtd.TIME, -1.*T_rtd.DOWN_RTD, '-b', T_rtd.TIME, T_rtd.FORECAST_ERROR_Brtd_Artd, '-r');
% 
% ax2 = subplot(2, 1, 2);
% plot(T_rtd.TIME, sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2), 'b');
% linkaxes([ax1(1), ax2],'x');

%% Load CAISO data RTPD
if ispc
    T_rtpd = readtable('.\tmp_excels\df_rtpd.csv', 'Delimiter',',');
elseif isunix
    T_rtpd = readtable('/home/bxl180002/git/SF2/tmp_excels/df_rtpd.csv', 'Delimiter',',');
end
T_rtpd.TIME = datetime(T_rtpd.Var1, 'InputFormat', 'yyyy-MM-dd'' ''HH:mm:ssXXX', 'TimeZone', 'UTC');
T_rtpd.TIME_START = T_rtpd.TIME - duration(0, dt_rtpd, 0); % This is the time of the start of an interval
T_rtpd.HOUR_START = datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, T_rtpd.TIME_START.Hour, 0, 0, 'TimeZone', 'UTC');
T_rtpd.local_time = datetime(T_rtpd.TIME, 'TimeZone', 'America/Los_Angeles');

% calculate pure solar forecast error
tmp_rtp_solar  = sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2);
tmp_rtpd_solar = sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2);
tmp_error_solar = -reshape(tmp_rtp_solar, 3, size(tmp_rtp_solar, 1)/3)' - (-repmat(tmp_rtpd_solar, 1, 3));

% Calculate net load forecast error
% tmp_rtd_nl_b  = T_rtd.LOAD_B_RTD - sum(T_rtd{:, {'Wind_NP15_B_RTD', 'Wind_SP15_B_RTD', 'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2);
% tmp_rtpd_nl_a = T_rtpd.LOAD_B_RTPD - sum(T_rtpd{:, {'Wind_NP15_RTPD', 'Wind_SP15_RTPD', 'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2); % Note we use binding forecast as advisory forecast since CAISO does not publish advisory forecast
tmp_rtd_nl_b  = - sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2);
tmp_rtpd_nl_a = - sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2); % Note we use binding forecast as advisory forecast since CAISO does not publish advisory forecast

tmp_error_nl = reshape(tmp_rtd_nl_b, 3, size(tmp_rtd_nl_b, 1)/3)' - (repmat(tmp_rtpd_nl_a, 1, 3));

T_rtpd.error_max = max(tmp_error_nl, [], 2);
T_rtpd.error_min = min(tmp_error_nl, [], 2);

% Demonstrate FRP requirements and binding interval forecasts, RTPD
% figure();
% ax1 = subplot(2, 1, 1);
% plot(T_rtpd.TIME, T_rtpd.UP_RTPD, '-b', T_rtpd.TIME, -1.*T_rtpd.DOWN_RTPD, '-b', T_rtpd.TIME, T_rtpd.error_max, '-r', T_rtpd.TIME, T_rtpd.error_min, '-r');
% 
% ax2 = subplot(2, 1, 2);
% plot(T_rtpd.TIME, sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2), 'b');
% linkaxes([ax1(1), ax2],'x');

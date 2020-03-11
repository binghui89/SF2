dt_rtd = 5; % min
dt_rtpd = 15; % min

%% Load CAISO data RTD
T_rtd = readtable('C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\tmp_excels\df_rtd.csv', 'Delimiter',',');
T_rtd.TIME = datetime(T_rtd.Var1, 'InputFormat', 'yyyy-MM-dd'' ''HH:mm:ssXXX', 'TimeZone', 'UTC'); % This is the time of the end of an interval
T_rtd.TIME_START = T_rtd.TIME - duration(0, dt_rtd, 0); % This is the time of the start of an interval
T_rtd.HOUR_START = datetime(T_rtd.TIME_START.Year, T_rtd.TIME_START.Month, T_rtd.TIME_START.Day, T_rtd.TIME_START.Hour, 0, 0, 'TimeZone', 'UTC');
% T_rtd.local_time = datetime(T_rtd.TIME, 'TimeZone', 'America/Los_Angeles');
T_rtd.local_time = datetime(T_rtd.TIME, 'TimeZone', '-08:00');

T_rtd.FORECAST_ERROR_Brtd_Artd_solar = (-T_rtd.Solar_NP15_B_RTD-T_rtd.Solar_ZP26_B_RTD-T_rtd.Solar_SP15_B_RTD) - (-T_rtd.Solar_NP15_A_RTD-T_rtd.Solar_ZP26_A_RTD-T_rtd.Solar_SP15_A_RTD); % Pure solar forecasting errors

% Load CAISO data RTPD
T_rtpd = readtable('C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\tmp_excels\df_rtpd.csv', 'Delimiter',',');
T_rtpd.TIME = datetime(T_rtpd.Var1, 'InputFormat', 'yyyy-MM-dd'' ''HH:mm:ssXXX', 'TimeZone', 'UTC');
T_rtpd.local_time = datetime(T_rtpd.TIME, 'TimeZone', 'America/Los_Angeles');

tmp_rtp_solar  = sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2);
tmp_rtpd_solar = sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2);
tmp_error_solar = -reshape(tmp_rtp_solar, 3, size(tmp_rtp_solar, 1)/3)' - (-repmat(tmp_rtpd_solar, 1, 3));
T_rtpd.error_max_solar = max(tmp_error_solar, [], 2);
T_rtpd.error_min_solar = min(tmp_error_solar, [], 2);

% %% Demonstrate FRP requirements and binding interval forecasts, RTD
% ax1 = subplot(2, 1, 1);
% plot(T_rtd.TIME, T_rtd.UP_RTD, '-b', T_rtd.TIME, -1.*T_rtd.DOWN_RTD, '-b', T_rtd.TIME, T_rtd.FORECAST_ERROR_Brtd_Artd, '-r');
% 
% ax2 = subplot(2, 1, 2);
% plot(T_rtd.TIME, sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2), 'b');
% linkaxes([ax1(1), ax2],'x');

% %% Demonstrate FRP requirements and binding interval forecasts, RTPD
% ax1 = subplot(2, 1, 1);
% plot(T_rtpd.TIME, T_rtpd.UP_RTPD, '-b', T_rtpd.TIME, -1.*T_rtpd.DOWN_RTPD, '-b', T_rtpd.TIME, T_rtpd.error_max, '-r', T_rtpd.TIME, T_rtpd.error_min, '-r');
% 
% ax2 = subplot(2, 1, 2);
% plot(T_rtpd.TIME, sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2), 'b');
% linkaxes([ax1(1), ax2],'x');

%% Explore full reference curve: Load CAISO's PV plants 
T_eia860_caiso = readtable('eia860_2018_caiso.csv', 'Delimiter',',');
% T_caiso_plts = grpstats(T_eia860_caiso(:, {'PlantCode', 'NameplateCapacity_MW_', 'Latitude', 'Longitude', 'AzimuthAngle', 'TiltAngle', 'Single_AxisTracking_'}), {'PlantCode', 'Latitude', 'Longitude'}, @sum);
T_caiso_plts = T_eia860_caiso;
T_caiso_plts.AzimuthAngle(isnan(T_caiso_plts.AzimuthAngle))=180; % Southward
T_caiso_plts.TiltAngle(isnan(T_caiso_plts.TiltAngle))=0; % Horizontal
%% Explore full reference curve: bottom-up method
ngen = size(T_caiso_plts, 1);

tarray_formean = repmat(T_rtd.TIME(:)', dt_rtd, 1) - datenum(repmat([dt_rtd:-1:1]'./24/60, 1, numel(T_rtd.TIME))); % All minutes during the past period, can be 15 or 5 min
tarray_formean = tarray_formean(:); % One minute time array
Time_formean.UTCOffset = zeros(size(tarray_formean, 1), 1); % We always use UTC time
% Time_formean.UTCOffset = tzoffset(tarray_formean); 
Time_formean.year   = tarray_formean.Year;
Time_formean.month  = tarray_formean.Month;
Time_formean.day    = tarray_formean.Day;
Time_formean.hour   = tarray_formean.Hour;
Time_formean.minute = tarray_formean.Minute;
Time_formean.second = tarray_formean.Second;
dayofyear = pvl_date2doy(Time_formean.year, Time_formean.month, Time_formean.day);

E_cs_formean = nan(size(tarray_formean, 1), ngen); % Total incident irradiance (W/m^2)
Albedo = 0.2;
PresPa=101325;

for i = 1:ngen
    Location = pvl_makelocationstruct(T_caiso_plts.Latitude(i), T_caiso_plts.Longitude(i));
    AxisTilt = T_caiso_plts.TiltAngle(i);
    AxisAzimuth = T_caiso_plts.AzimuthAngle(i);
%     [~, ~, AppSunEl_formean, ~] = pvl_ephemeris(Time_formean,Location);
    [SunAz_formean, SunEl_formean, AppSunEl_formean, ~] = pvl_ephemeris(Time_formean,Location);
%     ghi_cs_formean(:, i) = pvl_clearsky_haurwitz(90-AppSunEl_formean);

    if strcmp(T_caiso_plts.Single_AxisTracking_(i), 'Y')
        % Fixed tilt
        MaxAngle = 90;
        [~, AOI, ~, ~] = pvl_singleaxis(90-AppSunEl_formean, SunAz_formean, Location.latitude, AxisTilt, AxisAzimuth, MaxAngle);
    else
        AOI = pvl_getaoi(AxisTilt, AxisAzimuth, 90-AppSunEl_formean, SunAz_formean);
    end
    
    ghi_cs_formean = pvl_clearsky_haurwitz(90-AppSunEl_formean);
    EdiffGround = pvl_grounddiffuse(AxisTilt, ghi_cs_formean, Albedo);
    dni_cs_formean = pvl_disc(ghi_cs_formean,90-SunEl_formean, dayofyear,PresPa);
    dhi_cs_formean = ghi_cs_formean - cosd(90-SunEl_formean).*dni_cs_formean;
    Eb = 0*AOI; %Initiallize variable
    Eb(AOI<90) = dni_cs_formean(AOI<90).*cosd(AOI(AOI<90)); %Only calculate when sun is in view of the plane of array
    EdiffSky = pvl_isotropicsky(AxisTilt,dhi_cs_formean);
    E_cs_formean(:, i) = Eb + EdiffSky + EdiffGround; % Total incident irradiance (W/m^2)
    
end

ratio = E_cs_formean;
ratio(ratio<150) = ratio(ratio<150).^2./(150*1000);
ratio(ratio>=150) = ratio(ratio>=150)./1000;
ratio(ratio>=1) = 1;
Pac_formean = ratio*diag(T_caiso_plts.NameplateCapacity_MW_);
Pcap = repmat(T_caiso_plts.NameplateCapacity_MW_', size(tarray_formean, 1), 1);
% Pac_formean(Pac_formean>Pcap) = Pcap(Pac_formean>Pcap);
Pac = reshape(mean(reshape(Pac_formean, 5, numel(Pac_formean)/5), 1)', size(T_rtd, 1), ngen);
figure();
plot(T_rtd.TIME, sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2), 'k', T_rtd.TIME, sum(Pac, 2));

%% Explore reference curve: Explore CAISO's border points

caiso_esws = [
    41.031389, -121.4225;   % Northmost, plant ID: 58814
    32.586536, -117.004884; % Southmost, Plant ID: 61843
    38.773056, -123.017778; % Westmost, Plant ID: 58949
    33.332777, -112.911388; % Eastmost, Plant ID: 60307
    ];

figure();
cap = 100;
ghi_cs = 0;
Pac = 0;
for i = 1: size(caiso_esws, 1)
    ghi_cs = ghi_cs + mean_clearsky_ghi(caiso_esws(i, 1), caiso_esws(i, 2), T_rtd.TIME, dt_rtd);
    ratio = ghi_cs;
    ratio(ratio<150) = ratio(ratio<150).^2./(150*1000);
    ratio(ratio>=150) = ratio(ratio>=150)./1000;
    ratio(ratio >= 1) = 1;
    Pac = Pac + ratio*cap;
%     yyaxis left;
%     plot(T_rtd.TIME-duration(0, dt_rtd, 0), ghi_cs, '-b');
%     hold on;
%     yyaxis right;
%     plot(T_rtd.TIME-duration(0, dt_rtd, 0), Pac, '-k');
end
% plot(T_rtd.TIME-duration(0, dt_rtd, 0), Pac, '-k');
plot(T_rtd.TIME-duration(0, dt_rtd, 0), ghi_cs, '-k');
yyaxis right;
plot(T_rtd.TIME-duration(0, dt_rtd, 0), sum(T_rtd{:, {'Solar_NP15_RTD', 'Solar_SP15_RTD', 'Solar_ZP26_RTD'}}, 2), 'r');

%% Explore reference curve: Maximum data across a month

% h = [1/2 1/2];
% binomialCoeff = conv(h,h);
% for i = 1:4
%     binomialCoeff = conv(binomialCoeff,h);
% end

T_rtd.local_time_start = T_rtd.local_time - duration(0, 5, 0);
ref_curve = [];

for m = 1:12
    figure(); 
    hold on;
    for d = 1:eomday(2019, m)
        caiso_solar_ref_rtd(:, d) = sum(T_rtd{(T_rtd.local_time_start.Year==2019)&(T_rtd.local_time_start.Day==d)&(T_rtd.local_time_start.Month==m), {'Solar_NP15_RTD', 'Solar_SP15_RTD', 'Solar_ZP26_RTD'}}, 2);
        plot(caiso_solar_ref_rtd(:, d), 'k'); 
    end
%     plot(caiso_solar_ref_rtd(:, 1), 'r');
%     plot(caiso_solar_ref_rtd(:, end), 'b');
    max_daily = max(caiso_solar_ref_rtd, [], 2);
    plot(max_daily, 'g', 'LineWidth', 2);
    
%     binomialMA = filter(binomialCoeff, 1, max_daily);
%     fDelay = (length(binomialCoeff)-1)/2;
%     plot(binomialMA(fDelay+1:end), 'r', 'LineWidth', 2);

    hoursPerDay = 11;
    coeff24hMA = ones(1, hoursPerDay)/hoursPerDay;
    avg24hTempC = filter(coeff24hMA, 1, max_daily);
    fDelay = (length(coeff24hMA)-1)/2;
    plot(avg24hTempC(fDelay+1:end), 'r', 'LineWidth', 2);

%     [yupper,ylower] = envelope(max(caiso_solar_ref_rtd, [], 2));
%     plot(yupper, 'r', 'LineWidth', 2);
%     plot(ylower, 'r', 'LineWidth', 2);
    title(m);
    ref_curve = [ref_curve repmat([avg24hTempC(fDelay+1:end); zeros(fDelay, 1)], 1, eomday(2019, m))];
end

T_rtd.REF_CURVE = nan;
T_rtd.REF_CURVE(T_rtd(T_rtd.local_time_start.Year==2019, :), 1) = ref_curve(:);

T_rtd.Solar_B_RTD = sum(T_rtd{:, {'Solar_NP15_RTD', 'Solar_SP15_RTD', 'Solar_ZP26_RTD'}}, 2);
T_rtd.Solar_A_RTD = nan(size(T_rtd, 1), 1);
T_rtd.Solar_A_RTD(2:end) = T_rtd.Solar_B_RTD(1:end-1)./T_rtd.REF_CURVE(1:end-1).*T_rtd.REF_CURVE(2:end);
T_rtd.Solar_A_RTD(T_rtd.REF_CURVE==0) = 0;
T_rtd.Wind_B_RTD = sum(T_rtd{:, {'Wind_NP15_RTD', 'Wind_SP15_RTD'}}, 2);
T_rtd.Wind_A_RTD = sum(T_rtd{:, {'Wind_NP15_A_RTD', 'Wind_SP15_A_RTD'}}, 2);
T_rtd.FORECAST_ERROR_Brtd_Artd = -(T_rtd.Wind_B_RTD + T_rtd.Solar_B_RTD) + (T_rtd.Wind_A_RTD + T_rtd.Solar_A_RTD);
T_rtd.kp = T_rtd.Solar_B_RTD./T_rtd.REF_CURVE;
% T_rtd.kp(T_rtd.REF_CURVE==0) = 0;

% Select month
this_month = 12;

% Baseline
T_baseline_rtd = array2table(unique(T_rtd.HOUR_START(T_rtd.HOUR_START.Month==this_month)), 'VariableNames', {'HOUR_START'}); % Result container
% T_rtd_DATE = datetime(T_rtd.HOUR_START.Year, T_rtd.HOUR_START.Month, T_rtd.HOUR_START.Day, 'TimeZone', 'UTC');
% T_rtd_DATE_unique = unique(T_rtd_DATE);
for i = 1: size(T_baseline_rtd, 1)
    this_date = datetime(T_baseline_rtd.HOUR_START.Year(i), T_baseline_rtd.HOUR_START.Month(i), T_baseline_rtd.HOUR_START.Day(i), 'TimeZone', 'UTC');
    this_hour = T_baseline_rtd.HOUR_START.Hour(i);
    selected_days = (T_rtd.TIME_START<this_date)&(T_rtd.TIME_START>=this_date-days(30))&(T_rtd.HOUR_START.Hour==this_hour); % We use 30 previous days
%     selected_days = (T_rtd_DATE<which_date)&(T_rtd_DATE>=which_date-days(30))&(T_rtd.HOUR_START.Hour==which_hour); % We use 30 previous days
%     valid_days = T_rtd_DATE_unique((T_rtd_DATE_unique<which_date)&(T_rtd_DATE_unique>=which_date-days(30)));
%     selected_days = (ismember(T_rtd_DATE, valid_days))&(T_rtd.HOUR_START.Hour==which_hour);
    sample_error = T_rtd{selected_days, 'FORECAST_ERROR_Brtd_Artd'};
%     T_sample = T_pwr_hourly((T_pwr_hourly.DATE<which_date)&(T_pwr_hourly.HOUR_START.Hour==which_hour), :);
%     T_sample_sorted = T_sample(T_sample.DATE>=which_date-days(30), :); % We use 30 previous days
%     selected_days = ismember(datetime(T_rtd.TIME_START.Year, T_rtd.TIME_START.Month, T_rtd.TIME_START.Day), T_sample_sorted.DATE(1:30)); % 30 the nearest days
%     sample_error = T_rtd{selected_days&(T_rtd.TIME_START.Hour==which_hour), 'FORECAST_ERROR_Brtd_Artd'};
    [f,x] = ecdf(sample_error(:));
    T_baseline_rtd.FRU(i) = interp1(f, x, 0.975);
    T_baseline_rtd.FRD(i) = interp1(f, x, 0.025);
end

figure();
T_baseline_rtd.TIME = T_baseline_rtd.HOUR_START + duration(1, 0, 0); % We always use the end of an interval as time stamp

stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.UP_RTD, '-b');
hold on;
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), -1.*T_rtd.DOWN_RTD, '-b');
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.FORECAST_ERROR_Brtd_Artd, '-r')
stairs(T_baseline_rtd.TIME-duration(1, 0, 0), T_baseline_rtd.FRU, '-k');
stairs(T_baseline_rtd.TIME-duration(1, 0, 0), T_baseline_rtd.FRD, '-k');



%% Load IBM data, forecast 15-min
month_ibm = [
    201908;
    201909;
    201910;
    201911;
    201912;
    202001;
];

uniquegen     = {'gen55', 'gen56', 'gen59', 'gen60', 'gen61'};
uniquegensite = {'MNCC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz'};

allgen  = {'gen55', 'gen56', 'gen57', 'gen58', 'gen59', 'gen60', 'gen61',    'gen62', 'gen63', 'gen64'};
gensite = {'MNCC1', 'STFC1', 'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};
genlat = [34.31,34.12,34.12,34.12,37.41,35.53,35.38,34.31,34.31,35.53];
genlon = [-117.5,-117.94, -117.94, -117.94,-119.74, -118.63, -120.18,-117.5,-117.5,-118.63];

dirhome = pwd;
cell_pwr = cell(numel(uniquegen), 1);
dirwork = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\IBM\F_pwr.201908.15min';
for i = 1: numel(uniquegen)
    gen = uniquegen{i};
    for m = 1: numel(month_ibm)
        mstr = num2str(month_ibm(m));
        dirwork = strcat('C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\IBM\F_pwr.', mstr, '.15min');
        csvname = strcat(dirwork, '\', 'frcst_', gen, '.csv');
        cd(dirwork);
        if m == 1
            T_pwr = readtable(csvname, 'Delimiter', ',');
        else
            T_pwr = vertcat(T_pwr, readtable(csvname, 'Delimiter', ','));
        end
    end
    [~, ia, ~] = unique(T_pwr{:, {'Year', 'Month', 'Day', 'Hour', 'Minute'}}, 'rows');
    T_pwr = T_pwr(ia, :);
    T_pwr.TIME = datetime(T_pwr.Year, T_pwr.Month, T_pwr.Day, T_pwr.Hour, T_pwr.Minute, 0, 'TimeZone', 'UTC');
    
    % Calculate uncertainty bandwidth and variability
    T_pwr.kt_p025 = T_pwr.ghi_p025./T_pwr.ghi_cs;
    T_pwr.kt_p050 = T_pwr.ghi_p050./T_pwr.ghi_cs;
    T_pwr.kt_p075 = T_pwr.ghi_p075./T_pwr.ghi_cs;
    T_pwr.kt_width = T_pwr.kt_p075 - T_pwr.kt_p025;
    
    cell_pwr{i} = T_pwr;
end
cd(dirhome);

%% Load USCRN data, compare with IBM's nearest site

% s_uscrn = 'Santa Barbara';
% s_uscrn = 'Yosemite';
s_uscrn = 'Bodega';
switch s_uscrn
    case 'Santa Barbara'
        ibm_site = 'SBVC1';
        txtfile = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Santa_Barbara_11_W.txt';
        csv_f = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\IBM\F_ghi.201908.15min\IBM_processed_SBVC1.csv';
        csv_a = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\IBM\A_ghi.201908\IBM_processed_SBVC1.hourly.csv';
    case 'Yosemite'
        ibm_site = 'MIAC1';
        txtfile = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Yosemite_Village_12_W.txt';
        csv_f = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\IBM\F_ghi.201908.15min\IBM_processed_MIAC1.csv';
        csv_a = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\IBM\A_ghi.201908\IBM_processed_MIAC1.hourly.csv';
    case 'Bodega'
        ibm_site = 'RSAC1';
        txtfile = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Bodega_6_WSW.txt';
        csv_f = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\IBM\F_ghi.201908.15min\IBM_processed_RSAC1.csv';
        csv_a = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\IBM\A_ghi.201908\IBM_processed_RSAC1.hourly.csv';
end

T_uscrn = readtable(txtfile, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
T_uscrn.Properties.VariableNames = {'WBANNO', 'UTC_DATE', 'UTC_TIME', 'LST_DATE', 'LST_TIME', 'CRX_VN', 'LONGITUDE', 'LATITUDE', 'AIR_TEMPERATURE', 'PRECIPITATION', 'SOLAR_RADIATION', 'SR_FLAG', 'SURFACE_TEMPERATURE', 'ST_TYPE', 'ST_FLAG', 'RELATIVE_HUMIDITY', 'RH_FLAG', 'SOIL_MOISTURE_5', 'SOIL_TEMPERATURE_5', 'WETNESS', 'WET_FLAG', 'WIND_1_5', 'WIND_FLAG'};
T_uscrn.TIME = datetime(T_uscrn.UTC_DATE, 'ConvertFrom', 'yyyymmdd') + timeofday(datetime(num2str(T_uscrn.UTC_TIME, '%04d'), 'Format', 'HHmm'));
T_uscrn.TIME.TimeZone = 'UTC';

T_frcst = readtable(csv_f, 'Delimiter', ',');
T_frcst.TIME = datetime(T_frcst.Year, T_frcst.Month, T_frcst.Day, T_frcst.Hour, T_frcst.Minute, 0, 'TimeZone', 'UTC');

T_actual = readtable(csv_a, 'Delimiter', ',');
T_actual.TIME = datetime(T_actual.Year, T_actual.Month, T_actual.Day, T_actual.Hour, T_actual.Minute, 0, 'TimeZone', 'UTC');

ax1 = subplot(3, 1, 1);
stairs(T_frcst.TIME-duration(0, dt_rtpd, 0), T_frcst.p50, 'k');
hold on;
stairs(T_frcst.TIME-duration(0, dt_rtpd, 0), T_frcst.p25, 'b');
stairs(T_frcst.TIME-duration(0, dt_rtpd, 0), T_frcst.p75, 'b');
stairs(T_frcst.TIME-duration(0, dt_rtpd, 0), T_frcst.p5, 'g');
stairs(T_frcst.TIME-duration(0, dt_rtpd, 0), T_frcst.p95, 'g');
title(strcat('p-Forecast:', ibm_site));

ax2 = subplot(3, 1, 2);
stairs(T_actual.TIME-duration(1, 0, 0), T_actual.global_solar_irradiance, 'k');
hold on;
stairs(T_uscrn.TIME-duration(0, dt_rtd, 0), T_uscrn.SOLAR_RADIATION, 'b');
title(strcat('RAWS-', ibm_site, ' vs USCRN-', s_uscrn));
ylim([0, 1300]);

ax3 = subplot(3, 1, 3);
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2), 'b');
title('CAISO RTD Binding solar');


% linkaxes([ax1, ax2], 'x');
linkaxes([ax1, ax2, ax3], 'x');


%% Figures, clear-sky index, RTD 1
figure();
ax_ibm = nan(1, numel(cell_pwr));
for i = 1:numel(cell_pwr)
    ax_ibm(i) = subplot(3, 3, i);
    T_pwr = cell_pwr{i};
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p050./T_pwr.ghi_cs, 'k');
    hold on;
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p025./T_pwr.ghi_cs, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p075./T_pwr.ghi_cs, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p005./T_pwr.ghi_cs, 'g');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p095./T_pwr.ghi_cs, 'g');
    ylim([0, 1.2]);
    plot_title = strcat(uniquegen{i}, '-', uniquegensite{i});
    title(plot_title);
end

ax1 = subplot(3, 3, numel(cell_pwr)+1);
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.UP_RTD, '-b');
hold on;
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), -1.*T_rtd.DOWN_RTD, '-b');
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.FORECAST_ERROR_Brtd_Artd_solar, '-r')

ax2 = subplot(3, 3, numel(cell_pwr)+2);
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2), 'b');
title('RTD Binding solar');

linkaxes([ax1, ax2, ax_ibm],'x');
% suptitle('Clear-sky index');

%% Figure, GHI, RTD 2
figure();
ax_ibm = nan(1, numel(cell_pwr));
for i = 1:numel(cell_pwr)
    ax_ibm(i) = subplot(3, 3, i);
    T_pwr = cell_pwr{i};
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p050, 'k');
    hold on;
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p025, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p075, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p005, 'g');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p095, 'g');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_cs, 'r');
    plot_title = strcat(uniquegen{i}, '-', uniquegensite{i});
    title(plot_title);
end

ax1 = subplot(3, 3, numel(cell_pwr)+1);
% ax1 = subplot(2, 1, 1);
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.UP_RTD, '-b');
hold on;
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), -1.*T_rtd.DOWN_RTD, '-b');
% stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.FORECAST_ERROR_Brtd_Artd_solar, '-r')
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.FORECAST_ERROR_Brtd_Artd, '-r')
title('CAISO FRP');
ylabel('MW');

ax2 = subplot(3, 3, numel(cell_pwr)+2);
% ax2 = subplot(2, 1, 2);
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2), 'b');
title('CAISO RTD solar binding frcst');
ylabel('MW');

linkaxes([ax1, ax2, ax_ibm],'x');
% linkaxes([ax1, ax2],'x');
% suptitle('GHI');

%% Figure, normalized power by clear-sky power, RTD 3
figure();
ax_ibm = nan(1, numel(cell_pwr));
for i = 1:numel(cell_pwr)
    ax_ibm(i) = subplot(3, 3, i);
    T_pwr = cell_pwr{i};
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.pwr_p050./T_pwr.pwr_cs, 'k');
    hold on;
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.pwr_p025./T_pwr.pwr_cs, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.pwr_p075./T_pwr.pwr_cs, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.pwr_p005./T_pwr.pwr_cs, 'g');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.pwr_p095./T_pwr.pwr_cs, 'g');
    ylim([0, 1.2]);
    plot_title = strcat(uniquegen{i}, '-', uniquegensite{i});
    title(plot_title);
end

ax1 = subplot(3, 3, numel(cell_pwr)+1);
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.UP_RTD, '-b');
hold on;
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), -1.*T_rtd.DOWN_RTD, '-b');
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.FORECAST_ERROR_Brtd_Artd_solar, '-r')

ax2 = subplot(3, 3, numel(cell_pwr)+2);
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2), 'b');
title('RTD Binding solar');

linkaxes([ax1, ax2, ax_ibm],'x');
% suptitle('Normalized power');

%% Figure, power, RTD 4
figure();
ax_ibm = nan(1, numel(cell_pwr));
for i = 1:numel(cell_pwr)
    ax_ibm(i) = subplot(3, 3, i);
    T_pwr = cell_pwr{i};
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.pwr_p050, 'k');
    hold on;
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.pwr_p025, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.pwr_p075, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.pwr_p005, 'g');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.pwr_p095, 'g');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.pwr_cs, 'r');
    plot_title = strcat(uniquegen{i}, '-', uniquegensite{i});
    title(plot_title);
end

ax1 = subplot(3, 3, numel(cell_pwr)+1);
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.UP_RTD, '-b');
hold on;
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), -1.*T_rtd.DOWN_RTD, '-b');
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.FORECAST_ERROR_Brtd_Artd_solar, '-r')

ax2 = subplot(3, 3, numel(cell_pwr)+2);
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), sum(T_rtd{:, {'Solar_NP15_B_RTD', 'Solar_SP15_B_RTD', 'Solar_ZP26_B_RTD'}}, 2), 'b');
title('RTD Binding solar');

linkaxes([ax1, ax2, ax_ibm],'x');
% suptitle('Normalized power');

%% Figures, clear-sky index, RTPD 1
figure();
ax_ibm = nan(1, numel(cell_pwr));
for i = 1:numel(cell_pwr)
    ax_ibm(i) = subplot(3, 3, i);
    T_pwr = cell_pwr{i};
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p050./T_pwr.ghi_cs, 'k');
    hold on;
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p025./T_pwr.ghi_cs, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p075./T_pwr.ghi_cs, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p005./T_pwr.ghi_cs, 'g');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p095./T_pwr.ghi_cs, 'g');
    ylim([0, 1.2]);
    plot_title = strcat(uniquegen{i}, '-', uniquegensite{i});
    title(plot_title);
end

ax1 = subplot(3, 3, numel(cell_pwr)+1);
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), T_rtpd.UP_RTPD, '-b');
hold on;
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), -1.*T_rtpd.DOWN_RTPD, '-b');
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), T_rtpd.error_max_solar, '-m')
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), T_rtpd.error_min_solar, '-r')

ax2 = subplot(3, 3, numel(cell_pwr)+2);
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2), 'b');
title('RTPD Binding solar');

linkaxes([ax1, ax2, ax_ibm],'x');

%% Figures, GHI, RTPD 2
figure();
ax_ibm = nan(1, numel(cell_pwr));
for i = 1:numel(cell_pwr)
    ax_ibm(i) = subplot(3, 3, i);
    T_pwr = cell_pwr{i};
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p050, 'k');
    hold on;
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p025, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p075, 'b');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p005, 'g');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_p095, 'g');
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.ghi_cs, 'r');
    plot_title = strcat(uniquegen{i}, '-', uniquegensite{i});
    title(plot_title);
end

ax1 = subplot(3, 3, numel(cell_pwr)+1);
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), T_rtpd.UP_RTPD, '-b');
hold on;
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), -1.*T_rtpd.DOWN_RTPD, '-b');
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), T_rtpd.error_max_solar, '-r')
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), T_rtpd.error_min_solar, '-r')

ax2 = subplot(3, 3, numel(cell_pwr)+2);
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2), 'b');
title('RTPD Binding solar');

linkaxes([ax1, ax2, ax_ibm],'x');

%% Explore the relationship between RTD binding forecast and some indicators

% Figure comparing kt p50 fluctuation and RTD binding forecast
figure();
for i = 1:numel(cell_pwr)
    ax_ibm(i) = subplot(numel(cell_pwr)+1, 1, i);
    T_pwr = cell_pwr{i};
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), [nan;diff(T_pwr.kt_p050)], 'k');
    ylim([-0.2, 0.2]);
    plot_title = strcat(uniquegen{i}, '-', uniquegensite{i});
    title(plot_title);
end
ax2 = subplot(numel(cell_pwr)+1, 1, numel(cell_pwr)+1);
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2), 'b');
linkaxes([ax2, ax_ibm],'x');

% Figure comparing clear-sky index width and RTD binding forecast
figure();
for i = 1:numel(cell_pwr)
    ax_ibm(i) = subplot(numel(cell_pwr)+1, 1, i);
    T_pwr = cell_pwr{i};
    stairs(T_pwr.TIME-duration(0, dt_rtpd, 0), T_pwr.kt_width, 'k');
    ylim([-0.2, 0.2]);
    plot_title = strcat(uniquegen{i}, '-', uniquegensite{i});
    title(plot_title);
end
ax2 = subplot(numel(cell_pwr)+1, 1, numel(cell_pwr)+1);
stairs(T_rtpd.TIME-duration(0, dt_rtpd, 0), sum(T_rtpd{:, {'Solar_NP15_RTPD', 'Solar_SP15_RTPD', 'Solar_ZP26_RTPD'}}, 2), 'b');
linkaxes([ax2, ax_ibm],'x');

%% kNN: Explore RMSE of kp vs. FRP, RTD, i.e., only use delta kt as classifier
% Note that kp here is the ratio of total power over clear-sky power in CAISO

% Select month
this_month = 8;

% Baseline
T_baseline_rtd = array2table(unique(T_rtd.HOUR_START(T_rtd.HOUR_START.Month==this_month)), 'VariableNames', {'HOUR_START'}); % Result container
% T_rtd_DATE = datetime(T_rtd.HOUR_START.Year, T_rtd.HOUR_START.Month, T_rtd.HOUR_START.Day, 'TimeZone', 'UTC');
% T_rtd_DATE_unique = unique(T_rtd_DATE);
for i = 1: size(T_baseline_rtd, 1)
    this_date = datetime(T_baseline_rtd.HOUR_START.Year(i), T_baseline_rtd.HOUR_START.Month(i), T_baseline_rtd.HOUR_START.Day(i), 'TimeZone', 'UTC');
    this_hour = T_baseline_rtd.HOUR_START.Hour(i);
    selected_days = (T_rtd.TIME_START<this_date)&(T_rtd.TIME_START>=this_date-days(30))&(T_rtd.HOUR_START.Hour==this_hour); % We use 30 previous days
%     selected_days = (T_rtd_DATE<which_date)&(T_rtd_DATE>=which_date-days(30))&(T_rtd.HOUR_START.Hour==which_hour); % We use 30 previous days
%     valid_days = T_rtd_DATE_unique((T_rtd_DATE_unique<which_date)&(T_rtd_DATE_unique>=which_date-days(30)));
%     selected_days = (ismember(T_rtd_DATE, valid_days))&(T_rtd.HOUR_START.Hour==which_hour);
    sample_error = T_rtd{selected_days, 'FORECAST_ERROR_Brtd_Artd'};
%     T_sample = T_pwr_hourly((T_pwr_hourly.DATE<which_date)&(T_pwr_hourly.HOUR_START.Hour==which_hour), :);
%     T_sample_sorted = T_sample(T_sample.DATE>=which_date-days(30), :); % We use 30 previous days
%     selected_days = ismember(datetime(T_rtd.TIME_START.Year, T_rtd.TIME_START.Month, T_rtd.TIME_START.Day), T_sample_sorted.DATE(1:30)); % 30 the nearest days
%     sample_error = T_rtd{selected_days&(T_rtd.TIME_START.Hour==which_hour), 'FORECAST_ERROR_Brtd_Artd'};
    [f,x] = ecdf(sample_error(:));
    T_baseline_rtd.FRU(i) = interp1(f, x, 0.975);
    T_baseline_rtd.FRD(i) = interp1(f, x, 0.025);
end

T_rtd_hourly = grpstats(T_rtd(:, {'HOUR_START', 'kp'}), {'HOUR_START'}, 'std');
T_rtd_hourly.DATE = datetime(T_rtd_hourly.HOUR_START.Year, T_rtd_hourly.HOUR_START.Month, T_rtd_hourly.HOUR_START.Day, 'TimeZone', 'UTC');

% Test one month knn
T_stdbased_rtd = T_rtd_hourly(T_rtd_hourly.HOUR_START.Month==this_month, :); % Result container
for i = 1: size(T_stdbased_rtd, 1)
    this_date = datetime(T_stdbased_rtd.HOUR_START.Year(i), T_stdbased_rtd.HOUR_START.Month(i), T_stdbased_rtd.HOUR_START.Day(i), 'TimeZone', 'UTC');
    this_hour = T_stdbased_rtd.HOUR_START.Hour(i);
    T_sample = T_rtd_hourly((T_rtd_hourly.DATE<this_date)&(T_rtd_hourly.HOUR_START.Hour==this_hour), :);
    if any(isnan(T_stdbased_rtd{i, 'std_kp'}))
        T_sample_sorted = T_sample(T_sample.DATE>=this_date-days(30), :); % We use 30 previous days
    else
        T_sample.dist = abs(T_sample{:, 'std_kp'}-T_stdbased_rtd{i, 'std_kp'});
        T_sample_sorted = sortrows(T_sample, 'dist');
    end
    selected_days = ismember(datetime(T_rtd.TIME_START.Year, T_rtd.TIME_START.Month, T_rtd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:30)); % 30 the nearest days
    sample_error = T_rtd{selected_days&(T_rtd.TIME_START.Hour==this_hour), 'FORECAST_ERROR_Brtd_Artd'};
    [f,x] = ecdf(sample_error(:));
    T_stdbased_rtd.FRU(i) = interp1(f, x, 0.975);
    T_stdbased_rtd.FRD(i) = interp1(f, x, 0.025);
end

% Visualize the results, RTD
figure();
T_stdbased_rtd.TIME = T_stdbased_rtd.HOUR_START + duration(1, 0, 0); % We always use the end of an interval as time stamp
T_baseline_rtd.TIME = T_baseline_rtd.HOUR_START + duration(1, 0, 0); % We always use the end of an interval as time stamp

stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.UP_RTD, '-b');
hold on;
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), -1.*T_rtd.DOWN_RTD, '-b');
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.FORECAST_ERROR_Brtd_Artd, '-r')
stairs(T_stdbased_rtd.TIME-duration(1, 0, 0), T_stdbased_rtd.FRU, '-g');
stairs(T_stdbased_rtd.TIME-duration(1, 0, 0), T_stdbased_rtd.FRD, '-g');
stairs(T_baseline_rtd.TIME-duration(1, 0, 0), T_baseline_rtd.FRU, '-k');
stairs(T_baseline_rtd.TIME-duration(1, 0, 0), T_baseline_rtd.FRD, '-k');

% Use box plot to show changes of FRP, RTD
fru_compare_rtd = [T_baseline_rtd{(T_baseline_rtd.HOUR_START.Hour>=16)&(T_baseline_rtd.HOUR_START.Hour<=24), {'FRU'}} T_stdbased_rtd{(T_stdbased_rtd.HOUR_START.Hour>=16)&(T_stdbased_rtd.HOUR_START.Hour<=24), {'FRU'}}];
delta_fru_percent = (fru_compare_rtd(:, 2)-fru_compare_rtd(:, 1))./fru_compare_rtd(:, 1);

figure();
boxplot(reshape(delta_fru_percent, 8, size(delta_fru_percent, 1)/8)', 'Label', {'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
xlabel('Time (PST)');
ylabel('Relative change to baseline');
title('FRU');

frd_compare_rtd = [T_baseline_rtd{(T_baseline_rtd.HOUR_START.Hour>=16)&(T_baseline_rtd.HOUR_START.Hour<=24), {'FRD'}} T_stdbased_rtd{(T_stdbased_rtd.HOUR_START.Hour>=16)&(T_stdbased_rtd.HOUR_START.Hour<=24), {'FRD'}}];
delta_frd_percent = -(frd_compare_rtd(:, 2)-frd_compare_rtd(:, 1))./frd_compare_rtd(:, 1);

figure();
% boxplot(reshape(delta_frd_percent.*100, 8, size(delta_frd_percent, 1)/8)', 'Label', {'16', '17', '18', '19', '20', '21', '22', '23'});
boxplot(reshape(delta_frd_percent, 8, size(delta_frd_percent, 1)/8)', 'Label', {'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
xlabel('Time (PST)');
ylabel('Relative change to baseline');
title('FRD');

frp_need_rtd = T_rtd{ismember(T_rtd.HOUR_START, T_baseline_rtd.HOUR_START)&(T_rtd.HOUR_START.Hour>=16)&(T_rtd.HOUR_START.Hour<=24), 'FORECAST_ERROR_Brtd_Artd'};
frp_need_rtd = reshape(frp_need_rtd, 12, size(frp_need_rtd, 1)/12)';

fru_shortage_baseline_rtd = frp_need_rtd - fru_compare_rtd(:, 1);
fru_shortage_results_rtd  = frp_need_rtd - fru_compare_rtd(:, 2);
fru_n_shortage_hourly_results_rtd = sum(reshape(sum(fru_shortage_results_rtd>0, 2), 8, size(fru_shortage_results_rtd, 1)/8), 2);
fru_n_not_nan_hourly_results_rtd  = sum(reshape(~isnan(fru_shortage_results_rtd), 8, numel(fru_shortage_results_rtd)/8), 2);
fru_n_shortage_hourly_baseline_rtd = sum(reshape(sum(fru_shortage_baseline_rtd>0, 2), 8, size(fru_shortage_baseline_rtd, 1)/8), 2);
fru_n_not_nan_hourly_baseline_rtd  = sum(reshape(~isnan(fru_shortage_baseline_rtd), 8, numel(fru_shortage_baseline_rtd)/8), 2);
figure(); 
bar([fru_n_shortage_hourly_results_rtd./fru_n_not_nan_hourly_results_rtd fru_n_shortage_hourly_baseline_rtd./fru_n_not_nan_hourly_baseline_rtd]); % Percentage of violations each hour
legend({'Baseline', 'New'});
set(gca,'xticklabel',{'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
ylabel('Time of FRP shortage');
title('FRU');

frd_shortage_baseline_rtd = frp_need_rtd - frd_compare_rtd(:, 1);
frd_shortage_results_rtd  = frp_need_rtd - frd_compare_rtd(:, 2);
frd_n_shortage_hourly_results_rtd = sum(reshape(sum(frd_shortage_results_rtd<0, 2), 8, size(frd_shortage_results_rtd, 1)/8), 2);
frd_n_not_nan_hourly_results_rtd  = sum(reshape(~isnan(frd_shortage_results_rtd), 8, numel(frd_shortage_results_rtd)/8), 2);
frd_n_shortage_hourly_baseline_rtd = sum(reshape(sum(frd_shortage_baseline_rtd<0, 2), 8, size(frd_shortage_baseline_rtd, 1)/8), 2);
frd_n_not_nan_hourly_baseline_rtd  = sum(reshape(~isnan(frd_shortage_baseline_rtd), 8, numel(frd_shortage_baseline_rtd)/8), 2);
figure(); 
bar([frd_n_shortage_hourly_results_rtd./frd_n_not_nan_hourly_results_rtd frd_n_shortage_hourly_baseline_rtd./frd_n_not_nan_hourly_baseline_rtd]); % Percentage of violations each hour
legend({'Baseline', 'New'});
set(gca,'xticklabel',{'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
ylabel('Time of FRP shortage');
title('FRD');


%% kNN (kt and width of forecast, one site) starts here, RTD
T_pwr = cell_pwr{5}; % The 5th is CA_Topaz site
T_pwr.TIME_START = T_pwr.TIME - duration(0, dt_rtpd, 0);
T_pwr.HOUR_START = datetime(T_pwr.TIME_START.Year, T_pwr.TIME_START.Month, T_pwr.TIME_START.Day, T_pwr.TIME_START.Hour, 0, 0, 'TimeZone', 'UTC');
T_pwr.dkt = [nan; diff(T_pwr.kt_p050)]; % delta kT
T_pwr.dkt_sq = [nan; diff(T_pwr.kt_p050)].^2; % % (delta kT)^2
T_pwr_hourly = grpstats(T_pwr(:, {'HOUR_START', 'dkt_sq', 'kt_width'}), {'HOUR_START'}, 'mean'); 
T_pwr_hourly.u = sqrt(T_pwr_hourly.mean_dkt_sq); % This is variability within each hour, following Inman et al. 2013, section 2.6.1
T_pwr_hourly.DATE = datetime(T_pwr_hourly.HOUR_START.Year, T_pwr_hourly.HOUR_START.Month, T_pwr_hourly.HOUR_START.Day, 'TimeZone', 'UTC');

% Select month
this_month = 10;

% Baseline
T_baseline_rtd = array2table(unique(T_rtd.HOUR_START(T_rtd.HOUR_START.Month==this_month)), 'VariableNames', {'HOUR_START'}); % Result container
% T_rtd_DATE = datetime(T_rtd.HOUR_START.Year, T_rtd.HOUR_START.Month, T_rtd.HOUR_START.Day, 'TimeZone', 'UTC');
% T_rtd_DATE_unique = unique(T_rtd_DATE);
for i = 1: size(T_baseline_rtd, 1)
    this_date = datetime(T_baseline_rtd.HOUR_START.Year(i), T_baseline_rtd.HOUR_START.Month(i), T_baseline_rtd.HOUR_START.Day(i), 'TimeZone', 'UTC');
    this_hour = T_baseline_rtd.HOUR_START.Hour(i);
    selected_days = (T_rtd.TIME_START<this_date)&(T_rtd.TIME_START>=this_date-days(30))&(T_rtd.HOUR_START.Hour==this_hour); % We use 30 previous days
%     selected_days = (T_rtd_DATE<which_date)&(T_rtd_DATE>=which_date-days(30))&(T_rtd.HOUR_START.Hour==which_hour); % We use 30 previous days
%     valid_days = T_rtd_DATE_unique((T_rtd_DATE_unique<which_date)&(T_rtd_DATE_unique>=which_date-days(30)));
%     selected_days = (ismember(T_rtd_DATE, valid_days))&(T_rtd.HOUR_START.Hour==which_hour);
    sample_error = T_rtd{selected_days, 'FORECAST_ERROR_Brtd_Artd'};
%     T_sample = T_pwr_hourly((T_pwr_hourly.DATE<which_date)&(T_pwr_hourly.HOUR_START.Hour==which_hour), :);
%     T_sample_sorted = T_sample(T_sample.DATE>=which_date-days(30), :); % We use 30 previous days
%     selected_days = ismember(datetime(T_rtd.TIME_START.Year, T_rtd.TIME_START.Month, T_rtd.TIME_START.Day), T_sample_sorted.DATE(1:30)); % 30 the nearest days
%     sample_error = T_rtd{selected_days&(T_rtd.TIME_START.Hour==which_hour), 'FORECAST_ERROR_Brtd_Artd'};
    [f,x] = ecdf(sample_error(:));
    T_baseline_rtd.FRU(i) = interp1(f, x, 0.975);
    T_baseline_rtd.FRD(i) = interp1(f, x, 0.025);
end

% Test one month knn
T_results_rtd = T_pwr_hourly(T_pwr_hourly.DATE.Month==this_month, :); % Result container
% T_results_rtd = array2table(unique(T_rtd.HOUR_START(T_rtd.HOUR_START.Month==5)), 'VariableNames', {'HOUR_START'}); % Result container
for i = 1: size(T_results_rtd, 1)
    this_date = datetime(T_results_rtd.HOUR_START.Year(i), T_results_rtd.HOUR_START.Month(i), T_results_rtd.HOUR_START.Day(i), 'TimeZone', 'UTC');
    this_hour = T_results_rtd.HOUR_START.Hour(i);
    T_sample = T_pwr_hourly((T_pwr_hourly.DATE<this_date)&(T_pwr_hourly.HOUR_START.Hour==this_hour), :);
    if any(isnan(T_results_rtd{i, {'mean_kt_width', 'u'}}))
        T_sample_sorted = T_sample(T_sample.DATE>=this_date-days(30), :); % We use 30 previous days
    else
        T_sample.dist = sqrt(sum((T_sample{:, {'mean_kt_width', 'u'}}-T_results_rtd{i, {'mean_kt_width', 'u'}}).^2, 2));
        T_sample_sorted = sortrows(T_sample, 'dist');
    end
    selected_days = ismember(datetime(T_rtd.TIME_START.Year, T_rtd.TIME_START.Month, T_rtd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:30)); % 30 the nearest days
    sample_error = T_rtd{selected_days&(T_rtd.TIME_START.Hour==this_hour), 'FORECAST_ERROR_Brtd_Artd'};
    [f,x] = ecdf(sample_error(:));
    T_results_rtd.FRU(i) = interp1(f, x, 0.975);
    T_results_rtd.FRD(i) = interp1(f, x, 0.025);
end

% Visualize the results, RTD
figure();
T_results_rtd.TIME   = T_results_rtd.HOUR_START + duration(1, 0, 0); % We always use the end of an interval as time stamp
T_baseline_rtd.TIME = T_baseline_rtd.HOUR_START + duration(1, 0, 0); % We always use the end of an interval as time stamp

stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.UP_RTD, '-b');
hold on;
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), -1.*T_rtd.DOWN_RTD, '-b');
stairs(T_rtd.TIME-duration(0, dt_rtd, 0), T_rtd.FORECAST_ERROR_Brtd_Artd, '-r')
stairs(T_results_rtd.TIME-duration(1, 0, 0), T_results_rtd.FRU, '-g');
stairs(T_results_rtd.TIME-duration(1, 0, 0), T_results_rtd.FRD, '-g');
stairs(T_baseline_rtd.TIME-duration(1, 0, 0), T_baseline_rtd.FRU, '-k');
stairs(T_baseline_rtd.TIME-duration(1, 0, 0), T_baseline_rtd.FRD, '-k');

% Use box plot to show changes of FRP, RTD
fru_compare_rtd = [T_baseline_rtd{(T_baseline_rtd.HOUR_START.Hour>=16)&(T_baseline_rtd.HOUR_START.Hour<=24), {'FRU'}} T_results_rtd{(T_results_rtd.HOUR_START.Hour>=16)&(T_results_rtd.HOUR_START.Hour<=24), {'FRU'}}];
delta_fru_percent = (fru_compare_rtd(:, 2)-fru_compare_rtd(:, 1))./fru_compare_rtd(:, 1);

figure();
boxplot(reshape(delta_fru_percent, 8, size(delta_fru_percent, 1)/8)', 'Label', {'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
xlabel('Time (PST)');
ylabel('Relative change to baseline');
title('FRU');

frd_compare_rtd = [T_baseline_rtd{(T_baseline_rtd.HOUR_START.Hour>=16)&(T_baseline_rtd.HOUR_START.Hour<=24), {'FRD'}} T_results_rtd{(T_results_rtd.HOUR_START.Hour>=16)&(T_results_rtd.HOUR_START.Hour<=24), {'FRD'}}];
delta_frd_percent = -(frd_compare_rtd(:, 2)-frd_compare_rtd(:, 1))./frd_compare_rtd(:, 1);

figure();
% boxplot(reshape(delta_frd_percent.*100, 8, size(delta_frd_percent, 1)/8)', 'Label', {'16', '17', '18', '19', '20', '21', '22', '23'});
boxplot(reshape(delta_frd_percent, 8, size(delta_frd_percent, 1)/8)', 'Label', {'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
xlabel('Time (PST)');
ylabel('Relative change to baseline');
title('FRD');

frp_need_rtd = T_rtd{ismember(T_rtd.HOUR_START, T_baseline_rtd.HOUR_START)&(T_rtd.HOUR_START.Hour>=16)&(T_rtd.HOUR_START.Hour<=24), 'FORECAST_ERROR_Brtd_Artd'};
frp_need_rtd = reshape(frp_need_rtd, 12, size(frp_need_rtd, 1)/12)';

fru_shortage_baseline_rtd = frp_need_rtd - fru_compare_rtd(:, 1);
fru_shortage_results_rtd  = frp_need_rtd - fru_compare_rtd(:, 2);
fru_n_shortage_hourly_results_rtd = sum(reshape(sum(fru_shortage_results_rtd>0, 2), 8, size(fru_shortage_results_rtd, 1)/8), 2);
fru_n_not_nan_hourly_results_rtd  = sum(reshape(~isnan(fru_shortage_results_rtd), 8, numel(fru_shortage_results_rtd)/8), 2);
fru_n_shortage_hourly_baseline_rtd = sum(reshape(sum(fru_shortage_baseline_rtd>0, 2), 8, size(fru_shortage_baseline_rtd, 1)/8), 2);
fru_n_not_nan_hourly_baseline_rtd  = sum(reshape(~isnan(fru_shortage_baseline_rtd), 8, numel(fru_shortage_baseline_rtd)/8), 2);
figure(); 
bar([fru_n_shortage_hourly_results_rtd./fru_n_not_nan_hourly_results_rtd fru_n_shortage_hourly_baseline_rtd./fru_n_not_nan_hourly_baseline_rtd]); % Percentage of violations each hour
legend({'Baseline', 'New'});
set(gca,'xticklabel',{'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
ylabel('Time of FRP shortage');
title('FRU');

frd_shortage_baseline_rtd = frp_need_rtd - frd_compare_rtd(:, 1);
frd_shortage_results_rtd  = frp_need_rtd - frd_compare_rtd(:, 2);
frd_n_shortage_hourly_results_rtd = sum(reshape(sum(frd_shortage_results_rtd<0, 2), 8, size(frd_shortage_results_rtd, 1)/8), 2);
frd_n_not_nan_hourly_results_rtd  = sum(reshape(~isnan(frd_shortage_results_rtd), 8, numel(frd_shortage_results_rtd)/8), 2);
frd_n_shortage_hourly_baseline_rtd = sum(reshape(sum(frd_shortage_baseline_rtd<0, 2), 8, size(frd_shortage_baseline_rtd, 1)/8), 2);
frd_n_not_nan_hourly_baseline_rtd  = sum(reshape(~isnan(frd_shortage_baseline_rtd), 8, numel(frd_shortage_baseline_rtd)/8), 2);
figure(); 
bar([frd_n_shortage_hourly_results_rtd./frd_n_not_nan_hourly_results_rtd frd_n_shortage_hourly_baseline_rtd./frd_n_not_nan_hourly_baseline_rtd]); % Percentage of violations each hour
legend({'Baseline', 'New'});
set(gca,'xticklabel',{'8-9', '9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16'});
ylabel('Time of FRP shortage');
title('FRD');

%% RTPD requirements
% fru_rtpd = nan(length(1498:1506), 1);
% frd_rtpd = nan(length(1498:1506), 1);
% for i = 1498:1506
%     which_h = T_pwr_hourly{i, 'Hour'};
%     which_date = datetime(T_pwr_hourly{i, 'Year'}, T_pwr_hourly{i, 'Month'}, T_pwr_hourly{i, 'Day'});
%     T_tmp = T_pwr_hourly((T_pwr_hourly.DATE<which_date)&(T_pwr_hourly.Hour==16), :);
%     T_tmp.dist = sqrt(sum((T_tmp{:, {'mean_kt_width', 'u'}}-T_pwr_hourly{i, {'mean_kt_width', 'u'}}).^2, 2));
%     T_tmp_sorted = sortrows(T_tmp, 'dist');
%     selected_days = ismember(datetime(T_rtpd.TIME.Year, T_rtpd.TIME.Month, T_rtpd.TIME.Day), T_tmp_sorted.DATE(1:30)); % 30 the nearest days
%     sample_error_max = T_rtpd{selected_days, 'error_max'};
%     sample_error_min = T_rtpd{selected_days, 'error_min'};
%     [f_max,x_max] = ecdf(sample_error_max(:));
%     fru_rtpd(i-1498+1, 1) = interp1(f_max, x_max, 0.975);
%     [f_min,x_min] = ecdf(sample_error_min(:));
%     frd_rtpd(i-1498+1, 1) = interp1(f_min, x_min, 0.025);
% end
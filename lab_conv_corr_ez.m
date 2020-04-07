function lab_conv_corr_ez()

% prepare_data_uscrn();
% prepare_data_nsrdb();
% prepare_ibm_measurement(4);

paper_twosites();

% conv_corr_ez_nd();
% 
% test_cap_and_floor();

end

function paper_twosites()
addpath('ndhist\histcn'); % Add the path to the N-Dimensional histcount function

lat_all = [35.37	39.13	35.53	40.73	37.41	34.29	40.25	38.49	34.45	34.13];
lon_all = [-120.18	-123.06	-118.62	-123.94	-119.74	-117.5	-123.3	-122.7	-119.7	-117.94];
cell_sites = {'CA_Topaz', 'COWC1', 'DEMC1', 'KNNC1', 'MIAC1', 'MNCC1', 'RLKC1', 'RSAC1', 'SBVC1', 'STFC1'};

% Probabilistic forecast
cell_csvparam{1}  = 'C:\Users\bxl180002\Downloads\RampSolar\Cong\Results_1HA\CA_Topaz.param.Laplace.csv';
cell_csvparam{2}  = 'C:\Users\bxl180002\Downloads\RampSolar\Cong\Results_1HA\COWC1.param.Laplace.csv';
cell_csvparam{3}  = 'C:\Users\bxl180002\Downloads\RampSolar\Cong\Results_1HA\DEMC1.param.Laplace.csv';
cell_csvparam{4}  = 'C:\Users\bxl180002\Downloads\RampSolar\Cong\Results_1HA\KNNC1.param.Laplace.csv';
cell_csvparam{5}  = 'C:\Users\bxl180002\Downloads\RampSolar\Cong\Results_1HA\MIAC1.param.Laplace.csv';
cell_csvparam{6}  = 'C:\Users\bxl180002\Downloads\RampSolar\Cong\Results_1HA\MNCC1.param.Laplace.csv';
cell_csvparam{7}  = 'C:\Users\bxl180002\Downloads\RampSolar\Cong\Results_1HA\RLKC1.param.Laplace.csv';
cell_csvparam{8}  = 'C:\Users\bxl180002\Downloads\RampSolar\Cong\Results_1HA\RSAC1.param.Laplace.csv';
cell_csvparam{9}  = 'C:\Users\bxl180002\Downloads\RampSolar\Cong\Results_1HA\SBVC1.param.Laplace.csv';
cell_csvparam{10} = 'C:\Users\bxl180002\Downloads\RampSolar\Cong\Results_1HA\STFC1.param.Laplace.csv';

% Original NSRDB data
cell_csvnsrdb{1}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\CA_Topaz\96633_35.37_-120.18_2018.csv';
cell_csvnsrdb{2}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\COWC1\138221_39.13_-123.06_2018.csv';
cell_csvnsrdb{3}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\DEMC1\98317_35.53_-118.62_2018.csv';
cell_csvnsrdb{4}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\KNNC1\157537_40.73_-123.94_2018.csv';
cell_csvnsrdb{5}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\MIAC1\118489_37.41_-119.74_2018.csv';
cell_csvnsrdb{6}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\MNCC1\85923_34.29_-117.5_2018.csv';
cell_csvnsrdb{7}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\RLKC1\151675_40.25_-123.3_2018.csv';
cell_csvnsrdb{8}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\RSAC1\130670_38.49_-122.7_2018.csv';
cell_csvnsrdb{9}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\SBVC1\87449_34.45_-119.7_2018.csv';
cell_csvnsrdb{10} = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\STFC1\84346_34.13_-117.94_2018.csv';


% For testing purpose, we start by two sites only
% i_selected = [1, 3, 5, 6, 10];
i_selected = [1, 5];

lat_all = lat_all(i_selected);
lon_all = lon_all(i_selected);
cell_sites = cell_sites(i_selected);
cell_csvparam   = cell_csvparam(i_selected);
cell_csvnsrdb   = cell_csvnsrdb(i_selected);
N = numel(cell_csvparam);

% Load data, parameters of probabilistic forecasts
cell_Tparam = cell(N, 1);
for i = 1:N
    T = readtable(cell_csvparam{i});
    T.TIME = datetime(T.Year, T.Month, T.Day, T.Hour, T.Minute, 0, 'TimeZone', 'UTC');
    lat = lat_all(i);
    lon = lon_all(i);
    Location = pvl_makelocationstruct(lat, lon);
    Location.altitude = 0;
    Time.UTCOffset = zeros(size(T.TIME, 1), 1); % tarray must be UTC time
    Time.year   = T.TIME.Year;
    Time.month  = T.TIME.Month;
    Time.day    = T.TIME.Day;
    Time.hour   = T.TIME.Hour;
    Time.minute = T.TIME.Minute;
    Time.second = T.TIME.Second;
    [SunAz, SunEl, ApparentSunEl, SolarTime] = pvl_ephemeris(Time,Location);
    T.SunAz = SunAz;
    T.SunEl = SunEl;
    T.ApparentSunEl = ApparentSunEl;
    T.GHI_CS = pvl_clearsky_haurwitz(90-ApparentSunEl);
    T.k = T.GHI./T.GHI_CS;
    cell_Tparam{i} = T;
end

% Load data, original NSRDB data
cell_Tnsrdb = cell(N, 1);
for i = 1:numel(cell_csvnsrdb)
    T = readtable(cell_csvnsrdb{i}, 'HeaderLines', 2);
    T.TIME = datetime(T.Year, T.Month, T.Day, T.Hour, T.Minute, 0, 'TimeZone', 'UTC');
    lat = lat_all(i);
    lon = lon_all(i);
    Location = pvl_makelocationstruct(lat, lon);
    Location.altitude = 0;
    Time.UTCOffset = zeros(size(T.TIME, 1), 1); % tarray must be UTC time
    Time.year   = T.TIME.Year;
    Time.month  = T.TIME.Month;
    Time.day    = T.TIME.Day;
    Time.hour   = T.TIME.Hour;
    Time.minute = T.TIME.Minute;
    Time.second = T.TIME.Second;
    [SunAz, SunEl, ApparentSunEl, SolarTime] = pvl_ephemeris(Time,Location);
    T.SunAz = SunAz;
    T.SunEl = SunEl;
    T.ApparentSunEl = ApparentSunEl;
    T.GHI_CS = pvl_clearsky_haurwitz(90-ApparentSunEl);
    T.k = T.GHI./T.GHI_CS;
    cell_Tnsrdb{i} = T;
end

% Only use 0 < k < 5
tmp = nan(size(cell_Tnsrdb{i}, 1), N);
for i = 1: N
    tmp(:, i) = cell_Tnsrdb{i}.k;
end
i_selected = ~any((tmp<=0)|(tmp>=5)|isnan(tmp), 2);

% Marginal distributions of k
cell_gmmodel = cell(N, 1);
for i = 1: N
    T = cell_Tnsrdb{i};
    k_selected = T.k(T.SunEl>5);
    pd_gmmodel = fitgmdist(k_selected, 2); % Use two components
    cell_gmmodel{i} = pd_gmmodel;
end

% Fit Gaussian copuola function and generate random copula samples
ts_cdf = nan(sum(i_selected), N);
for i = 1:N
    ts_cdf(:, i) = cdf(cell_gmmodel{i}, cell_Tnsrdb{i}.k(i_selected));
end
ts_cdf(ts_cdf==1) = 1-eps;

% copula_type = 'Gaussian';
% copula_type = 'Frank';
% copula_type = 'Gumbel';
copula_type = 'Clayton';
paramhat = copulafit(copula_type, ts_cdf);
u_rnd = copularnd(copula_type, paramhat, 1E7); % Monte Carlo simulation


% Convolution start here
T = size(cell_Tparam{1}, 1);
crps_conv = nan(T, 1);
crps_conv_corr = nan(T, 1);
for t = 1:T
    disp(t);
    if cell_Tparam{i}.SunEl(t)<5 % We only look at sun elevation > 5 degree
        continue
    end
    cell_bincenter = cell(N, 1);
    cell_p = cell(N, 1);
    cell_edge_cdf = cell(N, 1); % Cumulative probability at the edges of grids
    power_observe = nan(N, 1);
    for i = 1:N
        T = cell_Tparam{i};
        power_observe(i) = T.GHI(t);
        G = 10:10:1300;
        % Laplace CDF
        MU = T.rf2(t);
        B  = T.Param(t);
        CDF = nan(numel(G), 1);
        CDF(G<=MU) = 1/2.*exp((G(G<=MU)-MU)./B);
        CDF(G>MU) = 1 - 1/2.*exp(-(G(G>MU)-MU)./B);

        % Discretization
        i_start = max(find(CDF<=0.001));
        i_end = min(find(CDF>=0.999));
        if isempty(i_start)
            i_start = 1;
        end
        p = diff(CDF(i_start:i_end));
        bincenter = (G(i_start:i_end-1) + G(i_start+1:i_end))/2;

        cell_bincenter{i} = bincenter; % Bin center locations
        cell_p{i} = p; % Probability in each bin
        cell_edge_cdf{i} = CDF(i_start:i_end);
    end
    power_observe(power_observe>=1000) = 1000; % Capped at 1000
    sumpower = sum(power_observe);

    % Calculate mean copula density over each hypercube
    switch N
        case 2
            [X1, X2] = ndgrid(cell_p{1}, cell_p{2}); % Grid size along each dimension
            grid_vol = X1.*X2; % Grid volume
            clear X1 X2; % Save some memory
            count_gaussian = histcn(u_rnd, cell_edge_cdf{1}, cell_edge_cdf{2});
            copula_pdf_nd_nocorr = ones(numel(cell_bincenter{1}), numel(cell_bincenter{2}));
        case 5
            [X1, X2, X3, X4, X5] = ndgrid(cell_p{1}, cell_p{2}, cell_p{3}, cell_p{4}, cell_p{5}); % Grid size along each dimension
            grid_vol = X1.*X2.*X3.*X4.*X5; % Grid volume
            clear X1 X2 X3 X4 X5; % Save some memory
            count_gaussian = histcn(u_rnd, cell_edge_cdf{1}, cell_edge_cdf{2}, cell_edge_cdf{3}, cell_edge_cdf{4}, cell_edge_cdf{5});
            copula_pdf_nd_nocorr = ones(numel(cell_bincenter{1}), numel(cell_bincenter{2}), numel(cell_bincenter{3}), numel(cell_bincenter{4}), numel(cell_bincenter{5}));
    end

    copula_pdf_nd_gaussian = count_gaussian./(sum(count_gaussian(:)).*grid_vol);
    [x_conv,      p_conv]  = conv_ez_corr(cell_bincenter, cell_p, 10, copula_pdf_nd_nocorr);
    [x_conv_corr, p_conv_corr] = conv_ez_corr(cell_bincenter, cell_p, 10, copula_pdf_nd_gaussian);
    crps_conv(t) = crps(x_conv, cumsum(p_conv), sumpower, 10);
    crps_conv_corr(t) = crps(x_conv_corr, cumsum(p_conv_corr), sumpower, 10);
end

end

function prepare_data_uscrn()
% Load USCRN data, compare with IBM's nearest site
add_pvlib;
dt_rtd = 5;

cell_txt = cell(7, 1);
cell_txt{1} = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Bodega_6_WSW.txt';
cell_txt{2} = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Fallbrook_5_NE.txt';
cell_txt{3} = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Merced_23_WSW.txt';
cell_txt{4} = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Redding_12_WNW.txt';
cell_txt{5} = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Santa_Barbara_11_W.txt';
cell_txt{6} = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Stovepipe_Wells_1_SW.txt';
cell_txt{7} = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Yosemite_Village_12_W.txt';

cell_uscrn = cell(7, 1);

for i = 1: numel(cell_txt)
    txtfile = cell_txt{i};
    T_uscrn = readtable(txtfile, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
    T_uscrn.Properties.VariableNames = {'WBANNO', 'UTC_DATE', 'UTC_TIME', 'LST_DATE', 'LST_TIME', 'CRX_VN', 'LONGITUDE', 'LATITUDE', 'AIR_TEMPERATURE', 'PRECIPITATION', 'SOLAR_RADIATION', 'SR_FLAG', 'SURFACE_TEMPERATURE', 'ST_TYPE', 'ST_FLAG', 'RELATIVE_HUMIDITY', 'RH_FLAG', 'SOIL_MOISTURE_5', 'SOIL_TEMPERATURE_5', 'WETNESS', 'WET_FLAG', 'WIND_1_5', 'WIND_FLAG'};
    T_uscrn.TIME = datetime(T_uscrn.UTC_DATE, 'ConvertFrom', 'yyyymmdd') + timeofday(datetime(num2str(T_uscrn.UTC_TIME, '%04d'), 'Format', 'HHmm'));
    T_uscrn.TIME.TimeZone = 'UTC';
    lat = unique(T_uscrn.LATITUDE);
    lon = unique(T_uscrn.LONGITUDE);
    T_uscrn.GHI_CS = mean_clearsky_ghi(lat, lon, T_uscrn.TIME, dt_rtd);
    T_uscrn.k = T_uscrn.SOLAR_RADIATION./T_uscrn.GHI_CS;
    cell_uscrn{i} = T_uscrn;
end


% s_uscrn = 'Santa Barbara';
% s_uscrn = 'Yosemite';
% s_uscrn = 'Bodega';
% switch s_uscrn
%     case 'Santa Barbara'
%         ibm_site = 'SBVC1';
%         txtfile = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Santa_Barbara_11_W.txt';
%     case 'Yosemite'
%         ibm_site = 'MIAC1';
%         txtfile = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Yosemite_Village_12_W.txt';
%     case 'Bodega'
%         ibm_site = 'RSAC1';
%         txtfile = 'C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\USCRN_2019\CRNS0101-05-2019-CA_Bodega_6_WSW.txt';
% end
% 
% T_uscrn = readtable(txtfile, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
% T_uscrn.Properties.VariableNames = {'WBANNO', 'UTC_DATE', 'UTC_TIME', 'LST_DATE', 'LST_TIME', 'CRX_VN', 'LONGITUDE', 'LATITUDE', 'AIR_TEMPERATURE', 'PRECIPITATION', 'SOLAR_RADIATION', 'SR_FLAG', 'SURFACE_TEMPERATURE', 'ST_TYPE', 'ST_FLAG', 'RELATIVE_HUMIDITY', 'RH_FLAG', 'SOIL_MOISTURE_5', 'SOIL_TEMPERATURE_5', 'WETNESS', 'WET_FLAG', 'WIND_1_5', 'WIND_FLAG'};
% T_uscrn.TIME = datetime(T_uscrn.UTC_DATE, 'ConvertFrom', 'yyyymmdd') + timeofday(datetime(num2str(T_uscrn.UTC_TIME, '%04d'), 'Format', 'HHmm'));
% T_uscrn.TIME.TimeZone = 'UTC';

tmp = [];
for i = 1: numel(cell_txt)
    T_uscrn = cell_uscrn{i};
    tmp = [tmp T_uscrn.k(T_uscrn.UTC_DATE>20190101)];
end

stairs(T_uscrn.TIME-duration(0, dt_rtd, 0), T_uscrn.SOLAR_RADIATION, 'b');
title(strcat('RAWS-', ibm_site, ' vs USCRN-', s_uscrn));
ylim([0, 1300]);

end

function prepare_data_nsrdb()
% Preliminary study on the NSRDB data, IBM sites
cell_nsrdb = cell(10, 1);

lat_all = [35.37	39.13	35.53	40.73	37.41	34.29	40.25	38.49	34.45	34.13];
lon_all = [-120.18	-123.06	-118.62	-123.94	-119.74	-117.5	-123.3	-122.7	-119.7	-117.94];

cell_names = {'CA_Topaz', 'COWC1', 'DEMC1', 'KNNC1', 'MIAC1', 'MNCC1', 'RLKC1', 'RSAC1', 'SBVC1', 'STFC1'};

cell_nsrdb{1}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\CA_Topaz\96633_35.37_-120.18_2018.csv';
cell_nsrdb{2}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\COWC1\138221_39.13_-123.06_2018.csv';
cell_nsrdb{3}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\DEMC1\98317_35.53_-118.62_2018.csv';
cell_nsrdb{4}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\KNNC1\157537_40.73_-123.94_2018.csv';
cell_nsrdb{5}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\MIAC1\118489_37.41_-119.74_2018.csv';
cell_nsrdb{6}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\MNCC1\85923_34.29_-117.5_2018.csv';
cell_nsrdb{7}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\RLKC1\151675_40.25_-123.3_2018.csv';
cell_nsrdb{8}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\RSAC1\130670_38.49_-122.7_2018.csv';
cell_nsrdb{9}  = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\SBVC1\87449_34.45_-119.7_2018.csv';
cell_nsrdb{10} = 'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\STFC1\84346_34.13_-117.94_2018.csv';

cell_Tnsrdb = cell(numel(cell_nsrdb), 1);

for i = 1:numel(cell_nsrdb)
    T = readtable(cell_nsrdb{i}, 'HeaderLines', 2);
    T.TIME = datetime(T.Year, T.Month, T.Day, T.Hour, T.Minute, 0, 'TimeZone', 'UTC');
    lat = lat_all(i);
    lon = lon_all(i);
    Location = pvl_makelocationstruct(lat, lon);
    Location.altitude = 0;
    Time.UTCOffset = zeros(size(T.TIME, 1), 1); % tarray must be UTC time
    Time.year   = T.TIME.Year;
    Time.month  = T.TIME.Month;
    Time.day    = T.TIME.Day;
    Time.hour   = T.TIME.Hour;
    Time.minute = T.TIME.Minute;
    Time.second = T.TIME.Second;
    [SunAz, SunEl, ApparentSunEl, SolarTime] = pvl_ephemeris(Time,Location);
    T.SunAz = SunAz;
    T.SunEl = SunEl;
    T.ApparentSunEl = ApparentSunEl;
    T.GHI_CS = pvl_clearsky_haurwitz(90-ApparentSunEl);
%     T.GHI_CS = mean_clearsky_ghi(lat, lon, T.TIME, dt);
    T.k = T.GHI./T.GHI_CS;
%     T.k = T.GHI./T.ClearskyGHI; % Clear-sky index using NSRDB's clear-sky GHI

    cell_Tnsrdb{i} = T;
end

tmp = [];

for i = 1: numel(cell_Tnsrdb)
    tmp = [tmp cell_Tnsrdb{i}.k];
end
i_selected = ~any((tmp<=0)|(tmp>=5)|isnan(tmp), 2); % Only use 0 < k < 5
corrcoef(tmp(i_selected, :))
% tmp(any(isnan(tmp), 2), :) = [];
% tmp(any(tmp==1, 2), :) = [];
% corrcoef(tmp)

% Try different marginal distributions
for i = 1: 10
    T = cell_Tnsrdb{i};
    k_selected = T.k(T.SunEl>5);
    pd_gmmodel   = fitgmdist(k_selected, 2);
    pd_gammodel  = fitdist(k_selected, 'Gamma');
    pd_lognormal = fitdist(k_selected, 'Lognormal');
    pd_gaussian  = fitdist(k_selected, 'Normal');
    [y_pdf, bincenter, x_edges] = return_pdf(k_selected, 0.02);
    
    figure();
    subplot(2, 1, 1); 
    plot(bincenter, pdf(pd_gmmodel, bincenter(:)), 'b');
    hold on;
    plot(bincenter, pdf(pd_gammodel, bincenter(:)), 'r');
    plot(bincenter, pdf(pd_lognormal, bincenter(:)), 'g');
    plot(bincenter, pdf(pd_gaussian, bincenter(:)), 'g--');
    legend('GMM', 'Gamma', 'Lognormal', 'Normal');
    bar(bincenter, y_pdf);
    
    subplot(2, 1, 2); % pp plot
    cdf_actual = cumsum(y_pdf.*0.02);
    plot([0, 1], [0, 1], 'k');
    hold on; 
    plot(cdf_actual, cdf(pd_gmmodel, bincenter(:)), 'b');
    plot(cdf_actual, cdf(pd_gammodel, bincenter(:)), 'r');
    plot(cdf_actual, cdf(pd_lognormal, bincenter(:)), 'g');
    plot(cdf_actual, cdf(pd_gaussian, bincenter(:)), 'g--');
    
    suptitle(cell_names{i});
end
end

function prepare_ibm_measurement(m)
% Preliminary study on IBM's measurement data, i.e., hourly average from RAWS stations
if m == 4
    dir_work = 'C:\Users\bxl180002\git\SF2\IBM\April\ghi_actual';
elseif m == 5
    dir_work = 'C:\Users\bxl180002\git\SF2\IBM\May\ghi_actual';
end
dir_home = pwd;

% IBMsitenames   = {'MNCC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz'};
% SiteLatitude  = [34.31,   34.12,    37.41,   35.53,   35.38];
% SiteLongitude = [-117.5,-117.94,  -119.74, -118.63, -120.18];

IBMsitenames  = {'CA_Topaz','RSAC1', 'RLKC1', 'SBVC1', 'KNNC1', 'MIAC1', 'MNCC1', 'STFC1', 'DEMC1',  'COWC1'};
SiteLatitude  = [35.38,       38.47,   40.25,   34.45,   40.71,   37.41,   34.31,   34.12,   35.53,   39.12];
SiteLongitude = [-120.18,   -122.71, -123.31, -119.70, -123.92, -119.74, -117.50, -117.94, -118.63, -123.07];

dir_home = pwd;
ghi_actual = [];
ghi_clearsky = [];

for k = 1:length(IBMsitenames)
    ibm_site = IBMsitenames{k};
    csvname_read  = strcat('IBM_processed_', ibm_site, '.hourly.csv');
    
    cd(dir_work);
    M = csvread(csvname_read, 1, 0);
    cd(dir_home);
    
    Location = pvl_makelocationstruct(SiteLatitude(k),SiteLongitude(k)); %Altitude is optional
    Time.UTCOffset(1:size(M,1),1) = zeros(size(M,1), 1); % Because we use UTC time, so utc offset is zero
    Time.year(1:size(M,1),1)   = M(:, 1);
    Time.month(1:size(M,1),1)  = M(:, 2);
    Time.day(1:size(M,1),1)    = M(:, 3);
    Time.hour(1:size(M,1),1)   = M(:, 4);
    Time.minute(1:size(M,1),1) = M(:, 5);
    Time.second(1:size(M,1),1) = zeros(size(M,1), 1);
    [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
    
    ghi_clearsky(:, k) = pvl_clearsky_haurwitz(90-AppSunEl); % Clear-sky GHI
    ghi_actual(:, k) = M(:, end);
end

tarray = datetime(M(:, 1), M(:, 2), M(:, 3), M(:, 4), M(:, 5), M(:, 6), 'TimeZone', 'UTC');
tarray.TimeZone = '-08:00'; % Let's just use PST, so that 12:00 is noon
tarray_pst = datetime(tarray.Year, tarray.Month, tarray.Day, tarray.Hour, tarray.Minute, tarray.Second);

kcs = ghi_actual./ghi_clearsky;
kcs_midday = kcs((tarray_pst.Hour>=7) & (tarray_pst.Hour<17), :);
kcs_other  = kcs((tarray_pst.Hour<7) | (tarray_pst.Hour>=17), :);

figure();
subplot(2, 1, 1);
hist(kcs_midday(:), 50);
title('7am to 5pm (PST)');
xlabel('Clear-sky index');
subplot(2, 1, 2);
hist(kcs_other(:), 50);
title('Other time');
xlabel('Clear-sky index');

dist_gc = nan(numel(IBMsitenames), numel(IBMsitenames));
for i = 1: size(kcs_midday, 2) - 1
    lat1 = SiteLatitude(i);
    lon1 = SiteLongitude(i);
    dist_gc(i, i) = 0;
    for j = i+1: size(kcs_midday, 2)
        lat2 = SiteLatitude(j);
        lon2 = SiteLongitude(j);
        d = DISTANCE(lat1/180*pi, lon1/180*pi, lat2/180*pi, lon2/180*pi);
        dist_gc(i, j) = d;
        dist_gc(j, i) = d;
    end
end
dist_gc(end) = 0;

corr_site = corrcoef(kcs_midday);
figure();
scatter(dist_gc(:), corr_site(:));
xlabel('Distance (km)');
ylabel('Linear correlation coefficient');

% Plot out scatter+histograms between any two pairs of sites
% for i = 1: size(kcs_midday, 2) - 1
%     lat1 = SiteLatitude(i);
%     lon1 = SiteLongitude(i);
%     for j = i+1: size(kcs_midday, 2)
%         figure();
%         lat2 = SiteLatitude(i);
%         lon2 = SiteLongitude(i);
%         scatterhist(kcs_midday(:, i), kcs_midday(:, j));
%         d = distance('gc', lat1, lon1, lat2, lon2);
%     end
% end

end

function conv_corr_ez_nd()
addpath('ndhist\histcn'); % Add the path to the N-Dimensional histcount function
load caiso;
all_years = [2016; 2017; 2018; 2019];
bin_width = 200; % MW

% Let's use GW for better visual effects
bin_width = bin_width./1E3;
for i = 1: 4
    btm{i} = btm{i}./1000;
    loads{i} = loads{i}./1000;
    wind{i} = wind{i}./1000;
    if i <= 3
        pv{i} = pv{i}./1E3;
    end
    if i <= 2
        st{i} = st{i}./1E3;
    end
    stpv2019 = stpv2019./1E3;
end

year_index = 2; % Change to 4 if want to run all four years
switch year_index
    case 1 % 2016
        ts = [loads{year_index}, -wind{year_index}, -st{year_index}, -btm{year_index}, -pv{year_index}];
    case 2 % 2017
        ts = [loads{year_index}, -wind{year_index}, -st{year_index}, -btm{year_index}, -pv{year_index}];
    case 3 % 2018
        ts = [loads{year_index}, -wind{year_index}, -btm{year_index}, -pv{year_index}];
    case 4 % 2019
        ts = [loads{year_index}, -wind{year_index}, -btm{year_index}, -stpv2019];
end

N = size(ts, 2); % Size of components
ts_cdf = nan(size(ts, 1), N); % Time series CDF
gr_pdf = cell(N, 1); % Probability density over each grid
gr_p   = cell(N, 1); % Probability of each grid;
gr_bincenter = cell(N, 1); % Bin center of each grid
gr_edges = cell(N, 1); % Bin edges of each grid, length + 1
gr_cdf = cell(N, 1); % Cumulative probability at the edges of grids

for i = 1:N
    [tmp_pdf, tmp_bincenter, tmp_edges] = return_pdf(ts(:, i), bin_width);
    gr_pdf{i} = tmp_pdf(:);
    gr_bincenter{i} = tmp_bincenter(:);
    gr_edges{i} = tmp_edges(:);
    gr_p{i} = gr_pdf{i}.*bin_width;
    gr_cdf{i} = [0; cumsum(gr_p{i})];
    gr_cdf{i}(end) = 1;
    ts_cdf(:, i) = interp1(gr_edges{i}, gr_cdf{i}, ts(:, i));
end

% Obtain the mean copula pdf (empirical copula)
ts_cdf(ts_cdf==1) = 1-eps; % 1 is the upper limit of the last bin, and if any data point falls on that boundary, histcn will add one more bin.

switch N
    case 4
        count = histcn(ts_cdf, gr_cdf{1}, gr_cdf{2}, gr_cdf{3}, gr_cdf{4});
        [X1, X2, X3, X4] = ndgrid(gr_pdf{1}.*bin_width, gr_pdf{2}.*bin_width, gr_pdf{3}.*bin_width, gr_pdf{4}.*bin_width); % Grid size along each dimension
        grid_vol = X1.*X2.*X3.*X4; % Grid volume
        clear X1 X2 X3 X4; % Save some memory
    case 5
        count = histcn(ts_cdf, gr_cdf{1}, gr_cdf{2}, gr_cdf{3}, gr_cdf{4}, gr_cdf{5});
        [X1, X2, X3, X4, X5] = ndgrid(gr_pdf{1}.*bin_width, gr_pdf{2}.*bin_width, gr_pdf{3}.*bin_width, gr_pdf{4}.*bin_width, gr_pdf{5}.*bin_width); % Grid size along each dimension
        grid_vol = X1.*X2.*X3.*X4.*X5; % Grid volume
        clear X1 X2 X3 X4 X5; % Save some memory
end

copula_pdf_nd = count./(sum(count(:)).*grid_vol); % This is the n-dim empirical copula
copula_pdf_nd(grid_vol==0) = 0; % Sometimes the grid volumn is 0 (when the pdf along any of the grid is 0).

rhohat_gaussian = copulafit('Gaussian', ts_cdf);
u = copularnd('Gaussian',rhohat_gaussian, 1E7); % Monte Carlo simulation
count_gaussian = histcn(u, gr_cdf{1}, gr_cdf{2}, gr_cdf{3}, gr_cdf{4}, gr_cdf{5});
copula_pdf_nd_gaussian = count_gaussian./(sum(count_gaussian(:)).*grid_vol);

% [rhohat_t, param2] = copulafit('t', ts_cdf);
% u = copularnd('t',rhohat_t, param2, 1E7);
% count_t = histcn(u, gr_cdf{1}, gr_cdf{2}, gr_cdf{3}, gr_cdf{4}, gr_cdf{5});
% copula_pdf_nd_t = count_t./(sum(count_t(:)).*grid_vol);

clear count count_gaussian grid_vol u;

[x_conv_corr, p_conv_corr] = conv_ez_corr(gr_bincenter, gr_p, bin_width, copula_pdf_nd);
[x_conv_corr_gaussian, p_conv_corr_gaussian] = conv_ez_corr(gr_bincenter, gr_p, bin_width, copula_pdf_nd_gaussian);
% [x_conv_corr_t, p_conv_corr_t] = conv_ez_corr(gr_bincenter, gr_p, bin_width, copula_pdf_nd_t);

% Actual distribution
e_conv = [x_conv_corr(1) - bin_width/2; x_conv_corr + bin_width/2]; % Bin edges of convolved results
count_actual = histcounts(sum(ts, 2), e_conv); % Actual counts

% Convolution, no correlation
for i = 2: N
    if i == 2
        p_conv = conv(gr_p{1}, gr_p{2});
    else
        p_conv = conv(p_conv, gr_p{i});
    end
end

% Chi^2
chi2_conv = sum((count_actual(:)-p_conv(:).*sum(count_actual(:))).^2./(p_conv(:).*sum(count_actual(:))), 'omitnan');
chi2_conv_corr = sum((count_actual(:)-p_conv_corr(:).*sum(count_actual(:))).^2./(p_conv_corr(:).*sum(count_actual(:))), 'omitnan');
chi2_conv_gaussian = sum((count_actual(:)-p_conv_corr_gaussian(:).*sum(count_actual(:))).^2./(p_conv_corr_gaussian(:).*sum(count_actual(:))), 'omitnan');

% PDF plot
figure();
bar(x_conv_corr ,count_actual./sum(count_actual));
hold on;
plot(x_conv_corr, p_conv_corr, 'k');
plot(x_conv_corr_gaussian, p_conv_corr_gaussian, 'r');
% plot(x_conv_corr_t, p_conv_corr_t, 'b');
plot(x_conv_corr, p_conv, 'b');
legend('Actual', 'Non-parametric', 'Parametric', 'No correlation');

% Calculate CDF
cdf_actual = cumsum(count_actual./sum(count_actual));
cdf_corr = cumsum(p_conv_corr);
cdf_gaussian = cumsum(p_conv_corr_gaussian);
cdf_nocorr = cumsum(p_conv);

% PP plot
figure()
plot(cdf_actual, cdf_corr, '-k.');
hold on;
plot(cdf_actual, cdf_gaussian, '-r.')
plot(cdf_actual, cdf_nocorr, '-b.');
plot([0, 1], [0, 1], 'g');
legend('Non-parametric', 'Parametric', 'No correlation', 'Diagonal');

% CoD
mdl_corr = fitlm(cdf_actual, cdf_corr);
mdl_gaussian = fitlm(cdf_actual, cdf_gaussian);
mdl_nocorr = fitlm(cdf_actual, cdf_nocorr);


end

function test_cap_and_floor()
cell_bincenter = cell(2, 1);
cell_bincenter{1} = [0; 1; 2;];
cell_bincenter{2} = [0; 1; 2;];
cell_p = cell(2, 1);
cell_p{1} = [1/3; 1/3; 1/3;];
cell_p{2} = [1/6; 1/2; 1/3;];
bin_width = 1;
copula_pdf_nd = [3 2/3 1/2; 1 1 1; 2 1/3 3/2;];
cell_lu = cell(2, 1);
cell_lu{1} = [0, 1];
cell_lu{2} = [0, 1];

[x_conv, p_conv] = conv_ez_corr(cell_bincenter, cell_p, bin_width, copula_pdf_nd, cell_lu);
end
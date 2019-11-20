function lab_ghi2power()
add_pvlib();

% dirwork = 'C:\Users\bxl180002\git\SF2\IBM\April\ghi_frcst';
% dirwrite = 'C:\Users\bxl180002\git\SF2\IBM\April\power_frcst';
% deltat = 15; % min
% allp = {'p005', 'mean', 'p095'}; % All percentiles, order should follow csv header
% [allgen, cell_frcst] = ghi2power_frcst(dirwork, deltat, allp, dirwrite);
% identify_violations(allgen, cell_frcst, {'ghi_p005', 'ghi_p095'});
% 
% dirwork = 'C:\Users\bxl180002\git\SF2\IBM\May\ghi_frcst';
% dirwrite = 'C:\Users\bxl180002\git\SF2\IBM\May\power_frcst';
% deltat = 15; % min
% allp = {'p005', 'mean', 'p095'}; % All percentiles, order should follow csv header
% [allgen, cell_frcst] = ghi2power_frcst(dirwork, deltat, allp, dirwrite);
% identify_violations(allgen, cell_frcst, {'ghi_p005', 'ghi_p095'});


dirwork = 'C:\Users\bxl180002\git\SF2\IBM\June.more_quantiles.5min\ghi_frcst';
dirwrite = 'C:\Users\bxl180002\git\SF2\IBM\June.more_quantiles.5min\power_frcst';
deltat = 5; % min
allp = {'p005', 'p025', 'p050', 'p075', 'p095', 'mean'}; % All percentiles, order should follow csv header
[allgen, cell_frcst] = ghi2power_frcst(dirwork, deltat, allp);
identify_violations(allgen, cell_frcst, {'ghi_p005', 'ghi_p025', 'ghi_p050', 'ghi_p075', 'ghi_p095'});


% dirwork = 'C:\Users\bxl180002\git\SF2\IBM\May.more_quantiles.5min\ghi_frcst';
% dirwrite = 'C:\Users\bxl180002\git\SF2\IBM\May.more_quantiles.5min\power_frcst';
% deltat = 5; % min
% allp = {'p005', 'p025', 'p050', 'p075', 'p095', 'mean'}; % All percentiles, order should follow csv header
% [allgen, cell_frcst] = ghi2power_frcst(dirwork, deltat, allp);
% identify_violations(allgen, cell_frcst, {'ghi_p005', 'ghi_p025', 'ghi_p050', 'ghi_p075', 'ghi_p095'});

% dirwork = 'C:\Users\bxl180002\git\SF2\IBM\April.more_quantiles.5min\ghi_frcst';
% dirwrite = 'C:\Users\bxl180002\git\SF2\IBM\April.more_quantiles.5min\power_frcst';
% deltat = 5; % min
% allp = {'p005', 'p025', 'p050', 'p075', 'p095', 'mean'}; % All percentiles, order should follow csv header
% [allgen, cell_frcst] = ghi2power_frcst(dirwork, deltat, allp);
% identify_violations(allgen, cell_frcst, {'ghi_p005', 'ghi_p025', 'ghi_p050', 'ghi_p075', 'ghi_p095'});


% ghi2power_frcst_morequantiles(4, 0);
% ghi2power_actual_hourly(4);
end

function [allgen, cell_frcst] = ghi2power_frcst(dirwork, deltat, allp, dirwrite)
if nargin == 3
    write_flag = false;
elseif nargin == 4
    write_flag = true;
end

% if m == 4
% %     dirwork = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_April\power_frcst';
%     dirwork = 'C:\Users\bxl180002\git\SF2\IBM\April\power_frcst'; % IBM's updated forecast
% elseif m == 5
% %     dirwork = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_May\power_frcst';
%     dirwork = 'C:\Users\bxl180002\git\SF2\IBM\May\power_frcst'; % IBM's updated forecast
% end
dirhome = pwd;

allgen  = {'gen55', 'gen56', 'gen57', 'gen58', 'gen59', 'gen60', 'gen61',    'gen62', 'gen63', 'gen64'};
allsite = {'MNCC1', 'STFC1', 'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};
alllat=[34.31,34.12,34.12,34.12,37.41,35.53,35.38,34.31,34.31,35.53];
alllon=[-117.5,-117.94, -117.94, -117.94,-119.74, -118.63, -120.18,-117.5,-117.5,-118.63];

siteforconv = unique(allsite); % We only interested in these sites because the other sites are not equipped with PV gens
% capacity_gen =  [232.44, 107.58, 116.05, 140.4, 151.32, ];
% capacity_site = [372.84, 134.88, 116.05, 338.52, 151.32];
% cap_scaler = capacity_site./capacity_gen;

% allp = {'p005', 'mean', 'p095'}; % All percentiles, order should follow csv header
strghi = cell(size(allp));
strpwr = cell(size(allp));
for ip = 1: length(allp)
    p = allp{ip};
    strghi{ip} = strcat('ghi', '_', p);
    strpwr{ip} = strcat('pwr', '_', p);
end

cell_frcst = cell(length(allgen), 1);
for i = 1: length(allgen)
    g = allgen{i};
    s = allsite{i};
    
    dirghi = fullfile(fileparts(dirwork), 'ghi_frcst');
    cd(dirghi);
    csvname = strcat('IBM_processed_', s, '.csv');
    T = readtable(csvname);
    cd(dirhome);
    T.Properties.VariableNames(end-length(allp)+1: end) = strghi; % We know the last length(allp) columns are the GHI data
    
    tarray = datetime(T.Year, T.Month, T.Day, T.Hour, T.Minute, zeros(size(T, 1), 1), 'TimeZone', 'UTC');
    Location = pvl_makelocationstruct(alllat(strcmp(allgen, g)),alllon(strcmp(allgen, g)));
        
    tarray_formean = repmat(tarray(:)', deltat, 1) - datenum(repmat([deltat:-1:1]'./24/60, 1, numel(tarray))); % All minutes during the past period, can be 15 or 5 min
    tarray_formean = tarray_formean(:);
    Time_formean.UTCOffset = zeros(size(tarray_formean, 1), 1); % Because IBM uses UTC time, so utc offset is zero
    Time_formean.year   = tarray_formean.Year;
    Time_formean.month  = tarray_formean.Month;
    Time_formean.day    = tarray_formean.Day;
    Time_formean.hour   = tarray_formean.Hour;
    Time_formean.minute = tarray_formean.Minute;
    Time_formean.second = tarray_formean.Second;
    [~, ~, AppSunEl_formean, ~] = pvl_ephemeris(Time_formean,Location);
    ghi_cs_formean = pvl_clearsky_haurwitz(90-AppSunEl_formean);
    ghi_cs_mean = mean(reshape(ghi_cs_formean, numel(ghi_cs_formean)/numel(tarray), numel(tarray)), 1)';
    power_cs_formean = ghi_to_ac_power(g, Time_formean.year, Time_formean.month, Time_formean.day, Time_formean.hour, Time_formean.minute, Time_formean.second, ghi_cs_formean);
    power_cs_mean = mean(reshape(power_cs_formean(:, end), numel(power_cs_formean(:, end))/numel(tarray), numel(tarray)), 1)';

    power_mean_allp = zeros(size(tarray, 1), length(allp));
    for ip = 1: length(allp)
        p = allp{ip};
        col = strghi{ip};
        ghi = T.(col);
        k_cs = ghi./ghi_cs_mean; % Clear-sky index
        ghi_formean = reshape( reshape(ghi_cs_formean, numel(ghi_cs_formean)/numel(tarray), numel(tarray))*diag(k_cs), numel(ghi_cs_formean), 1);
        ghi_formean(isnan(ghi_formean)) = 0; % nan means GHI_CS_mean is zero
        power_formean = ghi_to_ac_power(g, Time_formean.year, Time_formean.month, Time_formean.day, Time_formean.hour, Time_formean.minute, Time_formean.second, ghi_formean);
        power_mean = mean(reshape(power_formean(:, end), numel(power_formean(:, end))/numel(tarray), numel(tarray)), 1)';
        power_mean_allp(:, ip) = power_mean;
    end
    Tnew = [T array2table([ghi_cs_mean, power_mean_allp, power_cs_mean], 'VariableNames', [{'ghi_cs'} strpwr {'pwr_cs'}])];
    cell_frcst{i} = Tnew;

    x = (tarray(1)-duration(1, 0, 0)): datenum(1/24/60): (tarray(end)+datenum(1/24));
    xtime.UTCOffset = zeros(size(x, 1), 1);
    xtime.year   = x.Year;
    xtime.month  = x.Month;
    xtime.day    = x.Day;
    xtime.hour   = x.Hour;
    xtime.minute = x.Minute;
    xtime.second = x.Second;
    [~, ~, xAppSunEl, ~] = pvl_ephemeris(xtime, Location);
    xghi_cs = pvl_clearsky_haurwitz(90-xAppSunEl);
    xpwr_cs = ghi_to_ac_power(g, xtime.year(:), xtime.month(:), xtime.day(:), xtime.hour(:), xtime.minute(:), xtime.second(:), xghi_cs);

    figure();
    ax1 = subplot(2, 1, 1);
    plot(x, xghi_cs, '-k');
    hold on; 
    stairs(tarray-duration(0, deltat, 0), Tnew{:, strghi});
    stairs(tarray-duration(0, deltat, 0), Tnew{:, 'ghi_cs'}, 'k');
    hold off;
    ylabel('GHI');

    ax2 = subplot(2, 1, 2);
    plot(x, xpwr_cs(:, end), '-k');
    hold on; 
    stairs(tarray-duration(0, deltat, 0), Tnew{:, strpwr});
    stairs(tarray-duration(0, deltat, 0), Tnew{:, 'pwr_cs'}, 'k');
    hold off;
    ylabel('Power');

    suptitle(strcat(g, '-', s));
    linkaxes([ax1, ax2],'x');


    figure();
    ax1 = subplot(2, 1, 1);
    plot(tarray, Tnew{:, strghi}./Tnew{:, 'ghi_cs'});
    ylim([0, 2]);
    ylabel('Clear-sky index');

    ax2 = subplot(2, 1, 2);
    plot(tarray, Tnew{:, strpwr}./Tnew{:, 'pwr_cs'});
    ylim([0, 2]);
    ylabel('Normalized power');

    suptitle(strcat(g, '-', s));
    linkaxes([ax1, ax2],'x');

end

if write_flag
    cd(dirwrite);
    for i = 1: length(cell_frcst)
        g = allgen{i};
        csvname_write = strcat('frcst_',         g,      '.csv');
        writetable(cell_frcst{i}, csvname_write);
        fprintf('File %s write to %s!\n', csvname_write, dirwrite);
    end
    cd(dirhome);
end

end

function [allgen, cell_actual] = ghi2power_actual_hourly(m, dirwrite)
if nargin == 1
    write_flag = false;
else
    write_flag = true;
end
if m == 4
%     dirwork = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_April\power_frcst';
    dirwork = 'C:\Users\bxl180002\git\SF2\IBM\April\power_frcst'; % IBM's updated forecast
elseif m == 5
%     dirwork = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_May\power_frcst';
    dirwork = 'C:\Users\bxl180002\git\SF2\IBM\May\power_frcst'; % IBM's updated forecast
end
dirhome = pwd;

allgen  = {'gen55', 'gen56', 'gen57', 'gen58', 'gen59', 'gen60', 'gen61',    'gen62', 'gen63', 'gen64'};
allsite = {'MNCC1', 'STFC1', 'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};
alllat=[34.31,34.12,34.12,34.12,37.41,35.53,35.38,34.31,34.31,35.53];
alllon=[-117.5,-117.94, -117.94, -117.94,-119.74, -118.63, -120.18,-117.5,-117.5,-118.63];

siteforconv = unique(allsite); % We only interested in these sites because the other sites are not equipped with PV gens

cell_actual = cell(length(allgen), 1); % Measured GHI data
dirrealghi = fullfile(fileparts(dirwork), 'ghi_actual');

for i = 1: length(allgen)
    g = allgen{i};
    s = allsite{i};
    csvname = strcat('IBM_processed_', s, '.hourly.csv');
    cd(dirrealghi);
    T = readtable(csvname);
    cd(dirhome);
    T.Properties.VariableNames(end) = {'GHI'};
    T = [T array2table(zeros(size(T, 1), 3), 'VariableNames', {'GHI_CS', 'Power', 'Power_CS'})];

    tarray = datetime(T.Year, T.Month, T.Day, T.Hour, T.Minute, zeros(size(T, 1), 1), 'TimeZone', 'UTC');
    Time.UTCOffset = zeros(size(tarray, 1), 1); % Because we use UTC time, so utc offset is zero
    Time.year   = tarray.Year;
    Time.month  = tarray.Month;
    Time.day    = tarray.Day;
    Time.hour   = tarray.Hour;
    Time.minute = tarray.Minute;
    Time.second = tarray.Second;
    Location = pvl_makelocationstruct(alllat(strcmp(allgen, g)),alllon(strcmp(allgen, g)));

    [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
%     T.GHI_CS = pvl_clearsky_haurwitz(90-AppSunEl); % This is instantaneous clear-sky GHI value

    T.GHI(isnan(T.GHI)) = 0;

    tarray_formean = repmat(tarray(:)', 60, 1) - datenum(repmat([60:-1:1]'./24/60, 1, numel(tarray))); % All minutes during the past hour
    tarray_formean = tarray_formean(:);
    Time_formean.UTCOffset = zeros(size(tarray_formean, 1), 1);
    Time_formean.year   = tarray_formean.Year;
    Time_formean.month  = tarray_formean.Month;
    Time_formean.day    = tarray_formean.Day;
    Time_formean.hour   = tarray_formean.Hour;
    Time_formean.minute = tarray_formean.Minute;
    Time_formean.second = tarray_formean.Second;
    [~, ~, AppSunEl_formean, ~] = pvl_ephemeris(Time_formean,Location);
    ghi_cs_formean = pvl_clearsky_haurwitz(90-AppSunEl_formean);
    ghi_cs_mean = mean(reshape(ghi_cs_formean, numel(ghi_cs_formean)/numel(tarray), numel(tarray)), 1)';
    T.GHI_CS = ghi_cs_mean;

    % We assume the clear-sky index for the previous hour is the same.
    k_cs = T.GHI./T.GHI_CS; % Use average clear-sky index, nan means GHI_CS_mean is zero
    ghi_formean = reshape( reshape(ghi_cs_formean, numel(ghi_cs_formean)/numel(tarray), numel(tarray))*diag(k_cs), numel(ghi_cs_formean), 1);
    ghi_formean(isnan(ghi_formean)) = 0; % nan means GHI_CS_mean is zero

    power_formean = ghi_to_ac_power(g, Time_formean.year, Time_formean.month, Time_formean.day, Time_formean.hour, Time_formean.minute, Time_formean.second, ghi_formean);
    power_mean = mean(reshape(power_formean(:, end), numel(power_formean(:, end))/numel(tarray), numel(tarray)), 1)';
    power_cs_formean = ghi_to_ac_power(g, Time_formean.year, Time_formean.month, Time_formean.day, Time_formean.hour, Time_formean.minute, Time_formean.second, ghi_cs_formean);
    power_cs_mean = mean(reshape(power_cs_formean(:, end), numel(power_cs_formean(:, end))/numel(tarray), numel(tarray)), 1)';
    T.Power = power_mean;
    T.Power_CS = power_cs_mean;
%     T.Power = T.Power + power_mean;
%     T.Power_CS = T.Power_CS + power_cs_mean;

%     x = (tarray(1)-duration(1, 0, 0)): datenum(1/24/60): (tarray(end)+datenum(1/24));
%     xtime.UTCOffset = zeros(size(x, 1), 1);
%     xtime.year = x.Year;
%     xtime.month = x.Month;
%     xtime.day = x.Day;
%     xtime.hour = x.Hour;
%     xtime.minute =x.Minute;
%     xtime.second = x.Second;
%     Location = pvl_makelocationstruct(alllat(strcmp(allgen, g)),alllon(strcmp(allgen, g)));
%     [~, ~, xAppSunEl, ~] = pvl_ephemeris(xtime,Location);
%     xghi_cs = pvl_clearsky_haurwitz(90-xAppSunEl);
%     figure();
%     plot(x, xghi_cs, tarray_formean, ghi_formean, '-.', tarray, T.GHI, '+');

    cell_actual{i} = T; % Note this is UTC time
    figure();
    plot(tarray, T.Power, '-b.', tarray, T.Power_CS, '-ro')

end

if write_flag
    cd(dirwrite);
    for i = 1: length(cell_actual)
        g = allgen{i};
        csvname_write = strcat('actual_',         g,      '.csv');
        writetable(cell_actual{i}, csvname_write);
        fprintf('File %s write to %s!\n', csvname_write, dirwrite);
    end
    cd(dirhome);
end

end

function identify_violations(selected_gen, cell_frcst, allp)
% allgen is a cell of gen names, cell_frcst is a cell of tables, two
% variables should be the same length
% allp = {'ghi_p005', 'ghi_p095'};
allgen  = {'gen55', 'gen56', 'gen57', 'gen58', 'gen59', 'gen60', 'gen61',    'gen62', 'gen63', 'gen64'};
allsite = {'MNCC1', 'STFC1', 'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};

compare_ghi_greater   = zeros(numel(selected_gen), length(allp) - 1);
% compare_power_greater = zeros(numel(allgen), length(allp) - 1);
compare_ghi_equal     = zeros(numel(selected_gen), length(allp) - 1);
% compare_power_equal   = zeros(numel(allgen), length(allp) - 1);

title_compare = cell(1, length(allp)-1);
for j = 1: length(allp) - 1
    lp = allp{j};
    rp = allp{j+1};
    title_compare{j} = strcat(lp, '_ge_', rp);
end

for i = 1: length(selected_gen)
    Tnew = cell_frcst{i};
    for j = 1: length(allp) - 1
        lp = allp{j};
        rp = allp{j+1};
        compare_ghi_equal(i, j)   = sum((Tnew.(lp)==Tnew.(rp))&(Tnew.(lp)~=0)&(Tnew.(rp)~=0));
        compare_ghi_greater(i, j) = sum((Tnew.(lp)>Tnew.(rp))&(Tnew.(lp)~=0)&(Tnew.(rp)~=0));
    end
end
fprintf('GHI comparison');
[table(selected_gen(:), 'VariableNames', {'sitename'}) array2table(compare_ghi_greater+compare_ghi_equal, 'VariableNames', title_compare)]

end

function explore_solar_irradiance()
% Explore DNI, DHI, and GHI
clear;
sitenames={'gen55','gen56','gen57','gen58','gen59','gen60','gen61','gen62','gen63','gen64'};
IBMsitenames = {'MNCC1', 'STFC1', 'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};
% dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_old\ghi_frcst'; % Old data
% dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_April\ghi_frcst'; % April data
% dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_May\ghi_frcst'; % May data
dir_work = 'C:\Users\bxl180002\git\SF2\IBM\April\ghi_frcst'; % Updated April data
% dir_work = 'C:\Users\bxl180002\git\SF2\IBM\May\ghi_frcst'; % Updated May data
dir_home = pwd;
SiteLatitude  = [ 34.31,  34.12,   34.12,   34.12,   37.41,   35.53,   35.38,  34.31,  34.31,   35.53];
SiteLongitude = [-117.5,-117.94, -117.94, -117.94, -119.74, -118.63, -120.18, -117.5, -117.5, -118.63];

fprintf('%6s %8s %6s %6s %6s %6s %6s %6s %6s\n', 'gen', 'site', 'L>M(I)', 'M>U(I)', 'L>M(P)', 'M>U(P)', 'kt>1(L)', 'kt>1(M)', 'kt>1(H)');

for k = 1: length(sitenames)
    gen = sitenames{k};
    ibm_site = IBMsitenames{k};
    csvname_read  = strcat('IBM_processed_', ibm_site, '.csv');
    
    cd(dir_work);
    M = csvread(csvname_read, 1, 0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fix 1: Time shift, this is included in fix_ibm, added here for
    % plotting purpose
%     toff=2;
%     M0 = [M(1:size(M, 1)-toff, 1: 5), M(toff+1: end, 6:8)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd(dir_home);
    M0 = M;
%     M = fix_ibm(M, SiteLatitude(k), SiteLongitude(k));
    M(any(isnan(M(:, 6:8)), 2), 6:8) = 0; 
    
    site_name = gen;
    utc_year   = M(:, 1);
    utc_month  = M(:, 2);
    utc_day    = M(:, 3);
    utc_hour   = M(:, 4);
    utc_minute = M(:, 5);
    utc_second = zeros(size(M, 1), 1);

    modulepv=134; %Topaz uses first solar panels (FS272 is probably the oldest they have)
    inverterid=759; 
%     Site_tilt=[0,0,0,25,0,0,25,22.5,0,0];
    Site_tilt=[0,0,0,0,0,0,0,0,0,0]; % Single tracking
    modules_Series=11;
    modules_parallel=[2360,1870,0,2002,2096,2512,2319,2381,0,2353];
    ninvert=[149,97,0,23,82,90,97,90,0,127];

    %Other parameters
    SF=0.98;
    %Weather
    PresPa=101325;
    WIND=0;
    dryT=10;
    Albedo = 0.2;

    % clearvars -except sitenames SiteLatitude SiteLongitude  IBMsitenames k modulepv inverterid Site_tilt modules_Series modules_parallel ninvert SF PresPa WIND dryT Albedo
    Location = pvl_makelocationstruct(SiteLatitude(k),SiteLongitude(k)); %Altitude is optional
    %--- SOLAR FARM SPECS---
    %Define module
    ModuleParameters = pvl_sapmmoduledb(modulepv,'SandiaModuleDatabase_20120925.xlsx');
    %Define the inverter
    load('SandiaInverterDatabaseSAM2014.1.14.mat')
    Inverter = SNLInverterDB(inverterid);
    %Topaz uses power one inverters
    clear InverterNames SNLInverterDB
    %Define the array configuration
    Array.Tilt = Site_tilt(k); % Array tilt angle (deg)
    Array.Azimuth = 180; %Array azimuth (180 deg indicates array faces South)
    Array.Ms = modules_Series; %Number of modules in series
    Array.Mp = modules_parallel(k); %Number of paralell strings  
    %Location of site
    Array.a = -3.56;
    Array.b = -0.075;

    % Prepare for ac power calculation.
    Time.UTCOffset(1:size(M,1),1) = zeros(size(M,1), 1); % Because we use UTC time, so utc offset is zero
    Time.year(1:size(M,1),1)   = utc_year;
    Time.month(1:size(M,1),1)  = utc_month;
    Time.day(1:size(M,1),1)    = utc_day;
    Time.hour(1:size(M,1),1)   = utc_hour;
    Time.minute(1:size(M,1),1) = utc_minute;
    Time.second(1:size(M,1),1) = utc_second;
    ACPower(1:size(M,1),1:6) = M(:,1:6);

    %used for both forecast and actual
    dayofyear = pvl_date2doy(Time.year, Time.month, Time.day);
    [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Fix 2: GHI capped by GHI_cs
    M0 = M;
    ghi_clearsky = pvl_clearsky_haurwitz(90-AppSunEl); % Clear-sky GHI
%     fixrow = find(any(M(:, 6:8) > repmat(ghi_clearsky, 1, 3), 2));
%     for i = 1:length(fixrow)
%         ifixrow = fixrow(i);
%         ghi_max = max(M(ifixrow, 6:8));
%         M(ifixrow, 6:8) = M(ifixrow, 6:8).*ghi_clearsky(ifixrow)/ghi_max;
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if Site_tilt(k)==0    
        [TrkrTheta, AOI, SurfTilt, SurfAz] = pvl_singleaxis(90-AppSunEl, SunAz, Location.latitude, 0, 180, 45);
        Array.Tilt=SurfTilt;
        Array.Tilt(find(isnan(Array.Tilt)))=0;
    else
        AOI = pvl_getaoi(Array.Tilt, Array.Azimuth, 90-AppSunEl, SunAz);
    end
    Wspd=repmat(WIND,size(M,1));
    Drybulb=repmat(dryT,size(M,1));
    AMa = pvl_absoluteairmass(pvl_relativeairmass(90-AppSunEl),PresPa);
    F1 = max(0,polyval(ModuleParameters.a,AMa)); %Spectral loss function
    F2 = max(0,polyval(ModuleParameters.b,AOI)); % Angle of incidence loss function

    compare_result = zeros(1, 7);
    compare_result(1) = sum(M(:, 6)>M(:, 7));
    compare_result(2) = sum(M(:, 7)>M(:, 8));

    fig4in1 = figure();
    for j = 6: 8
        ghi = M(:, j);
    %     ghi = ReGHI(:, end);
        EdiffGround = pvl_grounddiffuse(Array.Tilt, ghi, Albedo);
        DNI_model = pvl_disc(ghi,90-SunEl, dayofyear,PresPa);
%         DNI_model = pvl_erbs(ghi, 90-SunEl, dayofyear);
        DHI_model = ghi - cosd(90-SunEl).*DNI_model;
        Eb = 0*AOI; %Initiallize variable
        Eb(AOI<90) = DNI_model(AOI<90).*cosd(AOI(AOI<90)); %Only calculate when sun is in view of the plane of array
        EdiffSky = pvl_isotropicsky(Array.Tilt,DHI_model);
        E = Eb + EdiffSky + EdiffGround; % Total incident irradiance (W/m^2)
        E0 = 1000; %Reference irradiance (1000 W/m^2)
        celltemp = pvl_sapmcelltemp(E, E0, Array.a, Array.b,Wspd(:,1), Drybulb(:,1), ModuleParameters.delT);
        Ediff = EdiffSky + EdiffGround; % Total diffuse incident irradiance (W/m^2)
        Ee = F1.*((Eb.*F2+ModuleParameters.fd.*Ediff)/E0)*SF; %Effective irradiance
        Ee(isnan(Ee))=0; % Set any NaNs to zero
        mSAPMResults = pvl_sapm(ModuleParameters, Ee, celltemp);
        aSAPMResults.Vmp = Array.Ms  *mSAPMResults.Vmp;
        aSAPMResults.Imp = Array.Mp  *mSAPMResults.Imp;
        aSAPMResults.Pmp = aSAPMResults.Vmp .* aSAPMResults.Imp;
        clear temp
        temp=find(ghi~=10000);
        ACPower(1:size(Time.hour,1), 7) =10000;
        ACPower(temp, 7)= pvl_snlinverter(Inverter, mSAPMResults.Vmp(temp)*Array.Ms, mSAPMResults.Pmp(temp)*Array.Ms*Array.Mp)*ninvert(k)/1000000;
        ACPower(ACPower(:,7)<0, end)=0;

        tarray = datetime(Time.year, Time.month, Time.day, Time.hour, Time.minute, Time.second, 'TimeZone', 'UTC');
        tarray.TimeZone = 'America/Los_Angeles';
        tarray_local = datetime(tarray.Year, tarray.Month, tarray.Day, tarray.Hour, tarray.Minute, tarray.Second);
        
        DayAngle = 2.*pi.*(dayofyear-1)./365;
        re = 1.00011 + 0.034221 .* cos(DayAngle) + (0.00128) .* sin(DayAngle)...
        +0.000719.*cos(2.*DayAngle) + (7.7E-5).*sin(2.*DayAngle);
        I0 = re.*1370;
        I0h= I0.*cosd(90-SunEl);
        I0h(I0h<0)=0;
        clearness_index = ghi./I0h;

        figure();
        ax1 = subplot(2, 1, 1);
        h2 = plot(tarray_local, ghi);
        ylabel('GHI');
        hold on;
%         h3 = plot(tarray_local, I0h);
        h3 = plot(tarray_local, ghi_clearsky);
        h5 = plot(tarray_local, M0(:, j));
        yyaxis right;
        h4 = plot(tarray_local, SunEl);
        ylabel('Elevation angle');
        legend([h2, h3, h4, h5], {'GHI', 'GHI (CS)', 'Elevation angle', 'GHI0'});
%         legend([h2, h3, h4, h5], {'GHI', 'GHI (E)', 'Elevation angle', 'GHI0'});
        
        ax2 = subplot(2, 1, 2);
        h1 = plot(tarray_local, clearness_index);
        ylabel('Clearness index (kt)');
%         legend([h1, h2, h3], {'Clearness index', 'GHI', 'I0h'});
        suptitle(strcat(gen, ',', ibm_site));

        linkaxes([ax1, ax2],'x');

%         switch j
%             case 6
%                 ax_l = subplot(4, 1, 1);
%                 ylabel('Lower');
%                 p_low = ACPower(:, end);
%             case 7
%                 ax_m = subplot(4, 1, 2);
%                 ylabel('Mean');
%                 p_mean = ACPower(:, end);
%             case 8
%                 ax_h = subplot(4, 1, 3);
%                 ylabel('Upper');
%                 p_up = ACPower(:, end);
%         end
%         plot(tarray_local, DNI_model, '--k', tarray_local, DHI_model, '-.k', tarray_local, ghi, '-k');
%         legend('DNI', 'DHI', 'GHI');
        
        switch j
            case 6
                p_low = ACPower(:, end);
                str_percentile = '5';
                linestyle = '-b';
            case 7
                p_mean = ACPower(:, end);
                str_percentile = '50';
                linestyle = '-k';
            case 8
                p_up = ACPower(:, end);
                str_percentile = '95';
                linestyle = '-r';
        end
        
        figure(fig4in1);
        ax_l = subplot(4, 1, 1);
        plot(tarray_local, DNI_model, linestyle);
        hold on;
        ylabel('DNI');

        ax_m = subplot(4, 1, 2);
        plot(tarray_local, DHI_model, linestyle);
        hold on;
        ylabel('DHI');

        ax_h = subplot(4, 1, 3);
        plot(tarray_local, ghi, linestyle);
        plot(tarray_local, M0(:, j), '--k');
        hold on;
        ylabel('GHI');
        
        switch j
            case 6
                compare_result(5) = sum(clearness_index>1);
            case 7
                compare_result(6) = sum(clearness_index>1);
            case 8 
                compare_result(7) = sum(clearness_index>1);
        end
        
    end

    compare_result(3) = sum(p_low>p_mean);
    compare_result(4) = sum(p_mean>p_up);

    fprintf('%6s %8s %6g %6g %6g %6g %6g %6g %6g\n', gen, ibm_site, compare_result(1), compare_result(2), compare_result(3), compare_result(4), compare_result(5), compare_result(6), compare_result(7));
    figure(fig4in1);
    suptitle(strcat(gen, ',', ibm_site));

    ax_p = subplot(4, 1, 4);
    plot(tarray_local, p_low, 'b', tarray_local, p_mean, 'k', tarray_local, p_up, 'r');
    ylabel('Power');
    legend('Low', 'Mean', 'High');
    linkaxes([ax_l, ax_m, ax_h, ax_p],'x');
    clear M M0;
        
end
end

function find_violations()
% Find out places where lower percentiles are higher than higher percentiles
clear;
IBMsitenames  = {'CA_Topaz','COWC1','DEMC1','KNNC1','MIAC1','MNCC1','RLKC1','RSAC1','SBVC1','STFC1'};
dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_old\ghi_frcst';
% dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_May\ghi_frcst';
% dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_April\ghi_frcst';
dir_home = pwd;

for k = 1:length(IBMsitenames)
    ibm_site = IBMsitenames{k};
    csvname_read  = strcat('IBM_processed_', ibm_site, '.csv');
    csvname_write = strcat('Errors_', ibm_site, '.csv');
    
    cd(dir_work);
    M = csvread(csvname_read, 1, 0);
    cd(dir_home);
    tarray = datetime( M(:, 1), M(:, 2), M(:, 3), M(:, 4), M(:, 5), zeros(size(M, 1), 1), 'TimeZone', 'UTC');
    
    cl_05 = M(:, 6);
    cl_50 = M(:, 7);
    cl_95 = M(:, 8);
    
    fprintf('site %s, 5-p > mean: %3g, mean > 95-p: %3g\n', ibm_site, sum(cl_05>cl_50), sum(cl_50>cl_95));
    idx = find((cl_05>cl_50) | (cl_50>cl_95));
    for i = 1: length(idx)
        isample = idx(i);
        fprintf('%s %6.2f  %6.2f  %6.2f\n', datestr(tarray(isample),'yyyy-mm-dd HH:MM'), cl_05(isample), cl_50(isample), cl_95(isample));
    end
    
%     cd(dir_work);
%     cHeader = {'Year' 'Month' 'Day' 'Hour' 'Minute' '5-p' 'Mean' '95-p'}; %dummy header
%     commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
%     commaHeader = commaHeader(:)';
%     textHeader = cell2mat(commaHeader); % cHeader in text with commas
%     fid = fopen(csvname_write,'w'); 
%     fprintf(fid,'%s\n',textHeader); % write header to file
%     fclose(fid);
%     dlmwrite(csvname_write,[M(idx, :)],'-append');
%     cd(dir_home);

end
end

function explore_ghi_vs_power()


gen = 'gen64';
PresPa=101325;

sitenames    = {'gen55', 'gen56', 'gen57', 'gen58', 'gen59', 'gen60', 'gen61',    'gen62', 'gen63', 'gen64'};
SiteLatitude=[34.31,34.12,34.12,34.12,37.41,35.53,35.38,34.31,34.31,35.53];
SiteLongitude=[-117.5,-117.94, -117.94, -117.94,-119.74, -118.63, -120.18,-117.5,-117.5,-118.63];

% figure; 
for h = -10:2
%     color_code = rand(1, 3);
    color_code = 'r';
    utc_year = 2019.*ones(1300, 1);
    utc_month = 5.*ones(1300, 1);
    utc_day = 2.*ones(1300, 1);
%     utc_hour = 2.*ones(1300, 1);
    utc_hour = h.*ones(1300, 1);
    utc_minute = 0.*ones(1300, 1);
    utc_second = zeros(1300, 1);

    Time.UTCOffset = zeros(size(1300,1), 1); % Because we use UTC time, so utc offset is zero
    Time.year   = utc_year;
    Time.month  = utc_month;
    Time.day    = utc_day;
    Time.hour   = utc_hour;
    Time.minute = utc_minute;
    Time.second = utc_second;
    dayofyear = pvl_date2doy(Time.year, Time.month, Time.day);

    k = strcmp(sitenames, gen);
    Location = pvl_makelocationstruct(SiteLatitude(k),SiteLongitude(k));

    [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);


    DayAngle = 2.*pi.*(dayofyear-1)./365;
    re = 1.00011 + 0.034221 .* cos(DayAngle) + (0.00128) .* sin(DayAngle)...
         +0.000719.*cos(2.*DayAngle) + (7.7E-5).*sin(2.*DayAngle);
    I0 = re.*1370;
    I0h= I0.*cosd(90-SunEl); % Extraterrestrial GHI
    ghi_clearsky = pvl_clearsky_haurwitz(90-AppSunEl); % Clear-sky GHI

    ghi_e = unique(I0h);
    ghi_cs = unique(ghi_clearsky);

    ghi = [1:round(ghi_e*1.5)]';
    n = min(size(utc_year, 1), size(ghi, 1));
    DNI_model = pvl_disc(ghi(1:n),90-SunEl(1:n), dayofyear(1:n),PresPa);
%     DNI_model = pvl_erbs(ghi, 90-SunEl(1:n), dayofyear(1:n));
    ACPower = ghi_to_ac_power(gen, utc_year(1:n), utc_month(1:n), utc_day(1:n), utc_hour(1:n), utc_minute(1:n), utc_second(1:n), ghi(1:n));

    figure();
    hax=axes; 
    plot(ghi(1:n), ACPower(:, end),'Color', 'b');
    line([ghi_e,  ghi_e], get(hax,'YLim'),'Color','r');
    line([ghi_cs, ghi_cs], get(hax,'YLim'),'Color','r');
    xlabel('GHI (W/m^2)');
    ylabel('Power (W)');
    title(strcat(gen, ', ', int2str(mod(h-7, 24)), ':00'));
%     yyaxis right;
%     plot(ghi, DNI_model);
%     ylabel('DNI (W/m^s)');
%     hold on;
end

end
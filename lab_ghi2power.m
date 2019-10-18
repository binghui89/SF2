function lab_ghi2power()
add_pvlib();

ghi2power_frcst_5min(4, 0);
end

function add_pvlib()
% Add path and Site ID mapping
addpath('C:\Users\bxl180002\git\MATLAB_PV_LIB');
addpath('C:\Users\bxl180002\git\MATLAB_PV_LIB\Example Data');
addpath('C:\Users\bxl180002\git\MATLAB_PV_LIB\Required Data');
end

function ghi2power_frcst_5min(m, write_flag)
% Convert GHI into power using Elina's script: 5-min forecast
if nargin==1
    write_flag = 0;
end

sitenames     = {'gen55', 'gen56', 'gen57', 'gen58', 'gen59', 'gen60', 'gen61',    'gen62', 'gen63', 'gen64'};
IBMsitenames  = {'MNCC1', 'STFC1', 'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};
SiteLatitude  = [34.31,   34.12,   34.12,   34.12,   37.41,   35.53,   35.38,  34.31, 34.31,    35.53];
SiteLongitude = [-117.5,-117.94, -117.94, -117.94, -119.74, -118.63, -120.18, -117.5, -117.5, -118.63];

if m == 5
    dir_work = 'C:\Users\bxl180002\git\SF2\IBM\May.more_quantiles.5min\ghi_frcst'; % 5-min May data
elseif m == 4
    dir_work = 'C:\Users\bxl180002\git\SF2\IBM\April.more_quantiles.5min\ghi_frcst'; % 5-min April data
end

dir_home = pwd;
fprintf('%6s %8s %6s %6s %6s %6s\n', 'gen', 'site', 'L>M(I)', 'M>U(I)', 'L>M(P)', 'M>U(P)');

compare_ghi_greater   = zeros(numel(sitenames), 4);
compare_power_greater = zeros(numel(sitenames), 4);
compare_ghi_equal   = zeros(numel(sitenames), 4);
compare_power_equal = zeros(numel(sitenames), 4);

for k = 1:numel(sitenames)
    gen = sitenames{k};
    ibm_site = IBMsitenames{k};
    csvname_read  = strcat('IBM_processed_', ibm_site, '.csv');
    csvname_write = strcat('power_',         gen,      '.csv');
    
    cd(dir_work);
%     M = csvread(csvname_read, 1, 0);
    T = readtable(csvname_read);
    cell_cols = {'p5', 'p25', 'p50', 'p75', 'p95', 'Mean'}; % This is the names of table columns
    cd(dir_home);
    
%     M = fix_ibm(M, SiteLatitude(k),SiteLongitude(k));
        
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Fix 1, time shift back by 30 min (2 x 15 min interval)
%     toff = 2;
%     Mp = [M(1:size(M, 1) - toff, 1: 5), M(toff + 1: end, 6:8)];
%     M = Mp;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear Time;
    Location = pvl_makelocationstruct(SiteLatitude(k),SiteLongitude(k)); %Altitude is optional
    Time.UTCOffset(1:size(T,1),1) = zeros(size(T,1), 1); % Because we use UTC time, so utc offset is zero
    Time.year(1:size(T,1),1)   = T.Year;
    Time.month(1:size(T,1),1)  = T.Month;
    Time.day(1:size(T,1),1)    = T.Day;
    Time.hour(1:size(T,1),1)   = T.Hour;
    Time.minute(1:size(T,1),1) = T.Minute;
    Time.second(1:size(T,1),1) = zeros(size(T,1), 1);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Fix 2, capped by clear-sky GHI
%     % Prepare for ac power calculation.
%     [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
%     ghi_clearsky = pvl_clearsky_haurwitz(90-AppSunEl); % Clear-sky GHI
%     fixrow = find(any(M(:, 6:8) > repmat(ghi_clearsky, 1, 3), 2));
%     for i = 1:length(fixrow)
%         ifixrow = fixrow(i);
%         ghi_max = max(M(ifixrow, 6:8));
%         M(ifixrow, 6:8) = M(ifixrow, 6:8).*ghi_clearsky(ifixrow)/ghi_max;
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    power_percentiles = nan(size(T, 1), numel(cell_cols));
    for i = 1: numel(cell_cols)
        thiscol = cell_cols{i};
        ac = ghi_to_ac_power(gen, T.Year, T.Month, T.Day, T.Hour, T.Minute, zeros(size(T, 1), 1), T{:, thiscol});
        power_percentiles(:, i) = ac(:, end);
    end

    tarray = datetime(Time.year, Time.month, Time.day, Time.hour, Time.minute, Time.second, 'TimeZone', 'UTC');
    tarray.TimeZone = 'America/Los_Angeles';
    tarray_local = datetime(tarray.Year, tarray.Month, tarray.Day, tarray.Hour, tarray.Minute, tarray.Second);
    array_power = [tarray_local.Year, tarray_local.Month, tarray_local.Day, tarray_local.Hour, tarray_local.Minute, power_percentiles];
    Tpower = array2table(array_power, 'VariableNames', T.Properties.VariableNames);
    
    compare_ghi_equal(k, 1) = sum((T.p5==T.p25)&(T.p5~=0)&(T.p25~=0));
    compare_ghi_equal(k, 2) = sum((T.p25==T.p50)&(T.p25~=0)&(T.p50~=0));
    compare_ghi_equal(k, 3) = sum((T.p50==T.p75)&(T.p50~=0)&(T.p75~=0));
    compare_ghi_equal(k, 4) = sum((T.p75==T.p95)&(T.p75~=0)&(T.p95~=0));
    logical_ghi_equal = ((T.p5==T.p25)&(T.p5~=0)&(T.p25~=0) | (T.p25==T.p50)&(T.p25~=0)&(T.p50~=0) | (T.p50==T.p75)&(T.p50~=0)&(T.p75~=0) | (T.p75==T.p95)&(T.p75~=0)&(T.p95~=0));
    csvname_ghi_equal = strcat('ghi.', int2str(m), '.', ibm_site, '.', 'equal.csv');
    if( (size(T(logical_ghi_equal, :), 1) > 0) && ~isfile(csvname_ghi_equal) )
        writetable(T(logical_ghi_equal, :), csvname_ghi_equal);
    end
    
    compare_power_equal(k, 1) = sum((Tpower.p5==Tpower.p25)&(Tpower.p5~=0)&(Tpower.p25~=0));
    compare_power_equal(k, 2) = sum((Tpower.p25==Tpower.p50)&(Tpower.p25~=0)&(Tpower.p50~=0));
    compare_power_equal(k, 3) = sum((Tpower.p50==Tpower.p75)&(Tpower.p50~=0)&(Tpower.p75~=0));
    compare_power_equal(k, 4) = sum((Tpower.p75==Tpower.p95)&(Tpower.p75~=0)&(Tpower.p95~=0));

    compare_ghi_greater(k, 1) = sum((T.p5>T.p25)&(T.p5~=0)&(T.p25~=0));
    compare_ghi_greater(k, 2) = sum((T.p25>T.p50)&(T.p25~=0)&(T.p50~=0));
    compare_ghi_greater(k, 3) = sum((T.p50>T.p75)&(T.p50~=0)&(T.p75~=0));
    compare_ghi_greater(k, 4) = sum((T.p75>T.p95)&(T.p75~=0)&(T.p95~=0));
    logical_ghi_greater = ((T.p5>T.p25)&(T.p5~=0)&(T.p25~=0) | (T.p25>T.p50)&(T.p25~=0)&(T.p50~=0) | (T.p50>T.p75)&(T.p50~=0)&(T.p75~=0) | (T.p75>T.p95)&(T.p75~=0)&(T.p95~=0));
    csvname_ghi_greater = strcat('ghi.', int2str(m), '.', ibm_site, '.', 'greater.csv');
    if( (size(T(logical_ghi_greater, :), 1) > 0) && ~isfile(csvname_ghi_greater) )
        writetable(T(logical_ghi_greater, :), csvname_ghi_greater);
    end

    compare_power_greater(k, 1) = sum((Tpower.p5>Tpower.p25)&(Tpower.p5~=0)&(Tpower.p25~=0));
    compare_power_greater(k, 2) = sum((Tpower.p25>Tpower.p50)&(Tpower.p25~=0)&(Tpower.p50~=0));
    compare_power_greater(k, 3) = sum((Tpower.p50>Tpower.p75)&(Tpower.p50~=0)&(Tpower.p75~=0));
    compare_power_greater(k, 4) = sum((Tpower.p75>Tpower.p95)&(Tpower.p75~=0)&(Tpower.p95~=0));

    if write_flag
        cd(dir_work);
        writetable(Tpower, csvname_write);
        cd(dir_home);
    end
 
%     figure();
%     h = plot(tarray_local, [ac_lb(:, end), ac_mean(:, end), ac_ub(:, end)]);
%     set(h, {'color'}, {'b'; 'k'; 'r'});
%     title(strcat(gen, ',', ibm_site));
%     ylabel('kW');
%     fprintf('%6s %8s %6g %6g %6g %6g\n', gen, ibm_site, compare_result(1), compare_result(2), compare_result(3), compare_result(4));

    
end

fprintf('GHI comparison');
[table(sitenames', IBMsitenames', 'VariableNames', {'sitename', 'ibmname'}) array2table(compare_ghi_greater+compare_ghi_equal, 'VariableNames', {'p5_ge_p25', 'p25_ge_p50', 'p50_ge_p75', 'p75_ge_p95'})]
fprintf('Power comparison');
[table(sitenames', IBMsitenames', 'VariableNames', {'sitename', 'ibmname'}) array2table(compare_power_greater+compare_power_equal, 'VariableNames', {'p5_ge_p25', 'p25_ge_p50', 'p50_ge_p75', 'p75_ge_p95'})]
end

function ghi2power_frcst_15min()
% Convert GHI into power using Elina's script: 15-min forecast
clear;
sitenames     = {'gen55', 'gen56', 'gen57', 'gen58', 'gen59', 'gen60', 'gen61',    'gen62', 'gen63', 'gen64'};
IBMsitenames  = {'MNCC1', 'STFC1', 'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};
SiteLatitude  = [34.31,   34.12,   34.12,   34.12,   37.41,   35.53,   35.38,  34.31, 34.31,    35.53];
SiteLongitude = [-117.5,-117.94, -117.94, -117.94, -119.74, -118.63, -120.18, -117.5, -117.5, -118.63];
% dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_old\ghi_frcst';
% dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_April\ghi_frcst';
% dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_May\ghi_frcst';
% dir_work = 'C:\Users\bxl180002\git\SF2\IBM\May\ghi_frcst'; % Updated May data
dir_work = 'C:\Users\bxl180002\git\SF2\IBM\April\ghi_frcst'; % Updated April data

dir_home = pwd;
fprintf('%6s %8s %6s %6s %6s %6s\n', 'gen', 'site', 'L>M(I)', 'M>U(I)', 'L>M(P)', 'M>U(P)');

for k = 1:length(sitenames)
    gen = sitenames{k};
    ibm_site = IBMsitenames{k};
    csvname_read  = strcat('IBM_processed_', ibm_site, '.csv');
    csvname_write = strcat('power_',         gen,      '.csv');
    
    cd(dir_work);
    M = csvread(csvname_read, 1, 0);
    cd(dir_home);
    
%     M = fix_ibm(M, SiteLatitude(k),SiteLongitude(k));
        
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Fix 1, time shift back by 30 min (2 x 15 min interval)
%     toff = 2;
%     Mp = [M(1:size(M, 1) - toff, 1: 5), M(toff + 1: end, 6:8)];
%     M = Mp;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear Time;
    Location = pvl_makelocationstruct(SiteLatitude(k),SiteLongitude(k)); %Altitude is optional
    Time.UTCOffset(1:size(M,1),1) = zeros(size(M,1), 1); % Because we use UTC time, so utc offset is zero
    Time.year(1:size(M,1),1)   = M(:, 1);
    Time.month(1:size(M,1),1)  = M(:, 2);
    Time.day(1:size(M,1),1)    = M(:, 3);
    Time.hour(1:size(M,1),1)   = M(:, 4);
    Time.minute(1:size(M,1),1) = M(:, 5);
    Time.second(1:size(M,1),1) = zeros(size(M,1), 1);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Fix 2, capped by clear-sky GHI
%     % Prepare for ac power calculation.
%     [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
%     ghi_clearsky = pvl_clearsky_haurwitz(90-AppSunEl); % Clear-sky GHI
%     fixrow = find(any(M(:, 6:8) > repmat(ghi_clearsky, 1, 3), 2));
%     for i = 1:length(fixrow)
%         ifixrow = fixrow(i);
%         ghi_max = max(M(ifixrow, 6:8));
%         M(ifixrow, 6:8) = M(ifixrow, 6:8).*ghi_clearsky(ifixrow)/ghi_max;
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    ac_lb   = ghi_to_ac_power(gen, M(:, 1), M(:, 2), M(:, 3), M(:, 4), M(:, 5), zeros(size(M, 1), 1), M(:, 6));
    ac_mean = ghi_to_ac_power(gen, M(:, 1), M(:, 2), M(:, 3), M(:, 4), M(:, 5), zeros(size(M, 1), 1), M(:, 7));
    ac_ub   = ghi_to_ac_power(gen, M(:, 1), M(:, 2), M(:, 3), M(:, 4), M(:, 5), zeros(size(M, 1), 1), M(:, 8));
    
    compare_result = zeros(1, 4);
    compare_result(1) = sum(M(:, 6)>M(:, 7));
    compare_result(2) = sum(M(:, 7)>M(:, 8));
    compare_result(3) = sum(ac_lb(:, end) > ac_mean(:, end));
    compare_result(4) = sum(ac_mean(:, end) > ac_ub(:, end));
    
    tarray = datetime(Time.year, Time.month, Time.day, Time.hour, Time.minute, Time.second, 'TimeZone', 'UTC');
    tarray.TimeZone = 'America/Los_Angeles';
    tarray_local = datetime(tarray.Year, tarray.Month, tarray.Day, tarray.Hour, tarray.Minute, tarray.Second);

    figure();
    h = plot(tarray_local, [ac_lb(:, end), ac_mean(:, end), ac_ub(:, end)]);
    set(h, {'color'}, {'b'; 'k'; 'r'});
    title(strcat(gen, ',', ibm_site));
    ylabel('kW');
    fprintf('%6s %8s %6g %6g %6g %6g\n', gen, ibm_site, compare_result(1), compare_result(2), compare_result(3), compare_result(4));
    
    ac_lb(:, 1:6)   = [tarray_local.Year, tarray_local.Month, tarray_local.Day, tarray_local.Hour, tarray_local.Minute, tarray_local.Second];
    ac_mean(:, 1:6) = [tarray_local.Year, tarray_local.Month, tarray_local.Day, tarray_local.Hour, tarray_local.Minute, tarray_local.Second];
    ac_ub(:, 1:6)   = [tarray_local.Year, tarray_local.Month, tarray_local.Day, tarray_local.Hour, tarray_local.Minute, tarray_local.Second];

    cd(dir_work);
    cHeader = {'Year' 'Month' 'Day' 'Hour' 'Minute' 'Second' 'Low' 'Mean' 'High'}; %dummy header
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); % cHeader in text with commas
    fid = fopen(csvname_write,'w'); 
    fprintf(fid,'%s\n',textHeader); % write header to file
    fclose(fid);
    dlmwrite(csvname_write,[ac_lb(:, 1:6), ac_lb(:, end), ac_mean(:, end), ac_ub(:, end)],'-append');
    cd(dir_home);

    
end

end

function ghi2power_actual_15min()
% Convert GHI into power using Elina's script: actual
clear;
sitenames     = {'gen55', 'gen56', 'gen57', 'gen58', 'gen59', 'gen60', 'gen61',    'gen62', 'gen63', 'gen64'};
IBMsitenames  = {'MNCC1', 'STFC1', 'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};
SiteLatitude  = [34.31,   34.12,   34.12,   34.12,   37.41,   35.53,   35.38,  34.31, 34.31,    35.53];
SiteLongitude = [-117.5,-117.94, -117.94, -117.94, -119.74, -118.63, -120.18, -117.5, -117.5, -118.63];

% dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_April\ghi_actual';
dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_April\ghi_actual';

dir_home = pwd;

for k = 1:length(sitenames)
    gen = sitenames{k};
    ibm_site = IBMsitenames{k};
    csvname_read  = strcat('IBM_processed_', ibm_site, '.hourly.csv');
    csvname_write = strcat('power_',         gen,      '.hourly.csv');
    
    cd(dir_work);
    M = csvread(csvname_read, 1, 0);
    cd(dir_home);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Fix 1: Time shift
%     toff=1;
%     Mp = [M(1:size(M, 1)-toff, 1: 5), M(toff+1: end, 6:end)];
%     M = Mp;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    ac_actual = ghi_to_ac_power(gen, M(:, 1), M(:, 2), M(:, 3), M(:, 4), M(:, 5), zeros(size(M, 1), 1), M(:, 6));
        
    % Calculate clear-sky GHI
    utc_year   = M(:, 1);
    utc_month  = M(:, 2);
    utc_day    = M(:, 3);
    utc_hour   = M(:, 4);
    utc_minute = M(:, 5);
    utc_second = zeros(size(M, 1), 1);    
    Time.UTCOffset(1:size(M,1),1) = zeros(size(M,1), 1); % Because we use UTC time, so utc offset is zero
    Time.year(1:size(M,1),1)   = utc_year;
    Time.month(1:size(M,1),1)  = utc_month;
    Time.day(1:size(M,1),1)    = utc_day;
    Time.hour(1:size(M,1),1)   = utc_hour;
    Time.minute(1:size(M,1),1) = utc_minute;
    Time.second(1:size(M,1),1) = utc_second;
    Location = pvl_makelocationstruct(SiteLatitude(k),SiteLongitude(k)); %Altitude is optional
    dayofyear = pvl_date2doy(Time.year, Time.month, Time.day);
    [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
    ghi_clearsky = pvl_clearsky_haurwitz(90-AppSunEl); % Clear-sky GHI

    % Plot them out
    tarray = datetime(Time.year, Time.month, Time.day, Time.hour, Time.minute, Time.second, 'TimeZone', 'UTC');
    tarray.TimeZone = 'America/Los_Angeles';
    tarray_local = datetime(tarray.Year, tarray.Month, tarray.Day, tarray.Hour, tarray.Minute, tarray.Second);

    figure();
    subplot(2, 1, 1);
    h1 = plot(tarray_local, M(:, 6), '-s');
%     h1 = plot(1: size(M, 1), ac_actual(:, end));
    set(h1, {'color'}, {'k'});
    hold on;
    h2 = plot(tarray_local, ghi_clearsky,  '-s');
    title(strcat(gen, ',', ibm_site));
    legend([h1, h2], 'Actual' ,'Clear-sky');
    ylabel('W/m^2');
    
    subplot(2, 1, 2);
    kcs_actual = M(:, 6)./ghi_clearsky;
    kcs_actual(ghi_clearsky==0) = nan;
    hist(kcs_actual, 50);
    
%     cd(dir_work);
%     cHeader = {'Year' 'Month' 'Day' 'Hour' 'Actual'}; %dummy header
%     commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
%     commaHeader = commaHeader(:)';
%     textHeader = cell2mat(commaHeader); %cHeader in text with commas
%     %write header to file
%     fid = fopen(csvname_write,'w'); 
%     fprintf(fid,'%s\n',textHeader);
%     fclose(fid);
%     dlmwrite(csvname_write,[M(:, 1:4), ac_actual],'-append');
%     cd(dir_home);

    
end
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

function explor_clearsky_index()

sitenames     = {'gen55', 'gen56',  'gen57', 'gen58', 'gen59', 'gen60',    'gen61', 'gen62', 'gen63', 'gen64'};
IBMsitenames  = {'MNCC1', 'STFC1',  'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};
SiteLatitude  = [34.31,     34.12,    34.12,   34.12,   37.41,   35.53,      35.38,   34.31,   34.31,  35.53];
SiteLongitude = [-117.5,   -117.94, -117.94, -117.94, -119.74, -118.63,    -120.18,  -117.5,  -117.5, -118.63];

for k = 1:length(sitenames)
    
    gen = sitenames{k};
    ibm_site = IBMsitenames{k};

    % dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_April\ghi_actual';
    % csvname_read  = strcat('IBM_processed_', ibm_site, '.hourly.csv');
    % cd(dir_work);
    % M = csvread(csvname_read, 1, 0);
    % cd(dir_home);
    % 
    % ghi = M(:, end);
    % 
    % Location = pvl_makelocationstruct(SiteLatitude(k),SiteLongitude(k)); %Altitude is optional
    % Time.UTCOffset(1:size(M,1),1) = zeros(size(M,1), 1); % Because we use UTC time, so utc offset is zero
    % Time.year(1:size(M,1),1)   = M(:, 1);
    % Time.month(1:size(M,1),1)  = M(:, 2);
    % Time.day(1:size(M,1),1)    = M(:, 3);
    % Time.hour(1:size(M,1),1)   = M(:, 4);
    % Time.minute(1:size(M,1),1) = zeros(size(M, 1), 1);
    % Time.second(1:size(M,1),1) = zeros(size(M, 1), 1);

    dir_work = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_April\ghi_frcst';
    csvname_read  = strcat('IBM_processed_', ibm_site, '.csv');
    cd(dir_work);
    M = csvread(csvname_read, 1, 0);
    cd(dir_home);

    ghi_05 = M(:, 6);
    ghi_50 = M(:, 7);
    ghi_95 = M(:, end);

    Location = pvl_makelocationstruct(SiteLatitude(k),SiteLongitude(k)); %Altitude is optional
    Time.UTCOffset(1:size(M,1),1) = zeros(size(M,1), 1); % Because we use UTC time, so utc offset is zero
    Time.year(1:size(M,1),1)   = M(:, 1);
    Time.month(1:size(M,1),1)  = M(:, 2);
    Time.day(1:size(M,1),1)    = M(:, 3);
    Time.hour(1:size(M,1),1)   = M(:, 4);
    Time.minute(1:size(M,1),1) = M(:, 5);
    Time.second(1:size(M,1),1) = zeros(size(M, 1), 1);


    [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location);
    ApparentZenith = 90-ApparentSunEl;

    % Clear-sky GHI: Model 1
    ghi_clearsky = pvl_clearsky_haurwitz(ApparentZenith);

    % Clear-sky GHI: Model 2
%     Location.altitude = 20; % Altitude is a must here, randomly selected
%     [ghi_clearsky, dni_clearsky, dhi_clearsky] = pvl_clearsky_ineichen(Time, Location);

    tarray = datetime(Time.year, Time.month, Time.day, Time.hour, Time.minute, Time.second, 'TimeZone', 'UTC');
    tarray.TimeZone = 'America/Los_Angeles';
    tarray_local = datetime(tarray.Year, tarray.Month, tarray.Day, tarray.Hour, tarray.Minute, tarray.Second);
    
%     figure();
%     x2 = [tarray_local', fliplr(tarray_local')];
%     inBetween = [ghi_05', fliplr(ghi_95')];
%     fill(x2, inBetween, 'k', 'FaceAlpha', 0.2);
%     hold on;
%     h1 = plot(tarray_local, ghi_clearsky, 'Color', 'r', 'LineWidth', 1.5);
%     title(strcat(gen, ',', ibm_site));
%     ylabel('GHI (W/m^-2)');
%     yyaxis right;
%     h2 = plot(tarray_local, clear_sky_index, 'b');
%     legend([h1, h2], {'Clear-sky GHI', 'Clear-sky index (CI 95%)'});
%     ylabel('Clear-sky index');
    
    figure();
    hist(clear_sky_index, 100);
    title(strcat(gen, ',', ibm_site));
    xlabel('clear sky index');
    ylabel('Frequency');
    fprintf('gen %s - site %s, max csi: %f\n', gen, ibm_site, max(clear_sky_index));

%     ac_power = ghi_to_ac_power(gen, Time.year, Time.month, Time.day, Time.hour, Time.minute, Time.second, ghi);
%     figure();
%     scatter(ghi, ac_power(:, end));
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
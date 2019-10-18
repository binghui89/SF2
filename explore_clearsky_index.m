
function explore_clearsky_index(dir_ghi_frcst)
dir_home = pwd;
% sitenames     = {'gen55', 'gen56',  'gen57', 'gen58', 'gen59', 'gen60',    'gen61', 'gen62', 'gen63', 'gen64'};
% IBMsitenames  = {'MNCC1', 'STFC1',  'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};
% SiteLatitude  = [34.31,     34.12,    34.12,   34.12,   37.41,   35.53,      35.38,   34.31,   34.31,  35.53];
% SiteLongitude = [-117.5,   -117.94, -117.94, -117.94, -119.74, -118.63,    -120.18,  -117.5,  -117.5, -118.63];

IBMsitenames  = {"CA_Topaz", "RSAC1", "RLKC1", "SBVC1", "KNNC1", "MIAC1", "MNCC1", "STFC1", "DEMC1", "COWC1"};
SiteLatitude  = [     35.38,   38.47,   40.25,   34.45,   40.71,   37.41,   34.31,   34.12,   35.53,   39.12];
SiteLongitude = [   -120.18, -122.71, -123.31, -119.70, -123.92, -119.74, -117.50, -117.94, -118.63, -123.07];


for k = 1:length(IBMsitenames)
    
    ibm_site = IBMsitenames{k};

    dir_work = dir_ghi_frcst;
    csvname_read  = strcat('IBM_processed_', ibm_site, '.csv');
    cd(dir_work);
    M = csvread(csvname_read, 1, 0);
    table_ghi = readtable(csvname_read);
    cd(dir_home);

    Location = pvl_makelocationstruct(SiteLatitude(k),SiteLongitude(k)); %Altitude is optional
    Time.UTCOffset(1:size(table_ghi,1),1) = zeros(size(table_ghi,1), 1); % Because we use UTC time, so utc offset is zero
    Time.year(1:size(table_ghi,1),1)   = table_ghi.Year;
    Time.month(1:size(table_ghi,1),1)  = table_ghi.Month;
    Time.day(1:size(table_ghi,1),1)    = table_ghi.Day;
    Time.hour(1:size(table_ghi,1),1)   = table_ghi.Hour;
    Time.minute(1:size(table_ghi,1),1) = table_ghi.Minute;
    Time.second(1:size(table_ghi,1),1) = zeros(size(table_ghi, 1), 1);


    [SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(Time, Location);
    ApparentZenith = 90-ApparentSunEl;

    ghi_clearsky = clear_sky_ghi(Time, Location);

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
    
%     figure();
%     hist(clear_sky_index, 100);
%     title(strcat(gen, ',', ibm_site));
%     xlabel('clear sky index');
%     ylabel('Frequency');
%     fprintf('gen %s - site %s, max csi: %f\n', gen, ibm_site, max(clear_sky_index));

%     ac_power = ghi_to_ac_power(gen, Time.year, Time.month, Time.day, Time.hour, Time.minute, Time.second, ghi);
%     figure();
%     scatter(ghi, ac_power(:, end));

    figure();
    ax1 = subplot(2, 1, 1);
    array_csi = table_ghi{:, 6:end}./repmat(ghi_clearsky, 1, size(table_ghi, 2)-5);
    plot(tarray_local, array_csi(:, 1:end-1), 'k', tarray_local, array_csi(:, end), 'b');
    ylim([0, 1.5]);
    ylabel('Clear-sky index');
    
    ax2 = subplot(2, 1, 2);
    plot(tarray_local, table_ghi{:, 6:end-1}, 'k', tarray_local, table_ghi{:, end}, 'b', tarray_local, ghi_clearsky, 'r');
    linkaxes([ax1, ax2],'x');
    ylabel('GHI (W/m^2)');
    
    suptitle(ibm_site);
    
end
end

function ghi_clearsky = clear_sky_ghi(time, location)
% Time is a struct with the following elements, which can be column vectors all of the same length.
% 
% Time.year = The year in the gregorian calendar.
% Time.month = the month of the year (January = 1 to December = 12).
% Time.day = the day of the month.
% Time.hour = the hour of the day.
% Time.minute = the minute of the hour.
% Time.second = the second of the minute.
% Time.UTCOffset = the UTC offset code, using the convention that a positive UTC offset is for time zones east of the prime meridian (e.g. EST = -5).
% Location is a struct with the following elements, which can be column vectors all of the same length.
% 
% Location.latitude = vector or scalar latitude in decimal degrees (positive is northern hemisphere).
% Location.longitude = vector or scalar longitude in decimal degrees (positive is east of prime meridian).
% Location.altitude = an optional component of the Location struct, used to calculate atmospheric pressure (see pvl_alt2pres).

[SunAz, SunEl, ApparentSunEl, SolarTime]=pvl_ephemeris(time, location);
ApparentZenith = 90-ApparentSunEl;

% Clear-sky GHI: Model 1
ghi_clearsky = pvl_clearsky_haurwitz(ApparentZenith);

% Clear-sky GHI: Model 2
% Location.altitude = 20; % Altitude is a must here, randomly selected
% [ghi_clearsky, dni_clearsky, dhi_clearsky] = pvl_clearsky_ineichen(Time, Location);

end

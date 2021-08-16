%% Load Carlo's IBM data, forecast 15-min
add_pvlib;

uniquegen     = {'gen55', 'gen56', 'gen59', 'gen60', 'gen61'};
uniquegensite = {'MNCC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz'};
uniquegenlat  = [34.31, 34.12, 37.41, 35.53, 35.38]; 
uniquegenlon  = [-117.50,-117.94,-119.74,-118.63,-120.18];

allgen  = {'gen55', 'gen56', 'gen57', 'gen58', 'gen59', 'gen60', 'gen61',    'gen62', 'gen63', 'gen64'};
allgensite = {'MNCC1', 'STFC1', 'STFC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz', 'MNCC1', 'MNCC1', 'DEMC1'};
allgenlat = [34.31,34.12,34.12,34.12,37.41,35.53,35.38,34.31,34.31,35.53];
allgenlon = [-117.5,-117.94, -117.94, -117.94,-119.74, -118.63, -120.18,-117.5,-117.5,-118.63];

dirhome = pwd;
cell_pwr = cell(numel(uniquegen), 1);
dirwork = 'C:\Users\bxl180002\Downloads\predictions\Processed_Predictions\F_pwr.15min';
for i = 1: numel(uniquegen)
    gen = uniquegen{i};
%     for m = 1: numel(month_ibm)
%         mstr = num2str(month_ibm(m));
%         dirwork = strcat('C:\Users\bxl180002\OneDrive\Tmp_RampSolar\Code\IBM\F_pwr.', mstr, '.15min');
%         csvname = strcat(dirwork, '\', 'frcst_', gen, '.csv');
%         cd(dirwork);
%         if m == 1
%             T_pwr = readtable(csvname, 'Delimiter', ',');
%         else
%             T_pwr = vertcat(T_pwr, readtable(csvname, 'Delimiter', ','));
%         end
%     end
    csvname = strcat(dirwork, '\', 'frcst_', gen, '.csv');
    T_pwr = readtable(csvname, 'Delimiter', ',');
    [~, ia, ~] = unique(T_pwr{:, {'Year', 'Month', 'Day', 'Hour', 'Minute'}}, 'rows');
    T_pwr = T_pwr(ia, :);
    T_pwr.TIME = datetime(T_pwr.Year, T_pwr.Month, T_pwr.Day, T_pwr.Hour, T_pwr.Minute, 0, 'TimeZone', 'UTC');
    
    % Calculate solar elevation to filter invalid clear-sky index
    Location = pvl_makelocationstruct(uniquegenlat(i), uniquegenlon(i));
    Time.year   = T_pwr.TIME.Year;
    Time.month  = T_pwr.TIME.Month;
    Time.day    = T_pwr.TIME.Day;
    Time.hour   = T_pwr.TIME.Hour;
    Time.minute = T_pwr.TIME.Minute;
    Time.second = T_pwr.TIME.Second;
    Time.UTCOffset = zeros(size(T_pwr, 1), 1);
    [SunAz, SunEl, ApparentSunEl, SolarTime] = pvl_ephemeris(Time, Location);
    T_pwr.ApparentSunEl = ApparentSunEl;
    T_pwr.SunEl = SunEl;
    T_pwr{T_pwr.ApparentSunEl<=3, 'ghi_cs'} = 0; % Consider only solar elevation > 3 degree
    T_pwr{T_pwr.ApparentSunEl<=3, 'pwr_cs'} = 0; % Consider only solar elevation > 3 degree
    
    % Calculate uncertainty bandwidth and variability of clear-sky index (k)
    T_pwr.k_p025 = T_pwr.ghi_p025./T_pwr.ghi_cs;
    T_pwr.k_p050 = T_pwr.ghi_p050./T_pwr.ghi_cs;
    T_pwr.k_p075 = T_pwr.ghi_p075./T_pwr.ghi_cs;
    T_pwr.k_width = T_pwr.k_p075 - T_pwr.k_p025;
    
    % Calculate uncertainty bandwidth and variability of clear-sky index of PV (k_PV)
    T_pwr.kpv_p025 = T_pwr.pwr_p025./T_pwr.pwr_cs;
    T_pwr.kpv_p050 = T_pwr.pwr_p050./T_pwr.pwr_cs;
    T_pwr.kpv_p075 = T_pwr.pwr_p075./T_pwr.pwr_cs;
    T_pwr.kpv_width = T_pwr.kpv_p075 - T_pwr.kpv_p025;
    
%     figure();
%     subplot(2, 1, 1);
%     hist(T_pwr.k_p050);
%     title('k');
%     subplot(2, 1, 2);
%     hist(T_pwr.kpv_p050);
%     title('kpv');
    
    cell_pwr{i} = T_pwr;
end
cd(dirhome);


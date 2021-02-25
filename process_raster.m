% CA poly example
fname = 'Solar forecast (PAIRS Team) (select sites)-GHI Quantile Forecast-horizon_60_quantile_50[solar_forecast_median]-02_15_2020T20_00_00.tiff';
[A1, R1] = readgeoraster(fname); % Example
info = georasterinfo(fname);
m = info.MissingDataIndicator;
A = standardizeMissing(A1,m);
figure();
imagesc(R1.LongitudeLimits(:), flipud(R1.LatitudeLimits(:)), A); % Flip because MATLAB array y axis points downward
set(gca,'YDir', 'normal'); % Set origin to the bottom left corner
colormap jet;

% Square example
fname = 'Solar forecast (PAIRS Team) (select sites)-GHI Quantile Forecast-horizon_60_quantile_50[solar_forecast_median]-02_01_2020T20_00_00.tiff';
[A1, R1] = readgeoraster(fname); % Example
info = georasterinfo(fname);
m = info.MissingDataIndicator;
A = standardizeMissing(A1,m);
figure();
imagesc(R1.LongitudeLimits(:), flipud(R1.LatitudeLimits(:)), A);
set(gca,'YDir', 'normal'); % Set origin to the bottom left corner
colormap jet;

%% Check file availability and data validity, run on Ganymede for better performance
dirhome = pwd;
if ispc
    dirwork = 'C:\Users\bxl180002\Downloads\RampSolar\IBM\CA_square_daily_202002\content';
elseif isunix
    dirwork = '/home/bxl180002/scratch/SF2/IBM/CA_square_daily_202002/content/'; % Path on Ganymede
end

% datetime_start = datetime(2020, 2, 1, 9, 0, 0, 'TimeZone', 'UTC');
% datetime_end = datetime(2020, 2, 2, 2, 50, 0, 'TimeZone', 'UTC');
deltat = duration(0, 10, 0);
ar_quantiles = [5, 25, 50, 75, 95];
ar_datetime = [];

for d = 1:27
    datetime_start = datetime(2020, 2, d, 9, 0, 0, 'TimeZone', 'UTC');
    datetime_end = datetime(2020, 2, d+1, 2, 50, 0, 'TimeZone', 'UTC');
    time_per_day = [datetime_start: deltat: datetime_end]';
    ar_datetime = [ar_datetime; time_per_day];
end

cell_data = cell(5, 1);
for j = 1: numel(ar_quantiles)
    cell_data{j} = nan(320, 480, numel(ar_datetime)); % 320 x 480 is my data size, hopefully all figures shoud have the same size
end

% Check file availability and read data
ar_istiff = zeros(numel(ar_datetime), numel(ar_quantiles));
tic;
cd(dirwork);
for i = 1: numel(ar_datetime)
    for j = 1: numel(ar_quantiles)
        tiff_name = strcat(... 
        'Solar forecast (PAIRS Team) (select sites)-GHI Quantile Forecast-horizon_60_quantile_', ...
        sprintf('%d', ar_quantiles(j)), ...
        '[solar_forecast_median]-', ...
        sprintf('%02d_%02d_%4dT%02d_%02d_%02d', ar_datetime(i).Month, ar_datetime(i).Day, ar_datetime(i).Year, ar_datetime(i).Hour, ar_datetime(i).Minute, ar_datetime(i).Second), ...
        '.tiff');
        if isfile(tiff_name)
            ar_istiff(i, j) = 1;
            [A, R] = readgeoraster(tiff_name); % Example
            info = georasterinfo(tiff_name);
            m = info.MissingDataIndicator;
            A = standardizeMissing(A,m);
            cell_data{j}(:, :, i) = A;
        else
            ar_istiff(i, j) = 0;
            fprintf('No file found: %s\n', tiff_name);
        end
    end
end
cd(dirhome);
toc;
disp('Data loaded!');

ar_datetime_local = ar_datetime;
ar_datetime_local.TimeZone = 'America/Los_Angeles';
T_istiff = [array2table(ar_datetime, 'VariableNames', {'Time'}), array2table(ar_istiff, 'VariableNames', {'q5', 'q25', 'q50', 'q75', 'q95'})];
T_istiff_local = [array2table(ar_datetime_local, 'VariableNames', {'Time'}), array2table(ar_istiff, 'VariableNames', {'q5', 'q25', 'q50', 'q75', 'q95'})];

% Check quantile crossing
tic;
ar_error = nan(numel(ar_datetime), numel(ar_quantiles)-1);
cd(dirwork);
for i = 1: numel(ar_datetime)
    for j = 1: numel(ar_quantiles)-1
        tiff_name_1 = strcat(... 
        'Solar forecast (PAIRS Team) (select sites)-GHI Quantile Forecast-horizon_60_quantile_', ...
        sprintf('%d', ar_quantiles(j)), ...
        '[solar_forecast_median]-', ...
        sprintf('%02d_%02d_%4dT%02d_%02d_%02d', ar_datetime(i).Month, ar_datetime(i).Day, ar_datetime(i).Year, ar_datetime(i).Hour, ar_datetime(i).Minute, ar_datetime(i).Second), ...
        '.tiff');
        tiff_name_2 = strcat(... 
        'Solar forecast (PAIRS Team) (select sites)-GHI Quantile Forecast-horizon_60_quantile_', ...
        sprintf('%d', ar_quantiles(j+1)), ...
        '[solar_forecast_median]-', ...
        sprintf('%02d_%02d_%4dT%02d_%02d_%02d', ar_datetime(i).Month, ar_datetime(i).Day, ar_datetime(i).Year, ar_datetime(i).Hour, ar_datetime(i).Minute, ar_datetime(i).Second), ...
        '.tiff');
        if isfile(tiff_name_1)&&isfile(tiff_name_2)
            [A1, R1] = readgeoraster(tiff_name_1); % Example
            [A2, R2] = readgeoraster(tiff_name_2);
            info1 = georasterinfo(tiff_name_1);
            m1 = info1.MissingDataIndicator;
            A1 = standardizeMissing(A1,m1);
            info2 = georasterinfo(tiff_name_2);
            m2 = info2.MissingDataIndicator;
            A2 = standardizeMissing(A2,m2);
            ar_error(i, j) = sum(A2(:) - A1(:)<0, 'omitnan');
        end
    end
end
cd(dirhome);
toc;

T_error = [array2table(ar_datetime, 'VariableNames', {'Time'}), array2table(ar_error, 'VariableNames', {'q25<q5', 'q50<q25', 'q75<q50', 'q95<q75'})];
T_error_local = [array2table(ar_datetime_local, 'VariableNames', {'Time'}), array2table(ar_error, 'VariableNames', {'q25<q5', 'q50<q25', 'q75<q50', 'q95<q75'})];

if false
    T_istiff_local.Properties.VariableNames{'Time'} = 'LocalTime';
    writetable([T_istiff(:, 'Time'), T_istiff_local], 'summary.xlsx', 'Sheet', 1);
    T_error_local.Properties.VariableNames{'Time'} = 'LocalTime';
    writetable([T_error(:, 'Time'), T_error_local], 'summary.xlsx', 'Sheet', 2, 'WriteMode','Append');
end

%% Step 1: Read all PV plants coordinates from EIA-860 2019 table
fname = 'Solar forecast (PAIRS Team) (select sites)-GHI Quantile Forecast-horizon_60_quantile_50[solar_forecast_median]-02_01_2020T20_00_00.tiff';
[A1, R1] = readgeoraster(fname); % Example, we need the R1 here to give estimate of lat and lon ranges
info = georasterinfo(fname);
m = info.MissingDataIndicator;
A = standardizeMissing(A1,m);

T_eia860_caiso = readtable('eia860_2019_caiso.csv', 'Delimiter',',');
cap = T_eia860_caiso.NameplateCapacity_MW_;

% Group PV generator capacity by plant code
T_pv = grpstats(T_eia860_caiso, 'PlantCode', {'sum', 'mean'}, 'DataVars', {'NameplateCapacity_MW_', 'Latitude', 'Longitude'});
grpstats(T_eia860_caiso, 'PlantCode', {'sum', 'mean'}, 'DataVars', {'NameplateCapacity_MW_', 'Latitude', 'Longitude'});
T_pv.Properties.VariableNames{'sum_NameplateCapacity_MW_'} = 'TotalCapacity';
T_pv.Properties.VariableNames{'mean_Latitude'} = 'Latitude';
T_pv.Properties.VariableNames{'mean_Longitude'} = 'Longitude';

lat_edge = R1.LatitudeLimits(1): R1.CellExtentInLatitude: R1.LatitudeLimits(2);
lon_edge = R1.LongitudeLimits(1): R1.CellExtentInLongitude: R1.LongitudeLimits(2);
lat_edge = flipud(lat_edge(:)); % Flip because MATLAB array points downward while Cartecian coordinate points upward
lon_edge = lon_edge(:);
lat_center = (lat_edge(1:end-1) + lat_edge(2:end))/2;
lon_center = (lon_edge(1:end-1) + lon_edge(2:end))/2;

% The boundary of available IBM raster, i.e., not nan
BOUNDARY = struct('WEST', nan, 'EAST', nan, 'SOUTH', nan, 'NORTH', nan); 
BOUNDARY.NORTH = lat_center(find(any(~isnan(A), 2), 1, 'first'));
BOUNDARY.SOUTH = lat_center(find(any(~isnan(A), 2), 1, 'last')); 
BOUNDARY.WEST  = lon_center(find(any(~isnan(A), 1), 1, 'first'));
BOUNDARY.EAST  = lon_center(find(any(~isnan(A), 1), 1, 'last'));

% Group PV plant capacities by cell locations in IBM raster forecast
T_pv{:, 'nx'} = nan;
T_pv{:, 'ny'} = nan;
for i = 1: height(T_pv)
    ar_nx = find(abs(lon_center - T_pv{i, 'Longitude'}) < R1.CellExtentInLongitude/2);
    if numel(ar_nx) == 1
        T_pv{i, 'nx'} = ar_nx;
    else
        disp('WARNING: CANNOT LOCATE THIS PV PLANT: LONGITUDE!');
    end
    ar_ny = find(abs(lat_center(:) - T_pv{i, 'Latitude'}) < R1.CellExtentInLatitude/2); 
    if numel(ar_ny) == 1
        T_pv{i, 'ny'} = ar_ny;
    else
        disp('WARNING: CANNOT LOCATE THIS PV PLANT: LATITUDE!');
    end
end
T_pvcell = grpstats(T_pv, {'nx', 'ny'}, {'sum'}, 'DataVars', {'TotalCapacity'}); % nx is the index of latitude of the PV plant, ny is the index of longitude of the PV plant
T_pvcell.Properties.VariableNames{'sum_TotalCapacity'} = 'TotalCapacity';
% T_pvcell{:, 'ind'} = sub2ind([size(cell_data{1}, 1), size(cell_data{1}, 2)], T_pvcell{:, 'ny'}, T_pvcell{:, 'nx'}); % Sub to index
T_pvcell{:, 'ind'} = sub2ind([size(A, 1), size(A, 2)], T_pvcell{:, 'ny'}, T_pvcell{:, 'nx'}); % Sub to index
T_pvcell{:, 'cell_lat'} = lat_center(T_pvcell{:, 'ny'}); % Note this is the center coordinate of the IBM cell, not the actual PV plant
T_pvcell{:, 'cell_lon'} = lon_center(T_pvcell{:, 'nx'}); % Note this is the center coordinate of the IBM cell, not the actual PV plant
T_pvcell{:, 'within_boundary'} = (T_pvcell{:, 'cell_lon'}>=BOUNDARY.WEST) & (T_pvcell{:, 'cell_lon'}<=BOUNDARY.EAST) & (T_pvcell{:, 'cell_lat'}>=BOUNDARY.SOUTH) & (T_pvcell{:, 'cell_lat'}<=BOUNDARY.NORTH);

% cell_pvcellghi = cell(5, 1);
% for j = 1: 5
%     tmp = reshape(cell_data{j}, size(cell_data{j}, 1)*size(cell_data{j}, 2), size(cell_data{j}, 3));
%     cell_pvcellghi{j} = tmp(T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'ind'}, :)'; % First dimension: Time, second dimension: Number of cells that include PV plants
% end

%% Step 2: For small memory machine, load data and process it based on CAISO PV plant locations
cell_dircontent = {'C:\Users\bxl180002\Downloads\RampSolar\IBM\CA_square_daily_202002\content', 'C:\Users\bxl180002\Downloads\RampSolar\IBM\CA_square_daily_202003\content'};
cell_datetime = cell(2, 1);
deltat = duration(0, 10, 0);

% All timestamps in Feb 2020
ar_datetime = [];
for d = 1:29
    datetime_start = datetime(2020, 2, d, 9, 0, 0, 'TimeZone', 'UTC');
    datetime_end = datetime(2020, 2, d, 2, 50, 0, 'TimeZone', 'UTC') + duration(24, 0, 0);
    time_per_day = [datetime_start: deltat: datetime_end]';
    ar_datetime = [ar_datetime; time_per_day];
end
cell_datetime{1} = ar_datetime;

% All time stamps in March 2020
ar_datetime = [];
for d = 1:31
    datetime_start = datetime(2020, 3, d, 9, 0, 0, 'TimeZone', 'UTC');
    datetime_end = datetime(2020, 3, d, 2, 50, 0, 'TimeZone', 'UTC') + duration(24, 0, 0);
    time_per_day = [datetime_start: deltat: datetime_end]';
    ar_datetime = [ar_datetime; time_per_day];
end
cell_datetime{2} = ar_datetime;

% Prepare file paths from all months
ar_quantiles = [5, 25, 50, 75, 95];
cell_filepath = {}; % First dim: Time, Second dim: Quantile
tic;
for m = 1:numel(cell_datetime)
    ar_datetime = cell_datetime{m};
    dircontent = cell_dircontent{m};
    for i = 1: size(ar_datetime, 1)
        cellrow = {};
        for j = 1: numel(ar_quantiles)
            tiff_name = strcat(... 
                'Solar forecast (PAIRS Team) (select sites)-GHI Quantile Forecast-horizon_60_quantile_', ...
                sprintf('%d', ar_quantiles(j)), ...
                '[solar_forecast_median]-', ...
                sprintf('%02d_%02d_%4dT%02d_%02d_%02d', ar_datetime(i).Month, ar_datetime(i).Day, ar_datetime(i).Year, ar_datetime(i).Hour, ar_datetime(i).Minute, ar_datetime(i).Second), ...
                '.tiff');
            cellrow = [cellrow fullfile(dircontent, tiff_name)];
        end
        cell_filepath = [cell_filepath; cellrow];
    end
end
toc;
ar_datetime = [cell_datetime{1}; cell_datetime{2}];
ar_datetime_local.TimeZone = 'America/Los_Angeles';
disp('File names ready!');

% Load data and only save those data from locations with PV plants
cell_pvcellghi = cell(numel(ar_quantiles), 1);
for j = 1: size(cell_filepath, 2)
    tic;
    data_quantile = nan(320, 480, numel(ar_datetime));
    for i = 1: size(cell_filepath, 1)
        tiff_name = string(cell_filepath(i, j));
        if isfile(tiff_name)
            ar_istiff(i, j) = 1;
            [A, R] = readgeoraster(tiff_name); % Example
            info = georasterinfo(tiff_name);
            m = info.MissingDataIndicator;
            A = standardizeMissing(A,m);
            data_quantile(:, :, i) = A;
        else
            ar_istiff(i, j) = 0;
            fprintf('No file found: %s\n', tiff_name);
        end
    end
    tmp = reshape(data_quantile, size(data_quantile, 1)*size(data_quantile, 2), size(data_quantile, 3));
    cell_pvcellghi{j} = tmp(T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'ind'}, :)'; % First dimension: Time, second dimension: Number of cells that include PV plants
    clear tmp;
    clear data_quantile;
    toc;
end

if false
    save('raster.mat', 'ar_datetime', 'ar_datetime_local', 'ar_quantiles', 'BOUNDARY', 'deltat', 'lat_center', 'lat_edge', 'lon_center', 'lon_edge', 'T_pv', 'T_pvcell');
end
%% Load CAISO data RTD and RTPD
load_caiso_data;
ar_datetime_hour = datetime(ar_datetime.Year, ar_datetime.Month, ar_datetime.Day, ar_datetime.Hour, 0, 0, 'TimeZone', 'UTC');
T_rtpd = T_rtpd((T_rtpd{:, 'HOUR_START'}>=min(ar_datetime))&(T_rtpd{:, 'HOUR_START'}<=max(ar_datetime)), :); % Only look at the hours with raster forecast

%% Save csv files for Cong
dir_cong = 'Cong';
if ~isfolder(dir_cong)
    mkdir(dir_cong);
end
dirhome = pwd;
cd(dir_cong);
% for i = 1: length(cell_pvcellghi)
%     T_ghi_forcong = [array2table(ar_datetime, 'VariableNames', {'TIME'}) array2table(cell_pvcellghi{i}, 'VariableNames', cellstr(num2str(T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'ind'})))];
%     writetable(T_ghi_forcong, strcat(num2str(i), '.csv'));
% end
T_rtpd_forcong = T_rtpd(:, {'HOUR_START', 'error_max', 'error_min'});
T_rtpd_forcong = grpstats(T_rtpd_forcong, 'HOUR_START', {'max', 'min'}, 'DataVars', {'error_max', 'error_min'});
T_rtpd_forcong = T_rtpd_forcong(:, {'HOUR_START', 'max_error_max', 'min_error_min'});
T_rtpd_forcong.Properties.VariableNames{'max_error_max'} = 'FRU';
T_rtpd_forcong.Properties.VariableNames{'min_error_min'} = 'FRD';
writetable(T_rtpd_forcong, 'rtpd.csv');
writetable(T_pvcell, 'T_pvcell.csv');
cd(dirhome);

%% Clear-sky GHI calculation
add_pvlib();

T_pvcell_within = T_pvcell(T_pvcell{:, 'within_boundary'}==1, :);
ar_ghics = nan(size(ar_datetime, 1), height(T_pvcell_within));
ar_AppSunEl = nan(size(ar_datetime, 1), height(T_pvcell_within));

for k = 1: height(T_pvcell_within)
    Location = pvl_makelocationstruct(T_pvcell_within{k, 'cell_lat'}, T_pvcell_within{k, 'cell_lon'});

    Time.UTCOffset = zeros(size(ar_datetime, 1), 1); % Because IBM uses UTC time, so utc offset is zero
    Time.year   = ar_datetime.Year;
    Time.month  = ar_datetime.Month;
    Time.day    = ar_datetime.Day;
    Time.hour   = ar_datetime.Hour;
    Time.minute = ar_datetime.Minute;
    Time.second = ar_datetime.Second;
    [~, ~, AppSunEl, ~] = pvl_ephemeris(Time, Location);
    ghi_cs = pvl_clearsky_haurwitz(90-AppSunEl);
    ar_ghics(:, k) = ghi_cs;
    ar_AppSunEl(:, k) = AppSunEl;
end


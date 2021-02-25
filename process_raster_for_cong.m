dirhome = pwd;
if ispc
    cell_dircontent = {'C:\Users\bxl180002\Downloads\RampSolar\IBM\CA_square_daily_202002\content', 'C:\Users\bxl180002\Downloads\RampSolar\IBM\CA_square_daily_202003\content'};
elseif isunix
    cell_dircontent = {...
        '/home/bxl180002/scratch/SF2/IBM/CA_square_daily_202002/content', ...
        '/home/bxl180002/scratch/SF2/IBM/CA_square_daily_202003/content' ...
        };
end

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
ar_datetime_local = ar_datetime;
ar_datetime_local.TimeZone = 'America/Los_Angeles';
disp('File names ready!');

% Prepare data container
cell_rawdata = cell(5, 1);
for j = 1: numel(ar_quantiles)
    cell_rawdata{j} = nan(320, 480, numel(ar_datetime)); 
end
ar_istiff = nan(numel(ar_datetime), numel(ar_quantiles));

% Load data
tic;
for i = 1: numel(ar_datetime)
    for j = 1: numel(ar_quantiles)
        tiff_name = string(cell_filepath(i, j));
        if isfile(tiff_name)
            ar_istiff(i, j) = 1;
            [A, R] = readgeoraster(tiff_name); % Example
            info = georasterinfo(tiff_name);
            m = info.MissingDataIndicator;
            A = standardizeMissing(A,m);
            cell_rawdata{j}(:, :, i) = A;
        else
            ar_istiff(i, j) = 0;
            fprintf('No file found: %s\n', tiff_name);
        end
    end
end
toc;

% Calculate clear-sky GHI
add_pvlib();

ghi_cs = nan(320, 480, numel(ar_datetime));

lat_edge = R.LatitudeLimits(1): R.CellExtentInLatitude: R.LatitudeLimits(2);
lon_edge = R.LongitudeLimits(1): R.CellExtentInLongitude: R.LongitudeLimits(2);
lat_edge = flipud(lat_edge(:)); % Flip because MATLAB array points downward while Cartecian coordinate points upward
lon_edge = lon_edge(:);
lat_center = (lat_edge(1:end-1) + lat_edge(2:end))/2;
lon_center = (lon_edge(1:end-1) + lon_edge(2:end))/2;

Time.UTCOffset = zeros(size(ar_datetime, 1), 1); % Because IBM uses UTC time, so utc offset is zero
Time.year   = ar_datetime.Year;
Time.month  = ar_datetime.Month;
Time.day    = ar_datetime.Day;
Time.hour   = ar_datetime.Hour;
Time.minute = ar_datetime.Minute;
Time.second = ar_datetime.Second;

tic;
for i = 1: length(lat_center)
    for j = 1: length(lon_center)
        Location = pvl_makelocationstruct(lat_center(i), lon_center(j));
        [~, ~, AppSunEl, ~] = pvl_ephemeris(Time, Location);
        ghi_cs(i, j, :) = pvl_clearsky_haurwitz(90-AppSunEl);
    end
end
toc;

% Calculate clear-sky index and prepare hourly classifier
cell_kcs = cell(5, 1);
for j = 1: numel(ar_quantiles)
    cell_kcs{j} = cell_rawdata{j}./ghi_cs; 
end
ar_hourstart = unique(datetime(ar_datetime.Year, ar_datetime.Month, ar_datetime.Day, ar_datetime.Hour, 0, 0, 'TimeZone', 'UTC'));

ar_kcs_mean = nan(320, 480, numel(ar_datetime)/6);
ar_kcs_std  = nan(320, 480, numel(ar_datetime)/6);
ar_kcs_vrb  = nan(320, 480, numel(ar_datetime)/6);
ar_kcs_wmean   = nan(320, 480, numel(ar_datetime)/6);

ar_kcs_p25 = cell_kcs{2};
ar_kcs_p50 = cell_kcs{3};
ar_kcs_p75 = cell_kcs{4};
for j = 1: numel(ar_hourstart)
    istart = 6*(j-1) + 1;
    iend   = 6*j;
    ar_kcs_mean(:, :, j) = mean(ar_kcs_p50(:, :, istart:iend), 3, 'omitnan');
    ar_kcs_std(:, :, j) = std(ar_kcs_p50(:, :, istart:iend), 0, 3, 'omitnan');
    ar_kcs_vrb(:, :, j) = sqrt(1/5.*sum(diff(ar_kcs_p50(:, :, istart:iend), 1, 3).^2, 3));
    ar_kcs_wmean(:, :, j) = mean(ar_kcs_p75(:, :, istart:iend) - ar_kcs_p25(:, :, istart:iend), 3, 'omitnan');
end

%% I. Write full square data
ar_hourstart_posix = posixtime(ar_hourstart);
dirhome = pwd;

cd('/home/bxl180002/scratch/SF2/');

h5create('full_square.h5','/kcs_mean',size(ar_kcs_mean));
h5write('full_square.h5', '/kcs_mean', ar_kcs_mean);

h5create('full_square.h5','/kcs_std',size(ar_kcs_std));
h5write('full_square.h5', '/kcs_std', ar_kcs_std);

h5create('full_square.h5','/kcs_vrb',size(ar_kcs_vrb));
h5write('full_square.h5', '/kcs_vrb', ar_kcs_vrb);

h5create('full_square.h5','/kcs_wmean',size(ar_kcs_wmean));
h5write('full_square.h5', '/kcs_wmean', ar_kcs_wmean);

h5create('full_square.h5','/timeutc_posix',size(ar_hourstart_posix));
h5write('full_square.h5', '/timeutc_posix', ar_hourstart_posix);

h5create('full_square.h5','/lat_center',size(lat_center));
h5write('full_square.h5', '/lat_center', lat_center);

h5create('full_square.h5','/lon_center',size(lon_center));
h5write('full_square.h5', '/lon_center', lon_center);

cd(dirhome);

%% CA poly, find out the coordinates of the boundary of the the minimum square 
fname = 'Solar forecast (PAIRS Team) (select sites)-GHI Quantile Forecast-horizon_60_quantile_50[solar_forecast_median]-02_15_2020T20_00_00.tiff';
[A1, R1] = readgeoraster(fname); % Example
info = georasterinfo(fname);
m = info.MissingDataIndicator;
A1 = standardizeMissing(A1,m);
figure();
imagesc(R1.LongitudeLimits(:), flipud(R1.LatitudeLimits(:)), A1); % Flip because MATLAB array y axis points downward
set(gca,'YDir', 'normal'); % Set origin to the bottom left corner
colormap jet;

lat_edge_1 = R1.LatitudeLimits(1): R1.CellExtentInLatitude: R1.LatitudeLimits(2);
lon_edge_1 = R1.LongitudeLimits(1): R1.CellExtentInLongitude: R1.LongitudeLimits(2);
lat_edge_1 = flipud(lat_edge_1(:)); % Flip because MATLAB array points downward while Cartecian coordinate points upward
lon_edge_1 = lon_edge_1(:);
lat_center_1 = (lat_edge_1(1:end-1) + lat_edge_1(2:end))/2;
lon_center_1 = (lon_edge_1(1:end-1) + lon_edge_1(2:end))/2;

% The boundary of the smallest square that contains CA 
BOUNDARY_CA = struct('WEST', nan, 'EAST', nan, 'SOUTH', nan, 'NORTH', nan); 
BOUNDARY_CA.NORTH = lat_center_1(find(any(~isnan(A1), 2), 1, 'first'));
BOUNDARY_CA.SOUTH = lat_center_1(find(any(~isnan(A1), 2), 1, 'last')); 
BOUNDARY_CA.WEST  = lon_center_1(find(any(~isnan(A1), 1), 1, 'first'));
BOUNDARY_CA.EAST  = lon_center_1(find(any(~isnan(A1), 1), 1, 'last'));

% Show the cropped CA picture
figure();
imagesc(...
    lon_center_1((lon_center_1<=BOUNDARY_CA.EAST)&(lon_center_1>=BOUNDARY_CA.WEST)), ...
    lat_center_1((lat_center_1<=BOUNDARY_CA.NORTH)&(lat_center_1>=BOUNDARY_CA.SOUTH)), ...
    A1((lat_center_1<=BOUNDARY_CA.NORTH)&(lat_center_1>=BOUNDARY_CA.SOUTH), (lon_center_1<=BOUNDARY_CA.EAST)&(lon_center_1>=BOUNDARY_CA.WEST)) ...
    ); % Flip because MATLAB array y axis points downward
set(gca,'YDir', 'normal'); % Set origin to the bottom left corner
colormap jet;

% Now, apply the boundary to the full square
cropped_lat = (lat_center<=BOUNDARY_CA.NORTH)&(lat_center>=BOUNDARY_CA.SOUTH);
cropped_lon = (lon_center<=BOUNDARY_CA.EAST)&(lon_center>=BOUNDARY_CA.WEST);

% Show the cropped full picture for visual inspection
figure();
A_tmp = squeeze(cell_rawdata{3}(:, :, 67)); % 67 is 12pm PST 2020-02-01
imagesc(lon_center(cropped_lon), lat_center(cropped_lat), A_tmp(cropped_lat, cropped_lon, :)); % Flip because MATLAB array y axis points downward
set(gca,'YDir', 'normal'); % Set origin to the bottom left corner
colormap jet;

%% II. Now, write the cropped square
cd('/home/bxl180002/scratch/SF2/');

h5create('cropped_square.h5','/kcs_mean',size(ar_kcs_mean(cropped_lat, cropped_lon, :)));
h5write('cropped_square.h5', '/kcs_mean', ar_kcs_mean(cropped_lat, cropped_lon, :));

h5create('cropped_square.h5','/kcs_std',size(ar_kcs_std(cropped_lat, cropped_lon, :)));
h5write('cropped_square.h5', '/kcs_std', ar_kcs_std(cropped_lat, cropped_lon, :));

h5create('cropped_square.h5','/kcs_vrb',size(ar_kcs_vrb(cropped_lat, cropped_lon, :)));
h5write('cropped_square.h5', '/kcs_vrb', ar_kcs_vrb(cropped_lat, cropped_lon, :));

h5create('cropped_square.h5','/kcs_wmean',size(ar_kcs_wmean(cropped_lat, cropped_lon, :)));
h5write('cropped_square.h5', '/kcs_wmean', ar_kcs_wmean(cropped_lat, cropped_lon, :));

h5create('cropped_square.h5','/timeutc_posix',size(ar_hourstart_posix));
h5write('cropped_square.h5', '/timeutc_posix', ar_hourstart_posix);

h5create('cropped_square.h5','/lat_center',size(lat_center(cropped_lat)));
h5write('cropped_square.h5', '/lat_center', lat_center(cropped_lat));

h5create('cropped_square.h5','/lon_center',size(lon_center(cropped_lon)));
h5write('cropped_square.h5', '/lon_center', lon_center(cropped_lon));

cd(dirhome);

%% Find the coordinates of the cells with PV plants in CAISO

% Load the coordinatees of all PV plants in CAISO
T_eia860_caiso = readtable('eia860_2019_caiso.csv', 'Delimiter',',');
cap = T_eia860_caiso.NameplateCapacity_MW_;

% Group PV generator capacity by plant code
T_pv = grpstats(T_eia860_caiso, 'PlantCode', {'sum', 'mean'}, 'DataVars', {'NameplateCapacity_MW_', 'Latitude', 'Longitude'});
grpstats(T_eia860_caiso, 'PlantCode', {'sum', 'mean'}, 'DataVars', {'NameplateCapacity_MW_', 'Latitude', 'Longitude'});
T_pv.Properties.VariableNames{'sum_NameplateCapacity_MW_'} = 'TotalCapacity';
T_pv.Properties.VariableNames{'mean_Latitude'} = 'Latitude';
T_pv.Properties.VariableNames{'mean_Longitude'} = 'Longitude';

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


%% III. Now, write the 340 site data
cd('/home/bxl180002/scratch/SF2/');

ind_sites340 = T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'ind'};

tmp = reshape(ar_kcs_mean, size(ar_kcs_mean, 1)*size(ar_kcs_mean, 2), size(ar_kcs_mean, 3));
tmp = tmp(ind_sites340, :)'; % First dimension: Time, second dimension: Number of cells that include PV plants
h5create('sites340.h5','/kcs_mean',size(tmp));
h5write('sites340.h5', '/kcs_mean', tmp);

tmp = reshape(ar_kcs_std, size(ar_kcs_std, 1)*size(ar_kcs_std, 2), size(ar_kcs_std, 3));
tmp = tmp(ind_sites340, :)'; % First dimension: Time, second dimension: Number of cells that include PV plants
h5create('sites340.h5','/kcs_std',size(tmp));
h5write('sites340.h5', '/kcs_std', tmp);

tmp = reshape(ar_kcs_vrb, size(ar_kcs_vrb, 1)*size(ar_kcs_vrb, 2), size(ar_kcs_vrb, 3));
tmp = tmp(ind_sites340, :)'; % First dimension: Time, second dimension: Number of cells that include PV plants
h5create('sites340.h5','/kcs_vrb',size(tmp));
h5write('sites340.h5', '/kcs_vrb', tmp);

tmp = reshape(ar_kcs_wmean, size(ar_kcs_wmean, 1)*size(ar_kcs_wmean, 2), size(ar_kcs_wmean, 3));
tmp = tmp(ind_sites340, :)'; % First dimension: Time, second dimension: Number of cells that include PV plants
h5create('sites340.h5','/kcs_wmean',size(tmp));
h5write('sites340.h5', '/kcs_wmean', tmp);

h5create('sites340.h5','/timeutc_posix',size(ar_hourstart_posix));
h5write('sites340.h5', '/timeutc_posix', ar_hourstart_posix);

tmp = T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'cell_lat'};
h5create('sites340.h5','/lat_center',size(tmp));
h5write('sites340.h5', '/lat_center', tmp);

tmp = T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'cell_lon'};
h5create('sites340.h5','/lon_center',size(tmp));
h5write('sites340.h5', '/lon_center', tmp);

tmp = T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'TotalCapacity'};
h5create('sites340.h5','/cell_capacity',size(tmp));
h5write('sites340.h5', '/cell_capacity', tmp);

cd(dirhome);


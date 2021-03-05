write_flag = input('Permission to write? (1 - Yes/0 - No): ');
dirhome = pwd;
if ispc
    cell_dircontent = {...
        'C:\Users\bxl180002\Downloads\RampSolar\IBM\CA_square_daily_202002\content', ...
        'C:\Users\bxl180002\Downloads\RampSolar\IBM\CA_square_daily_202003\content', ...
        'C:\Users\bxl180002\Downloads\RampSolar\IBM\CA_square_daily_202004\content' ...
        };
elseif isunix
    cell_dircontent = {...
        '/home/bxl180002/scratch/SF2/IBM/CA_square_daily_202002/content', ...
        '/home/bxl180002/scratch/SF2/IBM/CA_square_daily_202003/content', ...
        '/home/bxl180002/scratch/SF2/IBM/CA_square_daily_202004/content' ...
        };
end

cell_datetime = cell(3, 1);
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

% All time stamps in April 2020
ar_datetime = [];
for d = 1:30
    datetime_start = datetime(2020, 4, d, 9, 0, 0, 'TimeZone', 'UTC');
    datetime_end = datetime(2020, 4, d, 2, 50, 0, 'TimeZone', 'UTC') + duration(24, 0, 0);
    time_per_day = [datetime_start: deltat: datetime_end]';
    ar_datetime = [ar_datetime; time_per_day];
end
cell_datetime{3} = ar_datetime;

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
ar_datetime = [cell_datetime{1}; cell_datetime{2}; cell_datetime{3}];
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

%%
% Calculate clear-sky GHI
add_pvlib();

ghi_cs = nan(320, 480, numel(ar_datetime));

lat_edge = R.LatitudeLimits(1): R.CellExtentInLatitude: R.LatitudeLimits(2);
lon_edge = R.LongitudeLimits(1): R.CellExtentInLongitude: R.LongitudeLimits(2);
lat_edge = lat_edge(:);
lon_edge = lon_edge(:);
lat_center = (lat_edge(1:end-1) + lat_edge(2:end))/2;
lon_center = (lon_edge(1:end-1) + lon_edge(2:end))/2;

if false
    save('temp.mat', 'ar_datetime', 'ar_istiff', 'ar_quantiles', 'cell_rawdata', 'lat_center', 'lon_center', '-v7.3');
end

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

if false
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
end

cd(dirhome);

%% CA poly, find out the coordinates of the boundary of the the minimum square 
fname = 'Solar forecast (PAIRS Team) (select sites)-GHI Quantile Forecast-horizon_60_quantile_50[solar_forecast_median]-02_15_2020T20_00_00.tiff';
[A1, R1] = readgeoraster(fname); % Example
A1 = flipud(A1); % I think IBM has use the CA polygon upside down
info = georasterinfo(fname);
m = info.MissingDataIndicator;
A1 = standardizeMissing(A1,m); 
figure();
imagesc(R1.LongitudeLimits(:), R1.LatitudeLimits(:), A1);
set(gca,'YDir', 'normal'); % Set origin to the bottom left corner
colormap jet;
title('Note the GHI forecast is upside down. Time: 12pm PST 2/15/20');

lat_edge_1 = R1.LatitudeLimits(1): R1.CellExtentInLatitude: R1.LatitudeLimits(2);
lon_edge_1 = R1.LongitudeLimits(1): R1.CellExtentInLongitude: R1.LongitudeLimits(2);
lat_edge_1 = lat_edge_1(:); 
lon_edge_1 = lon_edge_1(:);
lat_center_1 = (lat_edge_1(1:end-1) + lat_edge_1(2:end))/2;
lon_center_1 = (lon_edge_1(1:end-1) + lon_edge_1(2:end))/2;

% The boundary of the smallest square that contains CA 
BOUNDARY_CA = struct('WEST', nan, 'EAST', nan, 'SOUTH', nan, 'NORTH', nan); 
BOUNDARY_CA.NORTH = lat_center_1(find(any(~isnan(A1), 2), 1, 'last'));
BOUNDARY_CA.SOUTH = lat_center_1(find(any(~isnan(A1), 2), 1, 'first')); 
BOUNDARY_CA.WEST  = lon_center_1(find(any(~isnan(A1), 1), 1, 'first'));
BOUNDARY_CA.EAST  = lon_center_1(find(any(~isnan(A1), 1), 1, 'last'));

% The boundary of the minimum non-nan square
BOUNDARY = struct('WEST', nan, 'EAST', nan, 'SOUTH', nan, 'NORTH', nan); 
BOUNDARY.NORTH = lat_center(find(any(~isnan(A), 2), 1, 'last'));
BOUNDARY.SOUTH = lat_center(find(any(~isnan(A), 2), 1, 'first')); 
BOUNDARY.WEST  = lon_center(find(any(~isnan(A), 1), 1, 'first'));
BOUNDARY.EAST  = lon_center(find(any(~isnan(A), 1), 1, 'last'));

% Show the cropped CA picture
figure();
imagesc(...
    lon_center_1((lon_center_1<=BOUNDARY_CA.EAST)&(lon_center_1>=BOUNDARY_CA.WEST)), ...
    lat_center_1((lat_center_1<=BOUNDARY_CA.NORTH)&(lat_center_1>=BOUNDARY_CA.SOUTH)), ...
    A1((lat_center_1<=BOUNDARY_CA.NORTH)&(lat_center_1>=BOUNDARY_CA.SOUTH), (lon_center_1<=BOUNDARY_CA.EAST)&(lon_center_1>=BOUNDARY_CA.WEST)) ...
    ); 
set(gca,'YDir', 'normal'); % Set origin to the bottom left corner
colormap jet;
title('Note the GHI forecast is upside down. Time: 12pm PST 2/15/20');

% Now, apply both the CA and the non-nan boundary to the full square
cropped_lat = (lat_center<=BOUNDARY_CA.NORTH)&(lat_center>=BOUNDARY_CA.SOUTH)&(lat_center<=BOUNDARY.NORTH)&(lat_center>=BOUNDARY.SOUTH);
cropped_lon = (lon_center<=BOUNDARY_CA.EAST)&(lon_center>=BOUNDARY_CA.WEST)&(lon_center<=BOUNDARY.EAST)&(lon_center>=BOUNDARY.WEST);

% Show the cropped full picture for visual inspection
figure();
A_tmp = squeeze(cell_rawdata{3}(:, :, 67)); % 67 is 12pm PST 2020-02-01
imagesc(lon_center(cropped_lon), lat_center(cropped_lat), A_tmp(cropped_lat, cropped_lon, :)); 
set(gca,'YDir', 'normal'); % Set origin to the bottom left corner
colormap jet;

% Extract the non-nan CA polygon and represent it using a matrix of logical values
is_california_A1 = A1;
is_california_A1(isnan(A1(:))) = 0;
is_california_A1(~isnan(A1(:))) = 1;
is_california_A = [is_california_A1 zeros(size(is_california_A1,  1), size(A, 2) - size(is_california_A1, 2))];
is_nonnan_A = zeros(size(A));
is_nonnan_A(cropped_lat, cropped_lon) = 1;
is_nonnan_ca_A = is_nonnan_A.*is_california_A;

% Show the non-nan CA polygon and compare it with Matlab's built-in CA
% polygon
figure();
axesm('mercator','Grid', 'on', 'MapLonLimit',[-125 -113], 'MapLatLimit',[32 43], 'GLineStyle', '.', 'Gcolor', [.5, .5, .5]);
alabamahi = shaperead('usastatehi', 'UseGeoCoords', true,...
'Selector',{@(name) strcmpi(name,'California'), 'Name'});
tmp = R;
tmp.ColumnsStartFrom = 'south';
latlim = tmp.LatitudeLimits;
lonlim = tmp.LongitudeLimits;
usamap(latlim, lonlim);
geoshow(is_nonnan_ca_A, tmp, 'DisplayType', 'surface', 'FaceAlpha', 0.5);
geoshow(alabamahi, 'FaceColor', [1, 0, 0]);
tightmap;


%% II. Now, write the cropped square 
cd('/home/bxl180002/scratch/SF2/');

if write_flag
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
end

cd(dirhome);
%% II.2 Write the cropped square, only keep CA region as non-nan
cd('/home/bxl180002/scratch/SF2/');

if write_flag
    ar_kcs_mean_croppedca = nan(size(ar_kcs_mean));
    for i = 1: size(ar_kcs_mean_croppedca, 3)
        layer = ar_kcs_mean(:, :, i);
        layer(is_nonnan_ca_A == 0) = nan;
        ar_kcs_mean_croppedca(:, :, i) = layer;
    end
    ar_kcs_mean_croppedca = ar_kcs_mean_croppedca(cropped_lat, cropped_lon, :);
    h5create('cropped_ca.h5','/kcs_mean',size(ar_kcs_mean_croppedca));
    h5write('cropped_ca.h5', '/kcs_mean', ar_kcs_mean_croppedca);

    ar_kcs_std_croppedca = nan(size(ar_kcs_std));
    for i = 1: size(ar_kcs_std_croppedca, 3)
        layer = ar_kcs_std(:, :, i);
        layer(is_nonnan_ca_A == 0) = nan;
        ar_kcs_std_croppedca(:, :, i) = layer;
    end
    ar_kcs_std_croppedca = ar_kcs_std_croppedca(cropped_lat, cropped_lon, :);
    h5create('cropped_ca.h5','/kcs_std',size(ar_kcs_std_croppedca));
    h5write('cropped_ca.h5', '/kcs_std', ar_kcs_std_croppedca);

    ar_kcs_vrb_croppedca = nan(size(ar_kcs_vrb));
    for i = 1: size(ar_kcs_vrb_croppedca, 3)
        layer = ar_kcs_vrb(:, :, i);
        layer(is_nonnan_ca_A == 0) = nan;
        ar_kcs_vrb_croppedca(:, :, i) = layer;
    end
    ar_kcs_vrb_croppedca = ar_kcs_vrb_croppedca(cropped_lat, cropped_lon, :);
    h5create('cropped_ca.h5','/kcs_vrb',size(ar_kcs_vrb_croppedca));
    h5write('cropped_ca.h5', '/kcs_vrb', ar_kcs_vrb_croppedca);

    ar_kcs_wmean_croppedca = nan(size(ar_kcs_wmean));
    for i = 1: size(ar_kcs_wmean_croppedca, 3)
        layer = ar_kcs_wmean(:, :, i);
        layer(is_nonnan_ca_A == 0) = nan;
        ar_kcs_wmean_croppedca(:, :, i) = layer;
    end
    ar_kcs_wmean_croppedca = ar_kcs_wmean_croppedca(cropped_lat, cropped_lon, :);
    h5create('cropped_ca.h5','/kcs_wmean',size(ar_kcs_wmean_croppedca));
    h5write('cropped_ca.h5', '/kcs_wmean', ar_kcs_wmean_croppedca);

    h5create('cropped_ca.h5','/timeutc_posix',size(ar_hourstart_posix));
    h5write('cropped_ca.h5', '/timeutc_posix', ar_hourstart_posix);

    h5create('cropped_ca.h5','/lat_center',size(lat_center(cropped_lat)));
    h5write('cropped_ca.h5', '/lat_center', lat_center(cropped_lat));

    h5create('cropped_ca.h5','/lon_center',size(lon_center(cropped_lon)));
    h5write('cropped_ca.h5', '/lon_center', lon_center(cropped_lon));
end

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

if write_flag
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
end

cd(dirhome);

%% IV. Load and write FRP requirements in hourly resolution
load_caiso_data;
ar_datetime_hour = datetime(ar_datetime.Year, ar_datetime.Month, ar_datetime.Day, ar_datetime.Hour, 0, 0, 'TimeZone', 'UTC');
T_rtpd = T_rtpd((T_rtpd{:, 'HOUR_START'}>=min(ar_datetime))&(T_rtpd{:, 'HOUR_START'}<=max(ar_datetime)), :); % Only look at the hours with raster forecast

dirhome = pwd;
cd('/home/bxl180002/scratch/SF2/');

% for i = 1: length(cell_pvcellghi)
%     T_ghi_forcong = [array2table(ar_datetime, 'VariableNames', {'TIME'}) array2table(cell_pvcellghi{i}, 'VariableNames', cellstr(num2str(T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'ind'})))];
%     writetable(T_ghi_forcong, strcat(num2str(i), '.csv'));
% end
T_rtpd_forcong = T_rtpd(:, {'HOUR_START', 'error_max', 'error_min'});
T_rtpd_forcong = grpstats(T_rtpd_forcong, 'HOUR_START', {'max', 'min'}, 'DataVars', {'error_max', 'error_min'});
T_rtpd_forcong = T_rtpd_forcong(:, {'HOUR_START', 'max_error_max', 'min_error_min'});
T_rtpd_forcong.Properties.VariableNames{'max_error_max'} = 'FRU';
T_rtpd_forcong.Properties.VariableNames{'min_error_min'} = 'FRD';

if write_flag
    writetable(T_rtpd_forcong, 'rtpd.csv');
    writetable(T_pvcell, 'T_pvcell.csv');
end
cd(dirhome);


% Save GHI of the 340 cells
cell_pvcellghi = cell(5, 1);
for j = 1: 5
    tmp = reshape(cell_rawdata{j}, size(cell_rawdata{j}, 1)*size(cell_rawdata{j}, 2), size(cell_rawdata{j}, 3));
    cell_pvcellghi{j} = tmp(T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'ind'}, :)'; % First dimension: Time, second dimension: Number of cells that include PV plants
end

% Save clear-sky GHI of the 340 cells
tmp = reshape(ghi_cs, size(ghi_cs, 1)*size(ghi_cs, 2), size(ghi_cs, 3));
ghi_cs_cell = tmp(T_pvcell{T_pvcell{:, 'within_boundary'}==1, 'ind'}, :)'; % First dimension: Time, second dimension: Number of cells that include PV plants

% Find the locations of all PV plant
T_pv{:, 'ind'} = sub2ind([size(A, 1), size(A, 2)], T_pv{:, 'ny'}, T_pv{:, 'nx'}); % Sub to index
T_eia860_caiso{:, 'nx'} = nan;
T_eia860_caiso{:, 'ny'} = nan;
for i = 1: height(T_eia860_caiso)
    ar_nx = find(abs(lon_center - T_eia860_caiso{i, 'Longitude'}) < R1.CellExtentInLongitude/2);
    if numel(ar_nx) == 1
        T_eia860_caiso{i, 'nx'} = ar_nx;
    else
        disp('WARNING: CANNOT LOCATE THIS PV PLANT: LONGITUDE!');
    end
    ar_ny = find(abs(lat_center(:) - T_eia860_caiso{i, 'Latitude'}) < R1.CellExtentInLatitude/2); 
    if numel(ar_ny) == 1
        T_eia860_caiso{i, 'ny'} = ar_ny;
    else
        disp('WARNING: CANNOT LOCATE THIS PV PLANT: LATITUDE!');
    end
end
T_eia860_caiso{:, 'ind'} = sub2ind([size(A, 1), size(A, 2)], T_eia860_caiso{:, 'ny'}, T_eia860_caiso{:, 'nx'}); % Sub to index
T_eia860_caiso{:, 'within_boundary'} = (T_eia860_caiso{:, 'Longitude'}>=BOUNDARY.WEST) & (T_eia860_caiso{:, 'Longitude'}<=BOUNDARY.EAST) & (T_eia860_caiso{:, 'Latitude'}>=BOUNDARY.SOUTH) & (T_eia860_caiso{:, 'Latitude'}<=BOUNDARY.NORTH);

%% V. Calculate PV power generation of each of the 340 cells
% Calculate PV power generation
dayofyear = pvl_date2doy(Time.year, Time.month, Time.day);
nT = size(ar_datetime, 1);
p_solar = nan(nT, height(T_eia860_caiso), 5); % Power output
p_solar_cs = nan(nT, height(T_eia860_caiso)); % Clear-sky power output

%Other parameters
SF=0.98;
%Weather
PresPa=101325;
WIND=0;
dryT=10;
Albedo = 0.2;

modulepv=134; %Topaz uses first solar panels (FS272 is probably the oldest they have)
inverterid=759; 
modules_Series=11;
modules_parallel=2360;
ninvert=149;

tic;
for i = 1: height(T_eia860_caiso)

    Location = pvl_makelocationstruct(T_eia860_caiso.Latitude(i), T_eia860_caiso.Longitude(i));
    %--- SOLAR FARM SPECS---
    %Define module
    ModuleParameters = pvl_sapmmoduledb(modulepv,'SandiaModuleDatabase_20120925.xlsx');
    %Define the inverter
    load('SandiaInverterDatabaseSAM2014.1.14.mat')
    Inverter = SNLInverterDB(inverterid);
    %Topaz uses power one inverters
    clear InverterNames SNLInverterDB

    %Define the array configuration
    if isnan(T_eia860_caiso.TiltAngle(i))
        Array.Tilt = 0;
    else
        Array.Tilt = T_eia860_caiso.TiltAngle(i); % Array tilt angle (deg)
    end
    if isnan(T_eia860_caiso.AzimuthAngle(i))
        Array.Azimuth = 180;
    else
        Array.Azimuth = T_eia860_caiso.AzimuthAngle(i); %Array azimuth (180 deg indicates array faces South)
    end
    Array.Ms = modules_Series; %Number of modules in series
    Array.Mp = modules_parallel; %Number of paralell strings  
    %Location of site
    Array.a = -3.56;
    Array.b = -0.075;

    [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
    if strcmp(T_eia860_caiso.Single_AxisTracking_(i), 'Y')
        % Fixed tilt
        MaxAngle = 90;
        [~, AOI, ~, ~] = pvl_singleaxis(90-AppSunEl, SunAz, Location.latitude, Array.Tilt, Array.Azimuth, MaxAngle);
    else
        AOI = pvl_getaoi(Array.Tilt, Array.Azimuth, 90-AppSunEl, SunAz);
    end

    Wspd=WIND.*ones(nT, 1);
    Drybulb=dryT.*ones(nT, 1);
    AMa = pvl_absoluteairmass(pvl_relativeairmass(90-AppSunEl),PresPa);
    F1 = max(0,polyval(ModuleParameters.a,AMa)); %Spectral loss function
    F2 = max(0,polyval(ModuleParameters.b,AOI)); % Angle of incidence loss function

    for k = 1: 5 % quantile
        if ~T_eia860_caiso{i, 'within_boundary'}
            p_solar(:, i, k) = 0;
            continue;
        else
            ghi = squeeze(cell_rawdata{k}(T_eia860_caiso{i, 'ny'}, T_eia860_caiso{i, 'nx'}, :)); % ny is the first dimension, nx is the second
            i_ghi_notnan = ~isnan(ghi);
            EdiffGround = nan(nT, 1);
            EdiffGround(i_ghi_notnan) = pvl_grounddiffuse(Array.Tilt, ghi(i_ghi_notnan), Albedo);
            DNI_model = nan(nT, 1);
            DNI_model(i_ghi_notnan) = pvl_disc(ghi(i_ghi_notnan),90-SunEl(i_ghi_notnan), dayofyear(i_ghi_notnan),PresPa);
            DHI_model = nan(nT, 1);
            DHI_model(i_ghi_notnan) = ghi(i_ghi_notnan) - cosd(90-SunEl(i_ghi_notnan)).*DNI_model(i_ghi_notnan);
            Eb = 0*AOI; %Initiallize variable
            Eb(AOI<90) = DNI_model(AOI<90).*cosd(AOI(AOI<90)); %Only calculate when sun is in view of the plane of array
            EdiffSky = nan(nT, 1);
            EdiffSky(i_ghi_notnan) = pvl_isotropicsky(Array.Tilt,DHI_model(i_ghi_notnan));
            E = Eb + EdiffSky + EdiffGround; % Total incident irradiance (W/m^2)
            E0 = 1000; %Reference irradiance (1000 W/m^2)
            celltemp = nan(nT, 1);
            celltemp(i_ghi_notnan) = pvl_sapmcelltemp(E(i_ghi_notnan), E0, Array.a, Array.b,Wspd((i_ghi_notnan),1), Drybulb((i_ghi_notnan),1), ModuleParameters.delT);
            Ediff = EdiffSky + EdiffGround; % Total diffuse incident irradiance (W/m^2)
            Ee = F1.*((Eb.*F2+ModuleParameters.fd.*Ediff)/E0)*SF; %Effective irradiance
            % Ee(isnan(Ee))=0; % Set any NaNs to zero
            mSAPMResults = pvl_sapm(ModuleParameters, Ee(i_ghi_notnan), celltemp(i_ghi_notnan));
            aSAPMResults.Vmp = nan(nT, 1);
            aSAPMResults.Imp = nan(nT ,1); 
            aSAPMResults.Pmp = nan(nT, 1);
            aSAPMResults.Vmp(i_ghi_notnan) = Array.Ms  *mSAPMResults.Vmp;
            aSAPMResults.Imp(i_ghi_notnan) = Array.Mp  *mSAPMResults.Imp;
            aSAPMResults.Pmp = aSAPMResults.Vmp .* aSAPMResults.Imp;
            % clear temp
            % temp=find(ghi~=10000);
            % ACPower(1:size(Time.hour,1), 7) =10000;
            ACPower = nan(nT, 1);
            ACPower(i_ghi_notnan)= pvl_snlinverter(Inverter, mSAPMResults.Vmp*Array.Ms, mSAPMResults.Pmp*Array.Ms*Array.Mp)*ninvert/1000000; % MW
            ACPower(ACPower<0, end)=0;
            p_solar(:, i, k) = ACPower;
        end
    end

    % Now the clear-sky power
    ghi = squeeze(ghi_cs(T_eia860_caiso{i, 'ny'}, T_eia860_caiso{i, 'nx'}, :)); % ny is the first dimension, nx is the second
    EdiffGround = pvl_grounddiffuse(Array.Tilt, ghi, Albedo);
    DNI_model = pvl_disc(ghi,90-SunEl, dayofyear,PresPa);
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
    % ACPower(1:size(Time.hour,1), 7) =10000;
    ACPower= pvl_snlinverter(Inverter, mSAPMResults.Vmp(temp)*Array.Ms, mSAPMResults.Pmp(temp)*Array.Ms*Array.Mp)*ninvert/1000000; % MW
    ACPower(ACPower<0, end)=0;
    p_solar_cs(:, i) = ACPower;
end
toc;

% Save breakpoint file
cd('/home/bxl180002/scratch/SF2/');
if write_flag
    save('breakpoint.mat', 'cell_pvcellghi', 'ghi_cs_cell', 'T_eia860_caiso', 'T_pv', 'T_pvcell', 'ar_datetime', 'A', 'A1', 'R', 'R1', 'lat_center', 'lon_center', 'lat_edge', 'p_solar', 'p_solar_cs', 'is_nonnan_ca_A');
end
cd(dirhome);


% %% Baseline run
% karray = 5:5:40; % the number of nearest neighbors
% ar_datetime_hour_unique = unique(ar_datetime_hour);
% testhour_end = datetime(2020, 4, 1, 2, 0, 0, 'TimeZone', 'UTC');
% testhour_start = datetime(2020, 3, 25, 9, 0, 0, 'TimeZone', 'UTC');
% T_baseline_rtpd = array2table(...
%     ar_datetime_hour_unique((ar_datetime_hour_unique>=testhour_start)&(ar_datetime_hour_unique<=testhour_end)), ...
%     'VariableNames', ...
%     {'HOUR_START'}); % Result container
% 
% fprintf('Baseline\n');
% for k = karray
%     T_baseline_rtpd{:, strcat('FRU_k', num2str(k))} = nan;
%     T_baseline_rtpd{:, strcat('FRD_k', num2str(k))} = nan;
%     for i = 1: size(T_baseline_rtpd, 1)
%         this_date = datetime(T_baseline_rtpd.HOUR_START.Year(i), T_baseline_rtpd.HOUR_START.Month(i), T_baseline_rtpd.HOUR_START.Day(i), 'TimeZone', 'UTC');
%         this_hour = T_baseline_rtpd.HOUR_START.Hour(i);
%         selected_days = (T_rtpd.TIME_START<this_date)&(T_rtpd.TIME_START>=this_date-days(k))&(T_rtpd.HOUR_START.Hour==this_hour); % We use 30 previous days
%         sample_error_max = T_rtpd{selected_days, 'error_max'};
%         sample_error_min = T_rtpd{selected_days, 'error_min'};
%         [f,x] = ecdf(sample_error_max(:));
%         T_baseline_rtpd{i, strcat('FRU_k', num2str(k))} = interp1(f, x, 0.975);
%         [f,x] = ecdf(sample_error_min(:));
%         T_baseline_rtpd{i, strcat('FRD_k', num2str(k))} = interp1(f, x, 0.025);
%     end
% 
%     fprintf('k = %g\n', k);
% end
% plot(T_baseline_rtpd{:, 'HOUR_START'}, T_baseline_rtpd{:, 2:end}, 'b')
% hold on;
% plot(...
%     T_rtpd_forcong{(T_rtpd_forcong{:, 'HOUR_START'}>=testhour_start)&(T_rtpd_forcong{:, 'HOUR_START'}<=testhour_end), 'HOUR_START'}, ...
%     T_rtpd_forcong{(T_rtpd_forcong{:, 'HOUR_START'}>=testhour_start)&(T_rtpd_forcong{:, 'HOUR_START'}<=testhour_end), 2:end}, ...
%     'r' ...
%     );
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Calculate classifier based on kpv
% sigma_psolar = nan(size(p_solar, 1), 5);
% % plot(sum(squeeze(p_solar(:, T_eia860_caiso{:, 'within_boundary'}==1, 3)*diag(T_eia860_caiso{T_eia860_caiso{:, 'within_boundary'}==1, 'NameplateCapacity_MW_'})), 2));
% for q = 1:5
%     sigma_psolar(:, q) = sum(squeeze(p_solar(:, T_eia860_caiso{:, 'within_boundary'}==1, 3)*diag(T_eia860_caiso{T_eia860_caiso{:, 'within_boundary'}==1, 'NameplateCapacity_MW_'})), 2);
% end
% sigma_psolarcs = sum(p_solar_cs(:, T_eia860_caiso{:, 'within_boundary'}==1)*diag(T_eia860_caiso{T_eia860_caiso{:, 'within_boundary'}==1, 'NameplateCapacity_MW_'}), 2);
% 
% kpv = sigma_psolar./repmat(sigma_psolarcs, 1, 5);
% 
% T_classifier_kpv = array2table(ar_datetime_hour_unique, 'VariableNames', {'HOUR_START'});
% T_classifier_kpv{:, 'mean_kpv'} = nan; % Mean of p50 kpv
% T_classifier_kpv{:, 'std_kpv'} = nan; % Standard deviation of p50 kpv
% T_classifier_kpv{:, 'vrb_kpv'} = nan; % Variability of p50 kpv
% T_classifier_kpv{:, 'wmean_kpv'} = nan; % Mean of p75 - p25
% for i = 1: height(T_classifier_kpv)
%     kpv_1h_p25 = kpv(ar_datetime_hour==T_classifier_kpv{i, 'HOUR_START'}, 2);
%     kpv_1h_p50 = kpv(ar_datetime_hour==T_classifier_kpv{i, 'HOUR_START'}, 3);
%     kpv_1h_p75 = kpv(ar_datetime_hour==T_classifier_kpv{i, 'HOUR_START'}, 4);
% 
%     T_classifier_kpv{i, 'mean_kpv'} = mean(kpv_1h_p50, 'omitnan');
%     T_classifier_kpv{i, 'std_kpv'} = std(kpv_1h_p50, 0, 'omitnan');
%     T_classifier_kpv{i, 'vrb_kpv'} = sqrt(1/5.*sum(diff(kpv_1h_p50).^2));
%     T_classifier_kpv{i, 'wmean_kpv'} = mean(kpv_1h_p75 - kpv_1h_p25, 'omitnan');
% end
% T_classifier_kpv.DATE = datetime(T_classifier_kpv.HOUR_START.Year, T_classifier_kpv.HOUR_START.Month, T_classifier_kpv.HOUR_START.Day, 'TimeZone', 'UTC');
% 
% fprintf('kNN\n');
% cell_knnkpv_rtpd = cell(4, 1);
% for c = 1:4
%     T_knnkpv_rtpd = T_classifier_kpv((ar_datetime_hour_unique>=testhour_start)&(ar_datetime_hour_unique<=testhour_end), :); % Result container
%     switch c
%         case 1
%             T_classifier_kpv{:, 'c'} = T_classifier_kpv{:, 'mean_kpv'};
%             T_knnkpv_rtpd{:, 'c'}    = T_knnkpv_rtpd{:, 'mean_kpv'};
%         case 2
%             T_classifier_kpv{:, 'c'} = T_classifier_kpv{:, 'std_kpv'};
%             T_knnkpv_rtpd{:, 'c'}    = T_knnkpv_rtpd{:, 'std_kpv'};
%         case 3
%             T_classifier_kpv{:, 'c'} = T_classifier_kpv{:, 'vrb_kpv'};
%             T_knnkpv_rtpd{:, 'c'}    = T_knnkpv_rtpd{:, 'vrb_kpv'};
%         case 4
%             T_classifier_kpv{:, 'c'} = T_classifier_kpv{:, 'wmean_kpv'};
%             T_knnkpv_rtpd{:, 'c'}    = T_knnkpv_rtpd{:, 'wmean_kpv'};
%     end
%     for k = karray
%         T_knnkpv_rtpd{:, strcat('FRU_k', num2str(k))} = nan;
%         T_knnkpv_rtpd{:, strcat('FRD_k', num2str(k))} = nan;
% 
%         for i = 1: size(T_knnkpv_rtpd, 1)
%             this_date = datetime(T_knnkpv_rtpd.HOUR_START.Year(i), T_knnkpv_rtpd.HOUR_START.Month(i), T_knnkpv_rtpd.HOUR_START.Day(i), 'TimeZone', 'UTC');
%             this_hour = T_knnkpv_rtpd.HOUR_START.Hour(i);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Can use the days AFTER this_date if results are not ideal
%             T_sample = T_classifier_kpv((T_classifier_kpv.DATE<this_date)&(T_classifier_kpv.HOUR_START.Hour==this_hour), :);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if any(isnan(T_knnkpv_rtpd{i, 'c'}))
%                 T_sample_sorted = T_sample(T_sample.DATE>=this_date-days(k), :); % We use 30 previous days, i.e., baseline, if data is nan
%                 selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); % 30 the nearest days for baseline
%             else
%                 T_sample.dist = sqrt(sum((T_sample{:, 'c'}-T_knnkpv_rtpd{i, 'c'}).^2, 2)); % Euclidean distance
%                 T_sample_sorted = sortrows(T_sample, 'dist');
%                 selected_days = ismember(datetime(T_rtpd.TIME_START.Year, T_rtpd.TIME_START.Month, T_rtpd.TIME_START.Day, 'TimeZone', 'UTC'), T_sample_sorted.DATE(1:k)); 
%             end
% 
%             sample_error_max = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_max'};
%             sample_error_min = T_rtpd{selected_days&(T_rtpd.TIME_START.Hour==this_hour), 'error_min'};
%             [f,x] = ecdf(sample_error_max(:));
%             T_knnkpv_rtpd{i, strcat('FRU_k', num2str(k))} = interp1(f, x, 0.975);
%             [f,x] = ecdf(sample_error_min(:));
%             T_knnkpv_rtpd{i, strcat('FRD_k', num2str(k))} = interp1(f, x, 0.025);
%         end
%         fprintf('c = %g, k = %g\n', c, k);
%     end
%     cell_knnkpv_rtpd{c} = T_knnkpv_rtpd;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Process results and visualize
% f_errormax_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_max'}; % 15-min
% f_errormin_rtpd = T_rtpd{ismember(T_rtpd.HOUR_START, T_baseline_rtpd.HOUR_START), 'error_min'}; % 15-min
% T_errormax_rtpd = array2table(reshape(f_errormax_rtpd, 4, numel(f_errormax_rtpd)/4)', 'VariableNames', {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'});
% T_errormin_rtpd = array2table(reshape(f_errormin_rtpd, 4, numel(f_errormin_rtpd)/4)', 'VariableNames', {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'});
% T_actual_error = [T_baseline_rtpd(:, 'HOUR_START') T_errormax_rtpd T_errormin_rtpd];
% ar_actual_error_fru = T_actual_error{:, {'error_max_1', 'error_max_2', 'error_max_3', 'error_max_4'}};
% ar_actual_error_frd = T_actual_error{:, {'error_min_1', 'error_min_2', 'error_min_3', 'error_min_4'}};
% 
% %% Baseline
% for k = karray
%     [oversupply_baseline_fru(karray==k), pshort_baseline_fru(karray==k)] = calculate_oversupply_and_pshort(ar_actual_error_fru, T_baseline_rtpd{:, strcat('FRU_k', num2str(k))});
%     [oversupply_baseline_frd(karray==k), pshort_baseline_frd(karray==k)] = calculate_oversupply_and_pshort(-ar_actual_error_frd, -T_baseline_rtpd{:, strcat('FRD_k', num2str(k))});
% end
% 
% % kNN
% for c = 1:4
%     T_knnkpv_rtpd = cell_knnkpv_rtpd{c};
%     for k = karray
%         [oversupply_knnkpv_fru(karray==k, c), pshort_knnkpv_fru(karray==k, c)] = calculate_oversupply_and_pshort(ar_actual_error_fru, T_knnkpv_rtpd{:, strcat('FRU_k', num2str(k))});
%         [oversupply_knnkpv_frd(karray==k, c), pshort_knnkpv_frd(karray==k, c)] = calculate_oversupply_and_pshort(-ar_actual_error_frd, -T_knnkpv_rtpd{:, strcat('FRD_k', num2str(k))});
%     end
% end
% 
% figure(); plot(karray, oversupply_baseline_frd, 'k', karray, oversupply_knnkpv_frd); title('Oversupply (FRU)');
% figure(); plot(karray, pshort_baseline_frd, 'k', karray, pshort_knnkpv_frd); title('Pshort (FRU)');
% figure(); plot(karray, oversupply_baseline_fru, 'k', karray, oversupply_knnkpv_fru); title('Oversupply (FRD)');
% figure(); plot(karray, pshort_baseline_fru, 'k', karray, pshort_knnkpv_fru); title('Pshort (FRD)');
% 
% 
% % Pareto frontier
% figure(); 
% h(5)= scatter(pshort_baseline_fru, oversupply_baseline_fru./1E3, 60, 'ko');
% hold on;
% for c = 1: 4
%     h(c) = scatter(pshort_knnkpv_fru(:, c), oversupply_knnkpv_fru(:, c)./1E3, 30, '^');
% end
% set(h(1), 'MarkerFaceColor', 'r');
% set(h(2), 'MarkerFaceColor', 'g');
% set(h(3), 'MarkerFaceColor', 'm');
% set(h(4), 'MarkerFaceColor', 'b');
% set(h(5), 'MarkerFaceColor', 'k');
% title('FRU');
% xline(pshort_baseline_fru(karray==30), 'Color', 'r', 'LineStyle', '--');
% yline(oversupply_baseline_fru(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
% 
% figure(); 
% h(5)= scatter(pshort_baseline_frd, oversupply_baseline_frd./1E3, 60, 'ko');
% hold on;
% for c = 1: 4
%     h(c) = scatter(pshort_knnkpv_frd(:, c), oversupply_knnkpv_frd(:, c)./1E3, 30, '^');
% end
% set(h(1), 'MarkerFaceColor', 'r');
% set(h(2), 'MarkerFaceColor', 'g');
% set(h(3), 'MarkerFaceColor', 'm');
% set(h(4), 'MarkerFaceColor', 'b');
% set(h(5), 'MarkerFaceColor', 'k');
% title('FRD');
% xline(pshort_baseline_frd(karray==30), 'Color', 'r', 'LineStyle', '--');
% yline(oversupply_baseline_frd(karray==30)./1E3, 'Color', 'r', 'LineStyle', '--');
% 

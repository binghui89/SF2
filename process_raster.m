% CA poly example
fname = 'Solar forecast (PAIRS Team) (select sites)-GHI Quantile Forecast-horizon_60_quantile_50[solar_forecast_median]-02_15_2020T20_00_00.tiff';
[A1, R1] = readgeoraster(fname); % Example
info = georasterinfo(fname);
m = info.MissingDataIndicator;
A = standardizeMissing(A1,m);
figure();imagesc(A);colormap jet;

% Square example
fname = 'Solar forecast (PAIRS Team) (select sites)-GHI Quantile Forecast-horizon_60_quantile_50[solar_forecast_median]-02_01_2020T20_00_00.tiff';
[A1, R1] = readgeoraster(fname); % Example
info = georasterinfo(fname);
m = info.MissingDataIndicator;
A = standardizeMissing(A1,m);
figure();imagesc(A);colormap jet;

%% Check file availability and data validity
dirhome = pwd;
dirwork = 'C:\Users\bxl180002\Downloads\RampSolar\IBM\CA_square_daily_202002\content';

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

% Check file availability
ar_istiff = zeros(numel(ar_datetime), numel(ar_quantiles));
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
        else
            ar_istiff(i, j) = 0;
        end
    end
end
cd(dirhome);

ar_datetime_local = ar_datetime;
ar_datetime.TimeZone = 'America/Los_Angeles';
T_istiff = [array2table(ar_datetime, 'VariableNames', {'Time'}), array2table(ar_istiff, 'VariableNames', {'q5', 'q25', 'q50', 'q75', 'q95'})];
T_istiff_local = [array2table(ar_datetime_local, 'VariableNames', {'Time'}), array2table(ar_istiff, 'VariableNames', {'q5', 'q25', 'q50', 'q75', 'q95'})];

% Check quantile crossing
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

T_error = [array2table(ar_datetime, 'VariableNames', {'Time'}), array2table(ar_error, 'VariableNames', {'q5>q25', 'q50>q25', 'q75>q50', 'q95>q75'})];
T_error_local = [array2table(ar_datetime_local, 'VariableNames', {'Time'}), array2table(ar_error, 'VariableNames', {'q5>q25', 'q50>q25', 'q75>q50', 'q95>q75'})];


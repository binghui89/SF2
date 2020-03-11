function ghi_cs_mean = mean_clearsky_ghi(lat, lon, tarray, deltat)

Location = pvl_makelocationstruct(lat, lon);
Location.altitude = 0;
tarray_formean = repmat(tarray(:)', deltat, 1) - datenum(repmat([deltat:-1:1]'./24/60, 1, numel(tarray))); % All minutes during the past period, can be 15 or 5 min
tarray_formean = tarray_formean(:);
Time_formean.UTCOffset = zeros(size(tarray_formean, 1), 1); % tarray must be UTC time
% Time_formean.UTCOffset = tzoffset(tarray_formean); 
Time_formean.year   = tarray_formean.Year;
Time_formean.month  = tarray_formean.Month;
Time_formean.day    = tarray_formean.Day;
Time_formean.hour   = tarray_formean.Hour;
Time_formean.minute = tarray_formean.Minute;
Time_formean.second = tarray_formean.Second;

% [ghi_cs_formean, ~, ~] = pvl_clearsky_ineichen(Time_formean, Location); % Let's use Ineichen model since it is recommended by Inman et al.
[~, ~, AppSunEl_formean, ~] = pvl_ephemeris(Time_formean,Location);
ghi_cs_formean = pvl_clearsky_haurwitz(90-AppSunEl_formean);

ghi_cs_mean = mean(reshape(ghi_cs_formean, numel(ghi_cs_formean)/numel(tarray), numel(tarray)), 1)';

end
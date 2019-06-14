function M_fixed = fix_ibm(M, lat, lon, plot_flag)
% M columns: year, month, day, hour, minute, second, lower, mean, upper
% Time should be in UTC time

if nargin == 3
    plot_flag = false;
end

Location = pvl_makelocationstruct(lat, lon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix 1: Time shift
toff=2;
Mp = [M(1:size(M, 1)-toff, 1: 5), M(toff+1: end, 6:8)];
M = Mp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following line is for data before April, not useful because the data 
% in December, 2018 is terrible
% M(any(M(:, 6:8)>1E4, 2), 6:8) = 0; 

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

[SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
ghi_clearsky = pvl_clearsky_haurwitz(90-AppSunEl); % Clear-sky GHI

M(ghi_clearsky==0, 6:8) = 0;
kcs = M(:, 6:8)./repmat(ghi_clearsky, 1, 3);
kcs(ghi_clearsky==0, :) = 0;

% Note we are using PST all year round.
tarray = datetime(Time.year, Time.month, Time.day, Time.hour, Time.minute, Time.second, 'TimeZone', 'UTC');
tarray.TimeZone = '-08:00'; % Let's just use PST, so that 12:00 is noon
tarray_pst = datetime(tarray.Year, tarray.Month, tarray.Day, tarray.Hour, tarray.Minute, tarray.Second);


% % The following section is to validate my assumption using histgrams
% figure()
% subplot(3, 1, 1);
% kcs_selected = kcs((tarray_pst.Hour>=7) & (tarray_pst.Hour<17), :);
% hist(kcs_selected(:), 50);
% subplot(3, 1, 2);
% kcs_selected = kcs((tarray_pst.Hour<7) | (tarray_pst.Hour>=17), :);
% hist(kcs_selected(:), 50);
% subplot(3, 1, 3);
% h = plot(tarray_pst, M(:, 6:8));
% set(h, {'color'}, {'b'; 'k'; 'r'});
% hold on;
% plot(tarray_pst, ghi_clearsky, 'm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix 2: GHI rescaled by maximum clear-sky index between 7 to 17, capped by
% clear-sky GHI otherwise.
M_pst = [tarray_pst.Year, tarray_pst.Month, tarray_pst.Day, tarray_pst.Hour, tarray_pst.Minute, M(:, 6:8)];
days_unique = unique(M_pst(:, 1:3), 'rows');
for iday = 1: size(days_unique, 1)
    thisday = days_unique(iday, :);
    i_selected = (M_pst(:, 1) == thisday(1)) & (M_pst(:, 2) == thisday(2)) & (M_pst(:, 3) == thisday(3));
    i_scale    = i_selected & (M_pst(:, 4)>=7) & (M_pst(:, 4)<17);
    kcs_scale = M_pst(i_scale, 6:8)./repmat(ghi_clearsky(i_scale), 1, 3);
    kcs_max = max(kcs_scale(:));
    if kcs_max > 1
        M_pst(i_selected, 6:8) = M_pst(i_selected, 6:8)./kcs_max;
    end
    
    i_cap = i_selected & ((M_pst(:, 4)<7) | (M_pst(:, 4)>=17));
    rows_cap = find(i_cap);
    for j = 1:length(rows_cap)
        irow = rows_cap(j);
        ghi_max = max(M_pst(irow, 6:8));
        if ghi_max > ghi_clearsky(irow)
            if ghi_clearsky(irow) == 0
                M_pst(irow, 6:8) = 0;
            else
                M_pst(irow, 6:8) = M_pst(irow, 6:8).*ghi_clearsky(irow)/ghi_max;
            end
        end
    end
end
M_pst(ghi_clearsky==0, 6:8) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_fixed = [M(:, 1: 5), M_pst(:, 6:8)];

if plot_flag
    figure()
    subplot(4, 1, 1);
    kcs_selected = kcs((tarray_pst.Hour>=7) & (tarray_pst.Hour<17), :);
    hist(kcs_selected(:), 50);
    subplot(4, 1, 2);
    kcs_selected = kcs((tarray_pst.Hour<7) | (tarray_pst.Hour>=17), :);
    hist(kcs_selected(:), 50);
    subplot(4, 1, 3);
    h1 = plot(tarray_pst, M_fixed(:, 6:8));
    set(h1, {'color'}, {'b'; 'k'; 'r'});
    hold on;
    h2 = plot(tarray_pst, M(:, 6:8));
    set(h2, {'color'}, {'b'; 'k'; 'r'});
    set(h2, {'linestyle'}, {'--'; '--'; '--'});
    plot(tarray_pst, ghi_clearsky, 'm');
    subplot(4, 1, 4);
    kcs_fix = M_fixed(:, 6:8)./repmat(ghi_clearsky, 1, 3);
    hist(kcs_fix(:), 50);
end

end
function lab_IBM()
% baseline_ibm(2019, 4, 30);
convolve_ibm_new(2019, 4, 30);

% convolve_ibm_old(4);
% CAISO_10_sites();
% explore_correlation(4);
end

function explore_correlation(m)
if m == 4
    dir_work = 'C:\Users\bxl180002\git\SF2\IBM\April\ghi_actual';
elseif m == 5
    dir_work = 'C:\Users\bxl180002\git\SF2\IBM\May\ghi_actual';
end
dir_home = pwd;

% IBMsitenames   = {'MNCC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz'};
% SiteLatitude  = [34.31,   34.12,    37.41,   35.53,   35.38];
% SiteLongitude = [-117.5,-117.94,  -119.74, -118.63, -120.18];

IBMsitenames  = {'CA_Topaz','RSAC1', 'RLKC1', 'SBVC1', 'KNNC1', 'MIAC1', 'MNCC1', 'STFC1', 'DEMC1',  'COWC1'};
SiteLatitude  = [35.38,       38.47,   40.25,   34.45,   40.71,   37.41,   34.31,   34.12,   35.53,   39.12];
SiteLongitude = [-120.18,   -122.71, -123.31, -119.70, -123.92, -119.74, -117.50, -117.94, -118.63, -123.07];

dir_home = pwd;
ghi_actual = [];
ghi_clearsky = [];

for k = 1:length(IBMsitenames)
    ibm_site = IBMsitenames{k};
    csvname_read  = strcat('IBM_processed_', ibm_site, '.hourly.csv');
    
    cd(dir_work);
    M = csvread(csvname_read, 1, 0);
    cd(dir_home);
    
    Location = pvl_makelocationstruct(SiteLatitude(k),SiteLongitude(k)); %Altitude is optional
    Time.UTCOffset(1:size(M,1),1) = zeros(size(M,1), 1); % Because we use UTC time, so utc offset is zero
    Time.year(1:size(M,1),1)   = M(:, 1);
    Time.month(1:size(M,1),1)  = M(:, 2);
    Time.day(1:size(M,1),1)    = M(:, 3);
    Time.hour(1:size(M,1),1)   = M(:, 4);
    Time.minute(1:size(M,1),1) = M(:, 5);
    Time.second(1:size(M,1),1) = zeros(size(M,1), 1);
    [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
    
    ghi_clearsky(:, k) = pvl_clearsky_haurwitz(90-AppSunEl); % Clear-sky GHI
    ghi_actual(:, k) = M(:, end);
end

tarray = datetime(M(:, 1), M(:, 2), M(:, 3), M(:, 4), M(:, 5), M(:, 6), 'TimeZone', 'UTC');
tarray.TimeZone = '-08:00'; % Let's just use PST, so that 12:00 is noon
tarray_pst = datetime(tarray.Year, tarray.Month, tarray.Day, tarray.Hour, tarray.Minute, tarray.Second);

kcs = ghi_actual./ghi_clearsky;
kcs_midday = kcs((tarray_pst.Hour>=7) & (tarray_pst.Hour<17), :);
kcs_other  = kcs((tarray_pst.Hour<7) | (tarray_pst.Hour>=17), :);

figure();
subplot(2, 1, 1);
hist(kcs_midday(:), 50);
title('7am to 5pm (PST)');
xlabel('Clear-sky index');
subplot(2, 1, 2);
hist(kcs_other(:), 50);
title('Other time');
xlabel('Clear-sky index');

dist_gc = nan(numel(IBMsitenames), numel(IBMsitenames));
for i = 1: size(kcs_midday, 2) - 1
    lat1 = SiteLatitude(i);
    lon1 = SiteLongitude(i);
    dist_gc(i, i) = 0;
    for j = i+1: size(kcs_midday, 2)
        lat2 = SiteLatitude(j);
        lon2 = SiteLongitude(j);
        d = DISTANCE(lat1/180*pi, lon1/180*pi, lat2/180*pi, lon2/180*pi);
        dist_gc(i, j) = d;
        dist_gc(j, i) = d;
    end
end
dist_gc(end) = 0;

corr_site = corrcoef(kcs_midday);
figure();
scatter(dist_gc(:), corr_site(:));
xlabel('Distance (km)');
ylabel('Linear correlation coefficient');

% Plot out scatter+histograms between any two pairs of sites
% for i = 1: size(kcs_midday, 2) - 1
%     lat1 = SiteLatitude(i);
%     lon1 = SiteLongitude(i);
%     for j = i+1: size(kcs_midday, 2)
%         figure();
%         lat2 = SiteLatitude(i);
%         lon2 = SiteLongitude(i);
%         scatterhist(kcs_midday(:, i), kcs_midday(:, j));
%         d = distance('gc', lat1, lon1, lat2, lon2);
%     end
% end

end

function baseline_ibm(YYYY, MM, DD)
% Follow CAISO's implementation using IBM data

TARGETDAY = datetime(YYYY, MM, DD);
nhistory = 29;

siteforconv = {'gen55', 'gen56', 'gen59', 'gen60', 'gen61', 'gen58', 'gen62', 'gen64'}; % These 8 gens have non-zero capacities
capacity = nan(length(siteforconv), 2);
cell_data5 = cell(length(siteforconv), 1);
cell_data15 = cell(length(siteforconv), 1);
dirhome = pwd;
dirwork5 ='C:\Users\bxl180002\git\SF2\IBM\April.more_quantiles.5min\power_frcst';
dirwork15 = 'C:\Users\bxl180002\git\SF2\IBM\April\power_frcst';

for i = 1: length(siteforconv)
    % Read 5-min data
    cd(dirwork5);
    csvname = strcat('frcst_', siteforconv{i}, '.csv');
    T = readtable(csvname);
    capacity(i, 1) = max(T.pwr_p095);
    cell_data5{i} = T;
    
    % Read 15-min data
    cd(dirwork15);
    csvname = strcat('frcst_', siteforconv{i}, '.csv');
    T = readtable(csvname);
    capacity(i, 2) = max(T.pwr_p095);
    cell_data15{i} = T;
end
cd(dirhome);
capacity = max(capacity, 2);

T5 = cell_data5{1};
tarray5_utc = datetime(T5.Year, T5.Month, T5.Day, T5.Hour, T5.Minute, zeros(size(T5, 1), 1), 'TimeZone', 'utc');
tarray5_local = tarray5_utc;
tarray5_local.TimeZone = 'America/Los_Angeles';

T15 = cell_data15{1};
tarray15_utc = datetime(T15.Year, T15.Month, T15.Day, T15.Hour, T15.Minute, zeros(size(T15, 1), 1), 'TimeZone', 'utc');
tarray15_local = tarray15_utc;
tarray15_local.TimeZone = 'America/Los_Angeles';

pmean_15 = zeros(length(tarray15_local), 1);
for i = 1: length(cell_data15)
    T = cell_data15{i};
    pmean_15 = pmean_15 + T{:, 'pwr_mean'};
end

pmean_5 = zeros(length(tarray5_local), 1);
for i = 1: length(cell_data5)
    T = cell_data5{i};
    pmean_5 = pmean_5 + T{:, 'pwr_mean'};
end

% Collect history information
pmean_history_5  = nan(12*24, nhistory);
pmean_history_15 = nan(4*24, nhistory);
for i = 1:nhistory
    thisday = TARGETDAY - days(i);
    tarray15_select = tarray15_local((tarray15_local.Year==thisday.Year)&(tarray15_local.Month==thisday.Month)&(tarray15_local.Day==thisday.Day));
%     iend_15 = find(tarray15_local==tarray15_select(end))+1; % This is assume each 15-min interval ended at 15 includes the 5-min intervals ending at 0, 5, and 15
%     istart_15 = find(tarray15_local==tarray15_select(1))+1; % This is assume each 15-min interval ended at 15 includes the 5-min intervals ending at 0, 5, and 15
    iend_15 = find(tarray15_local==tarray15_select(end)); 
    istart_15 = find(tarray15_local==tarray15_select(1));
    tarray5_select = tarray5_local((tarray5_local.Year==thisday.Year)&(tarray5_local.Month==thisday.Month)&(tarray5_local.Day==thisday.Day));
%     iend_5 = find(tarray5_local==tarray5_select(end))+1; % This is assume each 15-min interval ended at 15 includes the 5-min intervals ending at 0, 5, and 15
%     istart_5 = find(tarray5_local==tarray5_select(1))+1; % This is assume each 15-min interval ended at 15 includes the 5-min intervals ending at 0, 5, and 15
    iend_5 = find(tarray5_local==tarray5_select(end));
    istart_5 = find(tarray5_local==tarray5_select(1));

    pmean_history_5(:, nhistory+1-i) = pmean_5(istart_5: iend_5);
    pmean_history_15(:, nhistory+1-i) = pmean_15(istart_15: iend_15); 
end

pmean_delta5 = pmean_history_5 - reshape(repmat(pmean_history_15(:), 1, 3)', numel(pmean_history_15)*3/nhistory, nhistory);
pmean_delta15_max = reshape(max(reshape(pmean_delta5(:), 3, numel(pmean_delta5)/3), [], 1), 4*24, nhistory);
pmean_delta15_min = reshape(min(reshape(pmean_delta5(:), 3, numel(pmean_delta5)/3), [], 1), 4*24, nhistory);

fru_determ = nan(24, 1);
frd_determ = nan(24, 1);
for i = 1: 24
    rtpd_col = (i-1)*4+1: i*4;
    samples_fru = pmean_delta15_max(rtpd_col, :);
    [f,x] = ecdf(samples_fru(:));
    fru_determ(i) = interp1(f, x, 0.95);

    samples_frd = pmean_delta15_min(rtpd_col, :);
    [f,x] = ecdf(samples_frd(:));
    frd_determ(i) = interp1(f, x, 0.05);
end

end

function convolve_ibm_new(YYYY, MM, DD)
% We use MATLAB's conv function, not my version to complete the workflow,
% so no correlation is considered, but let's just start from the easy one.
siteforconv = {'gen55', 'gen56', 'gen59', 'gen60', 'gen61', 'gen58', 'gen62', 'gen64'}; % These 8 gens have non-zero capacities
capacity = nan(length(siteforconv), 2);
cell_data5 = cell(length(siteforconv), 1);
cell_data15 = cell(length(siteforconv), 1);
dirhome = pwd;
dirwork5 ='C:\Users\bxl180002\git\SF2\IBM\April.more_quantiles.5min\power_frcst';
dirwork15 = 'C:\Users\bxl180002\git\SF2\IBM\April\power_frcst';

for i = 1: length(siteforconv)
    % Read 5-min data
    cd(dirwork5);
    csvname = strcat('frcst_', siteforconv{i}, '.csv');
    T = readtable(csvname);
    capacity(i, 1) = max(T.pwr_p095);
    cell_data5{i} = T;
    
    % Read 15-min data
    cd(dirwork15);
    csvname = strcat('frcst_', siteforconv{i}, '.csv');
    T = readtable(csvname);
    capacity(i, 2) = max(T.pwr_p095);
    cell_data15{i} = T;
end
cd(dirhome);
capacity = max(capacity, [], 2);

T5 = cell_data5{1};
tarray5_utc = datetime(T5.Year, T5.Month, T5.Day, T5.Hour, T5.Minute, zeros(size(T5, 1), 1), 'TimeZone', 'utc');
tarray5_local = tarray5_utc;
tarray5_local.TimeZone = 'America/Los_Angeles';

T15 = cell_data15{1};
tarray15_utc = datetime(T15.Year, T15.Month, T15.Day, T15.Hour, T15.Minute, zeros(size(T15, 1), 1), 'TimeZone', 'utc');
tarray15_local = tarray15_utc;
tarray15_local.TimeZone = 'America/Los_Angeles';

actual_binwidth = 10; % MW

% Since IBM's timestamp indicates the 5/15 min period overwhich the GHI is
% averaged, the last timestamp of the three 5-min periods should equal to
% the timestamp of the corresponding 15-min period. Therefore, the
% timestamp of the last period of both resolution should be 00:00 of the
% next day.
tarray15_select = tarray15_local((tarray15_local.Year==YYYY)&(tarray15_local.Month==MM)&(tarray15_local.Day==DD));
% iend_15 = find(tarray15_local==tarray15_select(end))+1; % This is assume each 15-min interval ended at 15 includes the 5-min intervals ending at 0, 5, and 15
% istart_15 = find(tarray15_local==tarray15_select(1))+1; % This is assume each 15-min interval ended at 15 includes the 5-min intervals ending at 0, 5, and 15
iend_15 = find(tarray15_local==tarray15_select(end));
istart_15 = find(tarray15_local==tarray15_select(1));

tarray5_select = tarray5_local((tarray5_local.Year==YYYY)&(tarray5_local.Month==MM)&(tarray5_local.Day==DD));
% iend_5 = find(tarray5_local==tarray5_select(end))+1; % This is assume each 15-min interval ended at 15 includes the 5-min intervals ending at 0, 5, and 15
% istart_5 = find(tarray5_local==tarray5_select(1))+1; % This is assume each 15-min interval ended at 15 includes the 5-min intervals ending at 0, 5, and 15
iend_5 = find(tarray5_local==tarray5_select(end)); 
istart_5 = find(tarray5_local==tarray5_select(1)); 

rowindex15 = istart_15: iend_15;
rowindex5  = istart_5: iend_5;

% rowindex15 = 33:74; % April 1
% rowindex5  = 95:220; % April 1
% rowindex15 = 2907:2960; % May 31
% rowindex5  = 8717:8878; % May 31
% rowindex15 = 2814:2855; % April 30
% rowindex5  = 8438:8563; % April 30


% Result container for convoltuion
cell_conv_bincenter_15 = cell(length(rowindex15), 1);
cell_conv_binheight_15 = cell(length(rowindex15), 1);
cell_conv_bincenter_5 = cell(length(rowindex5), 1);
cell_conv_binheight_5 = cell(length(rowindex5), 1);
cell_conv_bincenter_5m15 = cell(length(rowindex5), 1);
cell_conv_binheight_5m15 = cell(length(rowindex5), 1);
cell_conv_binedge_5m15 = cell(length(rowindex5), 1);
cell_conv_edgecdf_5m15 = cell(length(rowindex5), 1);

% Result container for FRP requirements
fru_p95 = nan(length(rowindex5), 1);
frd_p05 = nan(length(rowindex5), 1);

for ik = 1:length(rowindex15)
    
    % Convolve 15-min data, step 1: Discretization
    k = rowindex15(ik);
    cell_bincenter15 = cell(length(cell_data15), 1);
    cell_binheight15 = cell(length(cell_data15), 1);
    for j = 1: length(cell_data15)
        T = cell_data15{j};
        cap = capacity(j);
        normalized_binwidth = 0.02;
        p = [0.05; 0.95;];
        xp = [T{k,'pwr_p005'}; T{k,'pwr_p095'}]./cap;
        x_bar = T{k,'pwr_mean'}./cap;
        [normalized_h_bin, normalized_binedge, flag_success, MU, SIGMA] = discretize_dist_logitnormal1(p, xp, x_bar, normalized_binwidth);
        if flag_success
            binedges = [0: floor(cap/actual_binwidth)].*actual_binwidth;
            cdf_edges = cdf_logit_normal(binedges./binedges(end), MU, SIGMA);
            bincenter = (binedges(1: end-1) + binedges(2: end))/2;
            binheight = (cdf_edges(2: end) - cdf_edges(1: end-1))./actual_binwidth;
            cell_bincenter15{j} = bincenter;
            cell_binheight15{j} = binheight;
        else
            cell_bincenter15{j} = [];
            cell_binheight15{j} = [];
        end
    end
    
    % Convolve 15-min data, step 2: Aggregation
    bincenter_conv15 = [];
    binheight_conv15 = [];
    for j = 1: length(cell_binheight15)
        if isempty(cell_binheight15{j})
            continue;
        end
        if isempty(binheight_conv15)
            bincenter_conv15 = cell_bincenter15{j};
            binheight_conv15 = cell_binheight15{j};
        else
            % Note we are doing discrete convolution so each element should
            % be the frequency of that discrete number (bin center), so we
            % should multiply the height (pdf) with the width to get the
            % area (frequency) in the bin.
            binheight_conv15 = conv(binheight_conv15.*actual_binwidth, cell_binheight15{j}.*actual_binwidth)./actual_binwidth;
            bincenter_conv15 = (bincenter_conv15(1) + cell_bincenter15{j}(1)):actual_binwidth:(bincenter_conv15(end) + cell_bincenter15{j}(end));
        end
    end
    cell_conv_bincenter_15{ik} = bincenter_conv15;
    cell_conv_binheight_15{ik} = binheight_conv15;

    % Now convolve 5-min data, step 1: Discretization
    for r = 1: 3 % Residual of mod(5-min index, 15-min index)
        i = rowindex5(3*(ik-1)+r);
        cell_bincenter5 = cell(length(cell_data5), 1);
        cell_binheight5 = cell(length(cell_data5), 1);
        for j = 1: length(cell_data5)
            T = cell_data5{j};
            cap = capacity(j);
            normalized_binwidth = 0.02;
            p = [0.05;0.25;0.50;0.75;0.95;];
            xp = [T{i,'pwr_p005'}; T{i,'pwr_p025'}; T{i,'pwr_p050'}; T{i,'pwr_p075'}; T{i,'pwr_p095'}]./cap;
            x_bar = T{i,'pwr_mean'}./cap;
            [normalized_h_bin, normalized_binedge, flag_success, MU, SIGMA] = discretize_dist_logitnormal1(p, xp, x_bar, normalized_binwidth);
            if flag_success
                binedges = [0: floor(cap/actual_binwidth)].*actual_binwidth;
                cdf_edges = cdf_logit_normal(binedges./binedges(end), MU, SIGMA);
                bincenter = (binedges(1: end-1) + binedges(2: end))/2;
                binheight = (cdf_edges(2: end) - cdf_edges(1: end-1))./actual_binwidth;
                cell_bincenter5{j} = bincenter;
                cell_binheight5{j} = binheight;
            else
                cell_bincenter5{j} = [];
                cell_binheight5{j} = [];
            end
        end
        
        % Convolve 5-min data, step 2: Aggregation
        bincenter_conv5 = [];
        binheight_conv5 = [];
        for j = 1: length(cell_binheight5)
            if isempty(cell_binheight5{j})
                continue;
            end
            if isempty(binheight_conv5)
                bincenter_conv5 = cell_bincenter5{j};
                binheight_conv5 = cell_binheight5{j};
            else
                % Note we are doing discrete convolution so each element should
                % be the frequency of that discrete number (bin center), so we
                % should multiply the height (pdf) with the width to get the
                % area (frequency) in the bin.
                binheight_conv5 = conv(binheight_conv5.*actual_binwidth, cell_binheight5{j}.*actual_binwidth)./actual_binwidth;
                bincenter_conv5 = (bincenter_conv5(1) + cell_bincenter5{j}(1)):actual_binwidth:(bincenter_conv5(end) + cell_bincenter5{j}(end));
            end
        end
        cell_conv_bincenter_5{3*(ik-1)+r} = bincenter_conv5;
        cell_conv_binheight_5{3*(ik-1)+r} = binheight_conv5;
    end
    
    % Now, convolution of (5-min PDF) - (15-min PDF)
    for r = 1:3
        f5_center = cell_conv_bincenter_5{3*(ik-1)+r}; % Forecast 5-min center
        f5_height = cell_conv_binheight_5{3*(ik-1)+r}; % Forecast 5-min pdf
        f15_center = -fliplr(cell_conv_bincenter_15{ik}); % Forecast 15-min center
        f15_height = fliplr(cell_conv_binheight_15{ik}); % Forecast 15-min pdf
        if isempty(f5_height) || isempty(f15_height)
            % Any one of the 5-min or the 15-min distribution is empty,
            % let's just assume the forecast error is zero?
            fru_p95(3*(ik-1)+r) = 0;
            frd_p05(3*(ik-1)+r) = 0;
            continue;
        end
        cell_conv_binheight_5m15{3*(ik-1)+r} = conv(f5_height.*actual_binwidth, f15_height.*actual_binwidth)./actual_binwidth;
        cell_conv_bincenter_5m15{3*(ik-1)+r} = (f5_center(1) + f15_center(1)):actual_binwidth:(f5_center(end) + f15_center(end));
        
        % Calculate 95th and 5th percentiles
        cell_conv_binedge_5m15{3*(ik-1)+r} = [cell_conv_bincenter_5m15{3*(ik-1)+r}(:) - actual_binwidth/2;cell_conv_bincenter_5m15{3*(ik-1)+r}(end) + actual_binwidth/2];
        cell_conv_edgecdf_5m15{3*(ik-1)+r} = [0; cumsum(cell_conv_binheight_5m15{3*(ik-1)+r}(:)).*actual_binwidth];
        [~, imin] = min(abs(cell_conv_edgecdf_5m15{3*(ik-1)+r} -0.95));
        fru_p95(3*(ik-1)+r) = cell_conv_binedge_5m15{3*(ik-1)+r}(imin);
        [~, imin] = min(abs(cell_conv_edgecdf_5m15{3*(ik-1)+r} -0.05));
        frd_p05(3*(ik-1)+r) = cell_conv_binedge_5m15{3*(ik-1)+r}(imin);
    end
end

% Uncomment this block if you want to see individual distributions.
% for ik = 1:length(rowindex15)
%     figure();
%     subplot(2, 1, 1);
%     k = rowindex15(ik);
%     plot(cell_conv_bincenter_15{ik}, cell_conv_binheight_15{ik}, 'k');
%     hold on;
%     for r = 1: 3
%         plot(cell_conv_bincenter_5{3*(ik-1)+r}, cell_conv_binheight_5{3*(ik-1)+r}, 'r');
%     end
%     xlabel('MW');
%     
%     subplot(2, 1, 2);
%     hold on;
%     for r = 1: 3
%         yyaxis left;
%         plot(cell_conv_bincenter_5m15{3*(ik-1)+r}, cell_conv_binheight_5m15{3*(ik-1)+r}, 'r');
%         yyaxis right;
%         plot(cell_conv_binedge_5m15{3*(ik-1)+r}, cell_conv_edgecdf_5m15{3*(ik-1)+r});
%     end
%     yyaxis left;
%     ylabel('PDF');
%     yyaxis right;
%     ylabel('CDF');
%     ylim([0, 1]);
%     xlabel('MW');
% end

end

function convolve_ibm_old(m, write_flag)
% Note this is used to produce convolution using my own conv_poly function
% back in June, 2019. 
if nargin == 1
    write_flag = false;
end
if m == 4
    dirwork = 'C:\Users\bxl180002\git\SF2\IBM\April\power_frcst.bak'; % IBM's updated forecast
elseif m == 5
    dirwork = 'C:\Users\bxl180002\git\SF2\IBM\May\power_frcst.bak'; % IBM's updated forecast
end
dirhome = pwd;
% siteforconv = {'gen55', 'gen56', 'gen59', 'gen60', 'gen61'}; % This 5 sites correspond to 5 unique IBM sites
siteforconv = {'gen55', 'gen56', 'gen59', 'gen60', 'gen61', 'gen58', 'gen62', 'gen64'}; % These 8 gens have non-zero capacities
ibm_sites   = {'MNCC1', 'STFC1', 'MIAC1', 'DEMC1', 'CA_Topaz'};
capacity_gen =  [232.44, 107.58, 116.05, 140.4, 151.32, ];
capacity_site = [372.84, 134.88, 116.05, 338.52, 151.32];
cap_scaler = capacity_site./capacity_gen;
csv_write = 'power_convolution.csv';

% Read data
capacity = nan(length(siteforconv), 1);
cell_data = cell(length(siteforconv), 1);
cd(dirwork);
for i = 1: length(siteforconv)
    csvname = strcat('power_', siteforconv{i}, '.csv');
    M = csvread(csvname, 1, 0);
%     M(:, 7:9) = M(:, 7:9).*cap_scaler(i);
    M(M(:, 2) ~=m, :) = []; % Remove march data
    tarray = datetime(M(:, 1), M(:, 2), M(:, 3), M(:, 4), M(:, 5), M(:, 6), 'TimeZone', 'America/Los_Angeles');
%     tarray.TimeZone = 'America/Los_Angeles';
%     M(:, 1) = tarray.Year;
%     M(:, 2) = tarray.Month;
%     M(:, 3) = tarray.Day;
%     M(:, 4) = tarray.Hour;
%     M(:, 5) = tarray.Minute;
%     M(:, 6) = tarray.Second;
    cell_data{i} = M;
    capacity(i) = max(M(:, end));
end
cd(dirhome);

nT = size(M, 1);

% Identify valid rows
i_valid   = true(nT, 1);
i_allzero = true(nT, 1);
for i = 1: length(cell_data)
    M = cell_data{i};
    i_valid   = i_valid & (M(:, 7) <= M(:, 8)) & (M(:, 8) <= M(:, 9)); 
    i_allzero = i_allzero & (M(:, 7) <=1E-3) & (M(:, 8) <=1E-3) & (M(:, 9) <=1E-3);
end
t_idx = 1: 1: nT;
% t_idx = t_idx(i_valid);

% Direct sum, if perfectly correlated
M_agg = zeros(size(M(:, 7:9)));
for i = 1: length(cell_data)
    M = cell_data{i};
    i_valid   = i_valid & (M(:, 7) <= M(:, 8)) & (M(:, 8) <= M(:, 9)); 
    i_allzero = i_allzero & (M(:, 7) <=1E-3) & (M(:, 8) <=1E-3) & (M(:, 9) <=1E-3);
    M_agg = M_agg + M(:, 7:9);
end

cell_a = cell(nT, 1);
cell_t = cell(nT, 1);
x_agg_950 = nan(nT, 1);
x_agg_050 = nan(nT, 1);
x_agg_bar = nan(nT, 1);
% Do convolution
% for t = 1: nT
for d = 1: 29 % Day 1 to 29
    x_agg_950_d = nan(24*4, 1);
    x_agg_050_d = nan(24*4, 1);
    x_agg_bar_d = nan(24*4, 1);
%     for t = 31:36 % Test samples
    for t = t_idx(i_valid&(M(:, 3)==d))
    % for t = 2466-28 % This is the row with negative distribution, 28 is the number of rows in March
        n = 0;
        a_convoluted = [];
        t_convoluted = [];
%         fig_distapprox = figure();
%         fig_convbystep = figure();
        x_agg_050_sum = 0;
        x_agg_950_sum = 0;
        x_agg_bar_sum = 0;
        for i = 1:length(cell_data)
            M = cell_data{i};
            cap = capacity(i);
            normalized_binwidth = 0.02;

            % Use logit-normal function to fit the distribution
            x500_normalized = M(t, 8)/cap;
            x050_normalized = M(t, 7)/cap;
            x950_normalized = M(t, 9)/cap;
            x_agg_050_sum = x_agg_050_sum + M(t, 7);
            x_agg_bar_sum = x_agg_bar_sum + M(t, 8);
            x_agg_950_sum = x_agg_950_sum + M(t, 9);
            [normalized_h_bin, normalized_binedge, flag_success] = discretize_dist_logitnormal1([0.05;0.95], [x050_normalized; x950_normalized], x500_normalized, normalized_binwidth);


%             % Visualize the approximation, plot out the given percentiles
%             % and mean (red lines), the logit-normal approximation (curve)
%             % and its percentiles and mean (green lines), and the
%             % discretized approximation (stairs) and its percentiles and
%             % mean (blue lines).
%             figure(fig_distapprox);
%             subplot(3, 3, i)
%             plot_logit(MU, SIGMA);
%             hold on;
%             stairs_conv_poly(fig_distapprox, [normalized_h_bin(:);0], normalized_binedge, 'b');
%             y1 = get(gca, 'ylim');
%             x500 = x500_normalized;
%             x050 = x050_normalized;
%             x950 = x950_normalized;
%             plot([x500 x500],y1 ,'r', 'LineWidth', 1.5)
%             plot([x050 x050],y1 ,'r')
%             plot([x950 x950],y1 ,'r')
%             p = percentile_logit([0.05;0.95], MU, SIGMA);
%             deltax = 0.01;
%             x_mc = 0.01:deltax:0.99;
%             p_mc = pdf_logit_normal(x_mc, MU, SIGMA);
%             x_bar = deltax*sum(x_mc(:).*p_mc(:));
%             plot([x_bar x_bar],y1 ,'g')
%             plot([p(1) p(1)],y1 ,'g')
%             plot([p(2) p(2)],y1 ,'g')
%             edge_cdf = cdf_logit_normal(normalized_binedge, MU, SIGMA);
%             x050_d = interp1(edge_cdf, normalized_binedge, 0.05); % Discretized
%             x950_d = interp1(edge_cdf, normalized_binedge, 0.95);
%             x_bar_d = normalized_h_bin*(normalized_binedge(2:end).^2 - normalized_binedge(1:length(normalized_binedge)-1).^2)'/2;
%             plot([x_bar_d x_bar_d],y1 ,'b')
%             plot([x050_d x050_d],y1 ,'b')
%             plot([x950_d x950_d],y1 ,'b')

            if flag_success
                h_bin = normalized_h_bin./cap;
                binedge = normalized_binedge.*cap;
                a_fitted = [h_bin(:); 0];
                t_fitted = binedge(:);

                if (i == 1) || isempty(t_convoluted)
                    a_convoluted = a_fitted;
                    t_convoluted = t_fitted;
                else
                    [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, a_fitted, t_fitted, normalized_binwidth*cap);
                end
                [a_convoluted, t_convoluted] = reformulate_poly(a_convoluted, t_convoluted);
                n = n + 1;
            else
                continue
            end

            
%             % Visualize convolution step by step
%             figure(fig_convbystep);
%             subplot(3, 3, i);
%             stairs_conv_poly(fig_convbystep, a_convoluted, t_convoluted, 'b');
%             if ~isempty(a_convoluted)
%                 cdf_agg = cdf_poly(a_convoluted, t_convoluted);
%                 if any(diff(cdf_agg)==0)
%                     % This case, we cannot use interp1 function
%                     [~, imin] = min(abs(cdf_agg - 0.05));
%                     x_agg_050_convstep = t_convoluted(imin);
%                     [~, imin] = min(abs(cdf_agg - 0.95));
%                     x_agg_950_convstep = t_convoluted(imin);
%                 else
%                     x_agg_050_convstep = interp1(cdf_agg, t_convoluted, 0.05);
%                     x_agg_950_convstep = interp1(cdf_agg, t_convoluted, 0.95);
%                 end
%                 x_agg_bar_convstep = a_convoluted(1:size(a_convoluted, 1)-1)'*(t_convoluted(2:end).^2 - t_convoluted(1:length(t_convoluted)-1).^2)/2;
%             end
%             hold on;
%             y2 = get(gca, 'ylim');
%             plot([x_agg_bar_convstep x_agg_bar_convstep],y2 ,'b')
%             plot([x_agg_050_convstep x_agg_050_convstep],y2 ,'b')
%             plot([x_agg_950_convstep x_agg_950_convstep],y2 ,'b')
%             plot([x_agg_bar_sum x_agg_bar_sum],y2 ,'r')
%             plot([x_agg_050_sum x_agg_050_sum],y2 ,'r')
%             plot([x_agg_950_sum x_agg_950_sum],y2 ,'r')

        end
        cell_a{t} = a_convoluted;
        cell_t{t} = t_convoluted;
        
        t_thatday = t - (d-1)*24*4;
        if ~isempty(a_convoluted)
            cdf_agg = cdf_poly(a_convoluted, t_convoluted);
            if any(diff(cdf_agg)<=0)
                % This case, we cannot use interp1 function
                [~, imin] = min(abs(cdf_agg - 0.05));
                x_agg_050(t) = t_convoluted(imin);
                x_agg_050_d(t_thatday) = t_convoluted(imin);
                [~, imin] = min(abs(cdf_agg - 0.95));
                x_agg_950(t) = t_convoluted(imin);
                x_agg_950_d(t_thatday) = t_convoluted(imin);
            else
                x_agg_050(t) = interp1(cdf_agg, t_convoluted, 0.05);
                x_agg_050_d(t_thatday) = x_agg_050(t);
                x_agg_950(t) = interp1(cdf_agg, t_convoluted, 0.95);
                x_agg_950_d(t_thatday) = x_agg_950(t);
            end
            x_agg_bar(t) = a_convoluted(1:size(a_convoluted, 1)-1)'*(t_convoluted(2:end).^2 - t_convoluted(1:length(t_convoluted)-1).^2)/2;
            x_agg_bar_d(t_thatday) = x_agg_bar(t);
        end
        if i_allzero(t)
            x_agg_050(t) = 0;
            x_agg_950(t) = 0;
            x_agg_050_d(t_thatday) = 0;
            x_agg_950_d(t_thatday) = 0;
        end
        fprintf('%4g/%g Done.\n', t, nT);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix: Assume all nan is zero, i.e., when low 95 CI and mean are 0, and up
% 95 CI is not zero, we assume it is zero
x_agg_050(isnan(x_agg_050)) = 0;
x_agg_950(isnan(x_agg_950)) = 0;
x_agg_bar(isnan(x_agg_bar)) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d = 1: 29
    % One day's snapshot
    figure();
    x_agg_950_d = x_agg_950((M(:, 2)==m)&(M(:, 3)==d));
    x_agg_050_d = x_agg_050((M(:, 2)==m)&(M(:, 3)==d));
    x_agg_bar_d = x_agg_bar((M(:, 2)==m)&(M(:, 3)==d));

    tarray_d = tarray((M(:, 2)==m)&(M(:, 3)==d)); % T array for the day
    h1 = plot(tarray_d, x_agg_050_d, 'b');
    hold on;
    h2 = plot(tarray_d, x_agg_950_d, 'r');
    h3 = plot(tarray_d, x_agg_bar_d, 'k');
    h4 = plot(tarray_d, M_agg(M(:, 3)==d, :));
    set(h4, {'linestyle'}, {'--'; '--'; '--'});
    set(h4, {'color'}, {'b'; 'k'; 'r'});
    legend([h4; h1; h2; h3], {'SUM 5%'; 'SUM MEAN'; 'SUM 95%'; 'CI 5%'; 'CI 95%'; 'CI MEAN'});
    ylabel('kW');
    title(d);

end

figure();
tarray = datetime(M(:, 1), M(:, 2), M(:,3), M(:, 4), M(:, 5), 0);
h1 = plot(tarray, x_agg_050, 'b');
hold on;
h2 = plot(tarray, x_agg_950, 'r');
h3 = plot(tarray, x_agg_bar, 'k');
legend([h1, h2, h3], {'CI 5%', 'CI 95%', 'CI MEAN'});
ylabel('kW');
% ylim([-10, 800]);
% sum(isnan(x_agg_050))/size(x_agg_050, 1);

if write_flag
    results = [M(:, 1: 6), x_agg_050, x_agg_bar, x_agg_950];
    cHeader = {'Year' 'Month' 'Day' 'Hour' 'Minute' 'Second' '5p' 'Mean' '95p'}; %dummy header
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); % cHeader in text with commas
    fid = fopen(csv_write,'w');
    fprintf(fid,'%s\n',textHeader); % write header to file
    fclose(fid);
    dlmwrite(csv_write, results,'-append');
end

end

function CAISO_10_sites()
% Test IBM's data
dirwork = 'C:\Users\bxl180002\Downloads\RampSolar\IBM_old\sample';
dirhome = pwd;
listcsv = {...
    'sample_35.38_-120.18.csv', ...
    'sample_38.47_-122.71.csv', ...
    'sample_40.25_-123.31.csv', ...
    'sample_34.45_-119.7.csv', ...
    'sample_40.71_-123.92.csv', ...
    'sample_37.41_-119.74.csv', ...
    'sample_34.31_-117.5.csv', ...
    'sample_34.12_-117.94.csv', ...
    'sample_35.53_-118.63.csv', ...
    'sample_39.12_-123.07.csv', ...
};

cap = 250; % kW

% Show fitted distributions
row = 36;
for i = 1:length(listcsv)
    fig = figure();
    csvname = listcsv{i};
    cd(dirwork);
    M = csvread(csvname, 1, 1);
    cd(dirhome);

    % Use logit-normal function to fit the distribution
    hour = mod(M(row, 4) - 8+24, 24);
    x500 = M(row, 6)/cap;
    x025 = M(row, 7)/cap;
    x975 = M(row, 8)/cap;
    if (x500~=0) && (x025~=0)
        param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.025-0.5)) 1]\[log(x500/(1-x500)); log(x025/(1-x025))];
    elseif (x025~=0) && (x975~=0)
        if x975==cap
            x975=x975-eps;
        end
        param = [sqrt(2)*erfinv(2*(0.025-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[log(x025/(1-x025)); log(x975/(1-x975))];
    elseif (x500~=0) && (x975~=0)
        if x975==cap
            x975=x975-eps;
        end
        if x975==cap
            x975=x975-eps;
        end
        param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[log(x500/(1-x500)); log(x975/(1-x975))];
    else
        continue % We need at least two non-zero points.
    end
    SIGMA = param(1);
    MU = param(2);
%     plot_logit(MU, SIGMA);
%     hold on;

    % Prepare for convolution
    binwidth = 10; % kW
    normalized_binedge = 0:binwidth/cap:1; % This is normalized bin edges
    edge_cdf = cdf_logit_normal(normalized_binedge, MU, SIGMA);
    normalized_h_bin = (edge_cdf(2:end) - edge_cdf(1:length(edge_cdf)-1)).*cap/binwidth;
    
    h_bin = normalized_h_bin./250;
    binedge = normalized_binedge.*cap;
    
    % Plot discretized and continuous figures
    stairs_conv_poly(fig, [h_bin(:); 0], binedge, 'k');
    hold on;
    plot(0:binwidth:cap, pdf_logit_normal(normalized_binedge, MU, SIGMA)./cap, 'b', 'LineWidth', 1.5);
    y1=get(gca,'ylim');
    plot([x500 x500].*cap,y1 ,'r', 'LineWidth', 1.5)
%     plot([x025 x025].*cap,y1)
%     plot([x975 x975].*cap,y1)
    patch([x025, x975, x975, x025].*cap, [y1(1), y1(1), y1(2), y1(2)], 'k', 'FaceAlpha', 0.2);
    title(i);
    xlabel('kW');
    set(gca, 'FontSize', 12');
    box on;
end

% See the distribution of aggretated solar power output
rows_selected = 30:76;
% rows_selected = 36;

for row = rows_selected
    fig = figure;
    n = 0; % Number of plants with non-negative output
    for i = 1:length(listcsv)
        csvname = listcsv{i};
        cd(dirwork);
        M = csvread(csvname, 1, 1);
        cd(dirhome);

        % Use logit-normal function to fit the distribution
        hour = mod(M(row, 4) - 8+24, 24);
        x500 = M(row, 6)/cap;
        x025 = M(row, 7)/cap;
        x975 = M(row, 8)/cap;
        if (x500~=0) && (x025~=0)
%             param = [sqrt(2)*erfinv(2*(0.050-0.5)) 1; sqrt(2)*erfinv(2*(0.025-0.5)) 1]\[x500/(1-x500); x025/(1-x025)];
            param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.025-0.5)) 1]\[log(x500/(1-x500)); log(x025/(1-x025))];
        elseif (x025~=0) && (x975~=0)
            if x975==cap
                x975=x975-eps;
            end
            param = [sqrt(2)*erfinv(2*(0.025-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[log(x025/(1-x025)); log(x975/(1-x975))];
%             param = [sqrt(2)*erfinv(2*(0.025-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[x025/(1-x025); x975/(1-x975)];
        elseif (x500~=0) && (x975~=0)
            if x975==cap
                x975=x975-eps;
            end
            if x975==cap
                x975=x975-eps;
            end
            param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[log(x500/(1-x500)); log(x975/(1-x975))];
%             param = [sqrt(2)*erfinv(2*(0.050-0.5)) 1; sqrt(2)*erfinv(2*(0.975-0.5)) 1]\[x500/(1-x500); x975/(1-x975)];
        else
            continue % We need at least two non-zero points.
        end
        n = n + 1;
        SIGMA = param(1);
        MU = param(2);
    %     plot_logit(MU, SIGMA);
    %     hold on;

        % Prepare for convolution
%         binwidth = 10; % kW
%         binedge = 0:binwidth/cap:1; % This is normalized bin edges
%         edge_cdf = cdf_logit_normal(binedge, MU, SIGMA);
%         binedge = binedge.*cap;
%         bin_pdf = (edge_cdf(2:end) - edge_cdf(1:length(edge_cdf)-1))./binwidth;
        
        binwidth = 10; % kW
        normalized_binedge = 0:binwidth/cap:1; % This is normalized bin edges
        edge_cdf = cdf_logit_normal(normalized_binedge, MU, SIGMA);
        normalized_h_bin = (edge_cdf(2:end) - edge_cdf(1:length(edge_cdf)-1)).*cap/binwidth;
        
        h_bin = normalized_h_bin./250;
        binedge = normalized_binedge.*cap;



        if i == 1
            a_convoluted = [h_bin(:); 0];
            t_convoluted = binedge;
        else
            [a_convoluted, t_convoluted] = conv_poly(a_convoluted, t_convoluted, [h_bin(:); 0], binedge, binwidth);
        end
        ax = subplot(3, 4, i);
        plot_conv_poly(fig, a_convoluted, t_convoluted, 'b');
        title(n);
    end
end

end

function stairs_conv_poly(fig, a, t, c)
figure(fig); hold on;
x_l_all = nan(size(a, 1)-1, 1);
y_l_all = nan(size(a, 1)-1, 1);
x_r_all = nan(size(a, 1)-1, 1);
y_r_all = nan(size(a, 1)-1, 1);
for i = 1:(size(a, 1)-1)
    x_l = t(i);
    x_r = t(i+1);
    x = linspace(x_l, x_r, 10);
    y = polyval(fliplr(a(i, :)), x);
%     h = stairs(x, y, 'color', c);
    x_l_all(i) = x_l;
    x_r_all(i) = x_r;
    y_l_all(i) = y(1);
    y_r_all(i) = y(end);
end
stairs(x_l_all, y_l_all, 'k')

end

function [normalized_h_bin, normalized_binedge, flag_success] = discretize_dist_logitnormal(x050, x500, x950, normalized_binwidth)
% Fit IBM's forecasted percentiles to distributions
% All inputs should be normalized, i.e., between 0 and 1
flag_lack_data = false;
if (x050~=0) && (x950~=0)
    if x950==1
        x950=x950-eps;
    end
    param = [sqrt(2)*erfinv(2*(0.050-0.5)) 1; sqrt(2)*erfinv(2*(0.950-0.5)) 1]\[log(x050/(1-x050)); log(x950/(1-x950))];
elseif (x500~=0) && (x050~=0)
    param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.050-0.5)) 1]\[log(x500/(1-x500)); log(x050/(1-x050))];
elseif (x500~=0) && (x950~=0)
    if x950==1
        x950=x950-eps;
    end
    param = [sqrt(2)*erfinv(2*(0.500-0.5)) 1; sqrt(2)*erfinv(2*(0.950-0.5)) 1]\[log(x500/(1-x500)); log(x950/(1-x950))];
else
    flag_lack_data = true; % We need at least two non-zero points, otherwise we consider it as zero
end

if ~flag_lack_data
    SIGMA = param(1);
    MU = param(2);
    %     plot_logit(MU, SIGMA);
    %     hold on;

    % Prepare for convolution        
    normalized_binedge = 0:normalized_binwidth:1; % This is normalized bin edges
    edge_cdf = cdf_logit_normal(normalized_binedge, MU, SIGMA);
    normalized_h_bin = (edge_cdf(2:end) - edge_cdf(1:length(edge_cdf)-1)).*1/normalized_binwidth;

    flag_success = true;
else
    normalized_binedge = nan;
    normalized_h_bin = nan;
    flag_success = false;
end

end

% function [normalized_h_bin, normalized_binedge, flag_success, MU, SIGMA] = discretize_dist_logitnormal1(x050, x500, x950, normalized_binwidth)
function [normalized_h_bin, normalized_binedge, flag_success, MU, SIGMA] = discretize_dist_logitnormal1(p, xp, x_bar, normalized_binwidth)
% Fit IBM's forecasted percentiles to distributions
% All inputs should be normalized, i.e., between 0 and 1
% p is the percentiles in ascending order, i.e., 0.05, 0.50, 0.95
% xp is the values associated with the percentiles in p
% x_bar is the mean
flag_lack_data = false;
% if (x050~=0) && (x950~=0)
%     if x950==1
%         x950=x950-eps;
%     end
% elseif (x500~=0) && (x950~=0)
%     if x050==0
%         x050=x050-eps;
%     end
if sum(xp~=0) < 2
    flag_lack_data = true; % We need at least two non-zero points, otherwise we consider it as zero
end

if ~flag_lack_data
    % Find the best mu and sigma by absolute difference of 5, 95 percentiles
    % and mean
    MU_all    = -3:0.05:3;
    SIGMA_all = 0.1: 0.05: 3;
    absdiff = nan(numel(SIGMA_all), numel(MU_all));
    for i = 1: length(SIGMA_all)
        for j = 1: length(MU_all)
            MU = MU_all(j);
            SIGMA = SIGMA_all(i);
            xp_fitted = percentile_logit(p, MU, SIGMA);
            deltax = 0.01;
            x_mc = 0.01:deltax:0.99;
            p_mc = pdf_logit_normal(x_mc, MU, SIGMA);
            x_bar_fitted = deltax*sum(x_mc(:).*p_mc(:));
            absdiff(i, j) = sum(abs([xp_fitted(:); x_bar_fitted] - [xp; x_bar]));
        end
    end
    [~, imin] = min(absdiff(:));
    [subi, subj] = ind2sub(size(absdiff), imin);
    SIGMA = SIGMA_all(subi);
    MU    = MU_all(subj);
    

    % Prepare for convolution        
    normalized_binedge = 0:normalized_binwidth:1; % This is normalized bin edges
    edge_cdf = cdf_logit_normal(normalized_binedge, MU, SIGMA);
    normalized_h_bin = (edge_cdf(2:end) - edge_cdf(1:length(edge_cdf)-1)).*1/normalized_binwidth;

    flag_success = true;
    
else
    normalized_binedge = nan;
    normalized_h_bin = nan;
    flag_success = false;
    SIGMA = nan;
    MU    = nan;
end



end



function y = cdf_logit_normal(x, MU, SIGMA)

y = 0.5.*( 1+erf((logit(x)-MU)./(sqrt(2)*SIGMA)) );
end

function y = pdf_logit_normal(x, MU, SIGMA)
logitx = logit(x);
y = 1/(SIGMA*sqrt(2*pi))*1./(x.*(1-x)).*exp(-(logitx-MU).^2./(2*SIGMA^2));
end

function x = percentile_logit(ALPHA, MU, SIGMA)
% Return the percentiles specified by ALPHA, MU and SIGMA are scalers
tmp = MU + sqrt(2)*SIGMA*erfinv(2.*ALPHA - 1);
x = exp(tmp)./(1+exp(tmp));

end

function plot_logit(MU, SIGMA)
x = 0.01:0.01:0.99;
% logitx = logit(x);
% y = 1/(SIGMA*sqrt(2*pi))*1./(x.*(1-x)).*exp(-(logitx-MU).^2./(2*SIGMA^2));
y = pdf_logit_normal(x, MU, SIGMA);
plot(x, y);
end

function y = logit(x)
y = log(x./(1-x));
end

function [d] = DISTANCE(lat1, lon1, lat2, lon2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return the distance between two point(lat1, lon1) and (lat2, lon2).
% Input: in rad; Output: in km.
% Reference: http://www.movable-type.co.uk/scripts/latlong.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It is worth noting that the transmission cable length is not necessarily
% equal to the distance between two points, there might be a coefficient to
% convert between straight-line distance and transmission distance. A
% literature review will be necessary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 6361; % Earth radius
dlat = lat2 - lat1;
dlon = lon2 - lon1;
a = sin(dlat./2).^2 + sin(dlon./2).^2.*cos(lat1).*cos(lat2);
d = R.*2.*atan2(sqrt(a), sqrt(1 - a));
% d = R.*c;

% Spherical Law of Cosines
% d = R*acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlat));
end

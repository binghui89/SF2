function ACPower_mean = reference_curve(lat, lon, tarray, deltat)
ACPower_mean = use_pvlib(lat, lon, tarray, deltat);
end

function ACPower_mean = use_piecewise(lat, lon, tarray, deltat)

end

function ACPower_mean = use_pvlib(lat, lon, tarray, deltat)

Location = pvl_makelocationstruct(lat, lon);
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
[~, ~, AppSunEl_formean, ~] = pvl_ephemeris(Time_formean,Location);
ghi_cs_formean = pvl_clearsky_haurwitz(90-AppSunEl_formean);
% ghi_cs_mean = mean(reshape(ghi_cs_formean, numel(ghi_cs_formean)/numel(tarray), numel(tarray)), 1)';


sitenames={'gen55','gen56','gen57','gen58','gen59','gen60','gen61','gen62','gen63','gen64'};
SiteLatitude=[34.31,34.12,34.12,34.12,37.41,35.53,35.38,34.31,34.31,35.53];
SiteLongitude=[-117.5,-117.94, -117.94, -117.94,-119.74, -118.63, -120.18,-117.5,-117.5,-118.63];
IBMsitenames={'MNCC1','STFC1','STFC1','STFC1','MIAC1','DEMC1','Topaz','MNCC1','MNCC1', 'DEMC1'};
modulepv=134; %Topaz uses first solar panels (FS272 is probably the oldest they have)
inverterid=759; 
% Site_tilt=[0,0,0,25,0,0,25,22.5,0,0];
Site_tilt=[0,0,0,0,0,0,0,0,0,0]; % Single tracking
modules_Series=11;
modules_parallel=[2360,1870,0,2002,2096,2512,2319,2381,0,2353];
ninvert=[149,97,0,23,82,90,97,90,0,127];


k = 1; % 1 to 10: 3 and 9 are 0

%Other parameters
SF=0.98;
%Weather
PresPa=101325;
WIND=0;
dryT=10;
Albedo = 0.2;

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
% ReGHI = [utc_year, utc_month, utc_day, utc_hour, utc_minute, utc_second, ghi];
% Time.UTCOffset(1:size(ReGHI,1),1) = zeros(size(ReGHI,1), 1); % Because we use UTC time, so utc offset is zero
% Time.year(1:size(ReGHI,1),1)   = utc_year;
% Time.month(1:size(ReGHI,1),1)  = utc_month;
% Time.day(1:size(ReGHI,1),1)    = utc_day;
% Time.hour(1:size(ReGHI,1),1)   = utc_hour;
% Time.minute(1:size(ReGHI,1),1) = utc_minute;
% Time.second(1:size(ReGHI,1),1) = utc_second;
% ACPower(1:size(Time_formean,1),1:6) = ReGHI(:,1:6);

%used for both forecast and actual
dayofyear = pvl_date2doy(Time_formean.year, Time_formean.month, Time_formean.day);
[SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time_formean,Location);

if Site_tilt(k)==0    
    [TrkrTheta, AOI, SurfTilt, SurfAz] = pvl_singleaxis(90-AppSunEl, SunAz, Location.latitude, 0, 180, 45);
    Array.Tilt=SurfTilt;
    Array.Tilt(find(isnan(Array.Tilt)))=0;
else
    AOI = pvl_getaoi(Array.Tilt, Array.Azimuth, 90-AppSunEl, SunAz);
end
Wspd=repmat(WIND,size(Time_formean,1), 1);
Drybulb=repmat(dryT,size(Time_formean,1), 1);
AMa = pvl_absoluteairmass(pvl_relativeairmass(90-AppSunEl),PresPa);
F1 = max(0,polyval(ModuleParameters.a,AMa)); %Spectral loss function
F2 = max(0,polyval(ModuleParameters.b,AOI)); % Angle of incidence loss function

EdiffGround = pvl_grounddiffuse(Array.Tilt, ghi_cs_formean, Albedo);
DNI_model = pvl_disc(ghi_cs_formean,90-SunEl, dayofyear,PresPa);
DHI_model = ghi_cs_formean - cosd(90-SunEl).*DNI_model;
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
% temp=find(ghi_cs_formean~=10000);
% ACPower(1:size(Time_formean,1)) =10000;
% ACPower_mean = pvl_snlinverter(Inverter, mSAPMResults.Vmp(temp)*Array.Ms, mSAPMResults.Pmp(temp)*Array.Ms*Array.Mp)*ninvert(k)/1000000;

ACPower_formean = pvl_snlinverter(Inverter, mSAPMResults.Vmp*Array.Ms, mSAPMResults.Pmp*Array.Ms*Array.Mp)*ninvert(k)/1000000;
ACPower_formean(ACPower_formean<0)=0;
ACPower_mean = mean(reshape(ACPower_formean, numel(ACPower_formean)/numel(tarray), numel(tarray)), 1)';

% DCPower_formean = (mSAPMResults.Vmp*Array.Ms).*(mSAPMResults.Pmp*Array.Ms*Array.Mp)*ninvert(k)/1000000;
% DCPower_formean(DCPower_formean<0)=0;
% ACPower_mean = mean(reshape(DCPower_formean, numel(DCPower_formean)/numel(tarray), numel(tarray)), 1)';


end
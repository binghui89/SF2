%mape for GHI
clear all
sitenames={'gen55','gen56','gen57','gen58','gen59','gen60','gen61','gen62','gen63','gen64'};
SiteLatitude=[34.31,34.12,34.12,34.12,37.41,35.53,35.38,34.31,34.31,35.53];
SiteLongitude=[-117.5,-177.94, -177.94, -177.94,-119.74, -118.63, -120.18,-117.5,-117.5,-118.63];
IBMsitenames={'MNCC1','STFC1','STFC1','STFC1','MIAC1','DEMC1','Topaz','MNCC1','MNCC1', 'DEMC1'};
modulepv=134; %Topaz uses first solar panels (FS272 is probably the oldest they have)
inverterid=759; 
Site_tilt=[0,0,0,25,0,0,25,22.5,0,0];
modules_Series=11;
modules_parallel=[2360,1870,0,2002,2096,2512,2319,2381,0,2353];
ninvert=[149,97,0,23,82,90,97,90,0,127];

%Other parameters
SF=0.98;
%Weather
PresPa=101325;
WIND=0;
dryT=10;
Albedo = 0.2;


for k = 1:size(sitenames,2)
    clearvars -except sitenames SiteLatitude SiteLongitude  IBMsitenames k modulepv inverterid Site_tilt modules_Series modules_parallel ninvert SF PresPa WIND dryT Albedo
    Location = pvl_makelocationstruct(SiteLatitude(k),SiteLongitude(k)); %Altitude is optional
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

    %--- Read GHI data from excel file
    ReGHI=xlsread(['March2019\',char(IBMsitenames(k)),'_Mar2019_5minGHIv4utc.xlsx'],'Reconstructed_GHI');
    DN = datenum(2019, 3,1):5/(24*60):datenum(2019, 3, 31, 23, 59, 59);
    Time5temp=char(tzoffset(datetime(cellstr(datestr(DN)), 'Timezone','America/Los_Angeles')));
    clear temppppp
    temppppp=find(not(isnan(tzoffset(datetime(cellstr(datestr(DN)), 'Timezone','America/Los_Angeles'))))); 
    Time.UTCOffset(1:size(ReGHI,1),1)=str2num(Time5temp(:,1:2));
    Time.year(1:size(ReGHI,1),1)=2019;
    Time.month(1:size(ReGHI,1),1)=3;
    Time.day(1:size(ReGHI,1),1)=ReGHI(:,2);
    Time.hour(1:size(ReGHI,1),1)=ReGHI(:,3);
    Time.minute(1:size(ReGHI,1),1)=ReGHI(:,4);
    Time.second(1:size(ReGHI,1),1)=0;
    ACPower(1:size(ReGHI,1),1:4)=ReGHI(:,1:4);
    
    %used for both forecast and actual
    dayofyear = pvl_date2doy(Time.year, Time.month, Time.day);
    [SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
    
    if Site_tilt(k)==0    
        [TrkrTheta, AOI, SurfTilt, SurfAz] = pvl_singleaxis(90-AppSunEl, SunAz, Location.latitude, 0, 180, 45);
        Array.Tilt=SurfTilt;
        Array.Tilt(find(isnan(Array.Tilt)))=0;
    else
        AOI = pvl_getaoi(Array.Tilt, Array.Azimuth, 90-AppSunEl, SunAz);
    end
    Wspd=repmat(WIND,size(ReGHI,1));
    Drybulb=repmat(dryT,size(ReGHI,1));
    AMa = pvl_absoluteairmass(pvl_relativeairmass(90-AppSunEl),PresPa);
    F1 = max(0,polyval(ModuleParameters.a,AMa)); %Spectral loss function
    F2 = max(0,polyval(ModuleParameters.b,AOI)); % Angle of incidence loss function
    
    GHI_forecast=ReGHI(:,8);
    GHI_actual=ReGHI(:,11);
    GHI_upperforecast=ReGHI(:,12);
    GHI_lowerforecast=ReGHI(:,13);
    
    %ACPower has the estimate based on actual measurement
    EdiffGround = pvl_grounddiffuse(Array.Tilt,GHI_actual, Albedo);
    DNI_model = pvl_disc(GHI_actual,90-SunEl, dayofyear,PresPa);
    DHI_model = GHI_actual - cosd(90-SunEl).*DNI_model;
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
    temp=find(GHI_actual~=10000);
    ACPower(1:size(Time.hour,1),5) =10000;
    ACPower (temp,5)= pvl_snlinverter(Inverter, mSAPMResults.Vmp(temp)*Array.Ms, mSAPMResults.Pmp(temp)*Array.Ms*Array.Mp)*ninvert(k)/1000000;

    %ACPower 6th colummn has the power based on mean forecast GHI
    EdiffGround2 = pvl_grounddiffuse(Array.Tilt,GHI_forecast(:,1), Albedo);
    DNI_model2 = pvl_disc(GHI_forecast(:,1),90-SunEl, dayofyear,PresPa);
    DHI_model2 = GHI_forecast(:,1) - cosd(90-SunEl).*DNI_model2;
    Eb2 = 0*AOI; %Initiallize variable
    Eb2(AOI<90) = DNI_model2(AOI<90).*cosd(AOI(AOI<90)); %Only calculate when sun is in view of the plane of array
    EdiffSky2 = pvl_isotropicsky(Array.Tilt,DHI_model2);
    E2= Eb2 + EdiffSky2 + EdiffGround2; % Total incident irradiance (W/m^2)
    E02 = 1000; %Reference irradiance (1000 W/m^2)
    celltemp2 = pvl_sapmcelltemp(E2, E02, Array.a, Array.b,Wspd(:,1), Drybulb(:,1), ModuleParameters.delT);
    Ediff2 = EdiffSky2 + EdiffGround2; % Total diffuse incident irradiance (W/m^2)
    Ee2 = F1.*((Eb2.*F2+ModuleParameters.fd.*Ediff2)/E02)*SF; %Effective irradiance
    Ee2(isnan(Ee2))=0; % Set any NaNs to zero
    mSAPMResults2 = pvl_sapm(ModuleParameters, Ee2, celltemp2);
    aSAPMResults2.Vmp = Array.Ms  *mSAPMResults2.Vmp;
    aSAPMResults2.Imp = Array.Mp  *mSAPMResults2.Imp;
    aSAPMResults2.Pmp = aSAPMResults2.Vmp .* aSAPMResults2.Imp;
    clear temp
    temp=find(GHI_forecast~=10000);
    ACPower(1:size(Time.hour,1),6) =10000;
    ACPower(temp,6) = pvl_snlinverter(Inverter, mSAPMResults2.Vmp(temp)*Array.Ms, mSAPMResults2.Pmp(temp)*Array.Ms*Array.Mp)*ninvert(k)/1000000;
    
    
    %ACPower 7th column has the power based on upper forecast GHI
    EdiffGround3 = pvl_grounddiffuse(Array.Tilt,GHI_upperforecast(:,1), Albedo);
    DNI_model3 = pvl_disc(GHI_upperforecast(:,1),90-SunEl, dayofyear,PresPa);
    DHI_model3 = GHI_upperforecast(:,1) - cosd(90-SunEl).*DNI_model3;
    Eb3 = 0*AOI; %Initiallize variable
    Eb3(AOI<90) = DNI_model3(AOI<90).*cosd(AOI(AOI<90)); %Only calculate when sun is in view of the plane of array
    EdiffSky3 = pvl_isotropicsky(Array.Tilt,DHI_model3);
    E3= Eb3 + EdiffSky3 + EdiffGround3; % Total incident irradiance (W/m^2)
    E03 = 1000; %Reference irradiance (1000 W/m^2)
    celltemp3 = pvl_sapmcelltemp(E3, E03, Array.a, Array.b,Wspd(:,1), Drybulb(:,1), ModuleParameters.delT);
    Ediff3 = EdiffSky3 + EdiffGround3; % Total diffuse incident irradiance (W/m^2)
    Ee3 = F1.*((Eb3.*F2+ModuleParameters.fd.*Ediff3)/E03)*SF; %Effective irradiance
    Ee3(isnan(Ee3))=0; % Set any NaNs to zero
    mSAPMResults3 = pvl_sapm(ModuleParameters, Ee3, celltemp3);
    aSAPMResults3.Vmp = Array.Ms  *mSAPMResults3.Vmp;
    aSAPMResults3.Imp = Array.Mp  *mSAPMResults3.Imp;
    aSAPMResults3.Pmp = aSAPMResults3.Vmp .* aSAPMResults3.Imp;
    clear temp
    temp=find(GHI_upperforecast~=10000);
    ACPower(1:size(Time.hour,1),7) =10000;
    ACPower(temp,7) =pvl_snlinverter(Inverter, mSAPMResults3.Vmp(temp)*Array.Ms, mSAPMResults3.Pmp(temp)*Array.Ms*Array.Mp)*ninvert(k)/1000000;
    
    %ACPower 8th column has the power based on lower forecast GHI
    EdiffGround4 = pvl_grounddiffuse(Array.Tilt,GHI_lowerforecast(:,1), Albedo);
    DNI_model4 = pvl_disc(GHI_lowerforecast(:,1),90-SunEl, dayofyear,PresPa);
    DHI_model4 = GHI_lowerforecast(:,1) - cosd(90-SunEl).*DNI_model4;
    Eb4 = 0*AOI; %Initiallize variable
    Eb4(AOI<90) = DNI_model4(AOI<90).*cosd(AOI(AOI<90)); %Only calculate when sun is in view of the plane of array
    EdiffSky4 = pvl_isotropicsky(Array.Tilt,DHI_model4);
    E4= Eb4 + EdiffSky4 + EdiffGround4; % Total incident irradiance (W/m^2)
    E04 = 1000; %Reference irradiance (1000 W/m^2)
    celltemp4 = pvl_sapmcelltemp(E4, E04, Array.a, Array.b,Wspd(:,1), Drybulb(:,1), ModuleParameters.delT);
    Ediff4 = EdiffSky4 + EdiffGround4; % Total diffuse incident irradiance (W/m^2)
    Ee4 = F1.*((Eb4.*F2+ModuleParameters.fd.*Ediff4)/E04)*SF; %Effective irradiance
    Ee4(isnan(Ee4))=0; % Set any NaNs to zero
    mSAPMResults4 = pvl_sapm(ModuleParameters, Ee4, celltemp4);
    aSAPMResults4.Vmp = Array.Ms  *mSAPMResults4.Vmp;
    aSAPMResults4.Imp = Array.Mp  *mSAPMResults4.Imp;
    aSAPMResults4.Pmp = aSAPMResults4.Vmp .* aSAPMResults4.Imp;
    clear temp
    temp=find(GHI_lowerforecast~=10000);
    ACPower(1:size(Time.hour,1),8) =10000;
    ACPower(temp,8) = pvl_snlinverter(Inverter, mSAPMResults4.Vmp(temp)*Array.Ms, mSAPMResults4.Pmp(temp)*Array.Ms*Array.Mp)*ninvert(k)/1000000;
    
    %it consumes energy when no sun (for convenience i turn it to 0 for
    %now)
    ACPower(find(ACPower(:,6)<0),6)=0;
    ACPower(find(ACPower(:,5)<0),5)=0;
    ACPower(find(ACPower(:,7)<0),7)=0;
    ACPower(find(ACPower(:,8)<0),8)=0;
    xlswrite(['March2019\',char(sitenames(k)),'_Mar2019_power_10_percent.xlsx'], ACPower,'ACPower');
end
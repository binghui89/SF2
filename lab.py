# from datetime import date as Pdate, datetime as Ptime, timedelta
import sys, json, requests, os, urllib, StringIO, datetime, platform, time
import pandas as pd, numpy as np, matplotlib.lines as mlines, pvlib
from os.path import splitext
from zipfile import ZipFile, BadZipfile
from collections import defaultdict, OrderedDict
from matplotlib import gridspec, pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF
from IPython import embed as IP


def greatcircle_distance(lat1, lon1, lat2, lon2):
    # Return the distance between two point(lat1, lon1) and (lat2, lon2).
    # Input: in degree; Output: in km.
    # Reference: http://www.movable-type.co.uk/scripts/latlong.html
    # Reference: https://github.com/binghui89/GulfStream/blob/master/LC_pa_2step.m

    from math import sin, cos, sqrt, atan2, radians

    R = 6373.0

    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat2)
    lon2 = radians(lon2)

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    distance = R * c
    return distance


def dist_tx2kb_spdis():
    # Calculate distances between pairs of TX2kB solar PV site and Solar Power
    # Data for Integration Study site (NREL).

    def return_Actual(members, scale=None):
        # Remove day-ahead and hour-ahead forecasted data.
        valid_members = list()
        if scale:
            for m in members:
                fname, fextention = splitext(m)
                if fname.split('_')[0] == 'Actual' and scale in fname.split('_'):
                    valid_members.append(m)
        else:
            for m in members:
                fname, fextention = splitext(m)
                if fname.split('_')[0] == 'Actual':
                    valid_members.append(m)

        return valid_members

    solar_TX2kB = [
        [29.739, -94.992],
        [29.316, -100.38],
        [29.741, -94.992],
        [29.742, -94.992],
        [31.347, -94.74],
        [28.554, -96.821],
        [32.254, -101.419],
        [28.548, -96.821],
        [30.551, -97.689],
        [30.549, -97.689],
        [30.548, -97.689],
        [30.547, -97.689],
        [30.546, -97.689],
        [29.584, -104.296],
        [31.816, -94.888],
        [30.552, -97.689],
        [32.257, -101.419],
        [29.746, -94.992],
        [28.551, -96.821],
        [29.743, -94.992],
        [28.545, -96.821],
        [29.744, -94.992],
    ] # Solar PV generators from the TX2kB system

    fname = 'C:\\Users\\bxl180002\\Downloads\\SolarRamp\\tx-pv-2006.zip'

    z = ZipFile(fname, 'r')

    members = return_Actual( z.namelist())
    site = 0
    list_lat  = list()
    list_lon  = list()
    list_cap  = list()

    for m in members:
        site += 1
        _, lat, lon, year, _, capacity, _, _ = splitext(m)[0].split('_')
        capacity = int( capacity[0:-2] ) # 'MW' is the last two chars
        lat      = float(lat)
        lon      = float(lon)
        annual_output  = 0
        print 'site # {:03}, {} MW, {} MWh, CF = {}'.format(site, capacity, 
                                                annual_output, 
                                                annual_output/capacity/8760.0)
        list_lat.append(lat)
        list_lon.append(lon)
        list_cap.append(capacity)

    z.close()

    for i in range(0, len(solar_TX2kB)):
        d_ij = list()
        lati, loni = solar_TX2kB[i]
        for j in range(0, len(list_lat)):
            latj, lonj = list_lat[j], list_lon[j]
            d_ij.append( greatcircle_distance(lati, loni, latj, lonj) )
        print min(d_ij)
    IP()


def read_nsrdb():
    # lat, lon, year = 33.2164, -97.1292, 2010
    lat, lon, year = 32.89, -96.9, 2010
    api_key = 'KUVPDJWCgREC9NbiO6uOBDSHBAEjCYDeCRwa9L2S'
    # attributes = 'ghi,dhi,dni,wind_speed_10m_nwp,surface_air_temperature_nwp,solar_zenith_angle'
    attributes = 'ghi,dhi,dni'
    # Choose year of data
    year = '2010'
    # Set leap year to true or false. True will return leap day data if present, false will not.
    leap_year = 'false'
    # Set time interval in minutes, i.e., '30' is half hour intervals. Valid intervals are 30 & 60.
    interval = '30'
    # Specify Coordinated Universal Time (UTC), 'true' will use UTC, 'false' will use the local time zone of the data.
    # NOTE: In order to use the NSRDB data in SAM, you must specify UTC as 'false'. SAM requires the data to be in the
    # local time zone.
    utc = 'false'
    # Your full name, use '+' instead of spaces.
    your_name = 'Binghui+Li'
    # Your reason for using the NSRDB.
    reason_for_use = 'beta+testing'
    # Your affiliation
    your_affiliation = 'UT+Dallas'
    # Your email address
    your_email = 'bxl180002@utdallas.edu'
    # Please join our mailing list so we can keep you up-to-date on new developments.
    mailing_list = 'false'

    # Declare url string
    url = (
        'http://developer.nrel.gov/api/solar/nsrdb_psm3_download.csv?'
        'wkt=POINT({lon}%20{lat})'
        '&names={year}&'
        'leap_day={leap}'
        '&interval={interval}'
        '&utc={utc}'
        '&full_name={name}'
        '&email={email}'
        '&affiliation={affiliation}'
        '&mailing_list={mailing_list}'
        '&reason={reason}'
        '&api_key={api}'
        '&attributes={attr}'
    ).format(
        year=year, lat=lat, lon=lon, leap=leap_year, interval=interval, 
        utc=utc, name=your_name, email=your_email, mailing_list=mailing_list, 
        affiliation=your_affiliation, reason=reason_for_use, api=api_key, 
        attr=attributes
    )
    # Return just the first 2 lines to get metadata:
    info = pd.read_csv(url, nrows=1)
    # See metadata for specified properties, e.g., timezone and elevation
    timezone, elevation = info['Local Time Zone'], info['Elevation']
    
    df = pd.read_csv(url, skiprows=2)
    df = df.set_index(
        pd.date_range(
            '1/1/{yr}'.format(yr=year),
            freq=interval+'Min',
            periods=525600/int(interval),
        )
    )

    IP()

def process_nsrdb():
    # Prepare for cong
    from pytz import timezone

    ls_csv = [
        r'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\CA_Topaz\96633_35.37_-120.18_2018.csv',
        # r'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\COWC1\138221_39.13_-123.06_2018.csv',
        r'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\DEMC1\98317_35.53_-118.62_2018.csv',
        # r'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\KNNC1\157537_40.73_-123.94_2018.csv',
        r'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\MIAC1\118489_37.41_-119.74_2018.csv',
        r'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\MNCC1\85923_34.29_-117.5_2018.csv',
        # r'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\RLKC1\151675_40.25_-123.3_2018.csv',
        # r'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\RSAC1\130670_38.49_-122.7_2018.csv',
        # r'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\SBVC1\87449_34.45_-119.7_2018.csv',
        r'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\STFC1\84346_34.13_-117.94_2018.csv',
    ]
    dirsave = r'C:\Users\bxl180002\Downloads\RampSolar\NSRDB\forCong'

    for csvname in ls_csv:
        df = pd.read_csv(csvname, header=2)
        df['target'] = np.concatenate([df['GHI'].values[2:], [np.nan]*2])
        df['Cl_GHI'] = df['GHI']/df['Clearsky GHI']
        df['timestamp'] = pd.to_datetime(
            df['Year'].map('{:4g}'.format) 
            + 
            '-' 
            + 
            df['Month'].map('{:02g}'.format) 
            + 
            '-' 
            + 
            df['Day'].map('{:02g}'.format) 
            + 
            '-' 
            + 
            df['Hour'].map('{:02g}'.format) 
            + 
            '-' 
            + 
            df['Minute'].map('{:02g}'.format), 
            format='%Y-%m-%d-%H-%M'
        ).dt.tz_localize('UTC')
        df['timestamp_local'] = df['timestamp'].dt.tz_convert('US/Pacific')
        ts_sampl = datetime.datetime(2018, 6, 23, 0, 0, tzinfo=timezone('US/Pacific')) 
        ts_frcst = datetime.datetime(2018, 12, 23, 0, 0, tzinfo=timezone('US/Pacific')) 
        te_frcst = datetime.datetime(2018, 12, 30, 0, 0, tzinfo=timezone('US/Pacific')) 

        save_columns = [
            'Month',
            'Hour',
            'DHI',
            'DNI',
            'GHI',
            'Clearsky DHI',
            'Clearsky DNI',
            'Clearsky GHI',
            'Dew Point',
            'Solar Zenith Angle',
            'Wind Speed',
            'Relative Humidity',
            'Temperature',
            'Pressure',
            'Cl_GHI',
        ]

        csvdir = os.path.dirname(csvname)

        csv_sampl = os.path.sep.join(
            [
                dirsave,
                os.path.basename(csvdir)+'_train.csv',
            ]
        )
        df.loc[
            (df['timestamp_local']>=ts_sampl)&(df['timestamp_local']<ts_frcst),
            save_columns+['target'],
        ].to_csv(csv_sampl, index=False)

        csv_frcst = os.path.sep.join(
            [
                dirsave,
                os.path.basename(csvdir)+'_frcst.csv',
            ]
        )
        df.loc[
            (df['timestamp_local']>=ts_frcst)&(df['timestamp_local']<te_frcst),
            save_columns,
        ].to_csv(csv_frcst, index=False)

###########################################################
# IBM paris data starts here.

def site_name(lat, lon):
    coord_caiso = {
        "CA_Topaz": [35.38, -120.18],
        "RSAC1":    [38.47, -122.71],
        "RLKC1":    [40.25, -123.31],
        "SBVC1":    [34.45, -119.70],
        "KNNC1":    [40.71, -123.92],
        "MIAC1":    [37.41, -119.74],
        "MNCC1":    [34.31, -117.50],
        "STFC1":    [34.12, -117.94],
        "DEMC1":    [35.53, -118.63],
        "COWC1":    [39.12, -123.07],
    }

    df_caiso = pd.DataFrame(coord_caiso).T
    df_caiso.columns = ['lat', 'lon']

    site = df_caiso.loc[(df_caiso['lat']==lat) & (df_caiso['lon']==lon), :].index[0]
    return site

def process_paris(dir_work, write_flag=False):
    '''
    Process IBM's raw data, extract all layers, fill missing data with NA
    This is for the 15-min data
    '''

    # Sadly we cannot use the load_dir function, bummer.

    ls_csv, ls_df = load_dir(dir_work)
    os.chdir(dir_work)

    df_all = pd.concat(ls_df, ignore_index=True)
    df_all = df_all[['latitude', 'longitude', 'timestamp', 'layerName', 'value']].drop_duplicates().reset_index(drop=True)
    df_all.dropna(how='any', subset=['latitude', 'longitude'], inplace=True)
    # Stupid datatime only works for INDEX, and the 'dt' accessor must be used for Series
    df_all.loc[:, 'timestamp'] = pd.to_datetime(df_all['timestamp']).dt.tz_localize('UTC') #.dt.tz_convert(timezone)
    locations = df_all[['latitude', 'longitude']].drop_duplicates().values
    df_all.loc[:, 'TIME_STR'] = df_all['timestamp'].dt.strftime('%Y-%m-%d %H:%M')
    layers = df_all['layerName'].unique()
    layers_wanted = [
        'dswr_95_lower_ci',
        'dswr_mean',
        'dswr_95_upper_ci',
        # 'sp_mean',
        # 'sp_95_lower_ci',
        # 'sp_95_upper_ci', 
    ]
    ls_df_results = list()
    for lat, lon in locations:
        df_site = df_all.loc[
            (df_all['latitude'] == lat)
            &
            (df_all['longitude'] == lon)
            &
            (
                (df_all['layerName'] == 'dswr_mean')
                |
                (df_all['layerName'] == 'dswr_95_lower_ci')
                |
                (df_all['layerName'] == 'dswr_95_upper_ci')
            )
            ,
            :
        ] #.set_index('TIME_STR')
        index_datetime = pd.date_range(
            df_site['timestamp'].min(), 
            df_site['timestamp'].max(), 
            freq='15min',
        ) # .tz_localize('UTC').tz_convert(timezone)
        df_results = pd.DataFrame(
            index=index_datetime.strftime('%Y-%m-%d %H:%M')
        )
        df_results.loc[:, 'timestamp'] = index_datetime
        df_results.loc[:, 'Year']      = index_datetime.year
        df_results.loc[:, 'Month']     = index_datetime.month
        df_results.loc[:, 'Day']       = index_datetime.day
        df_results.loc[:, 'Hour']      = index_datetime.hour
        df_results.loc[:, 'Minute']    = index_datetime.minute
        # for l in layers:
        #     df_tmp = df_site.loc[df_site['layerName']==l, :].groupby(['TIME_STR']).mean()
        #     df_results.loc[:, l] = df_tmp['value']
        # IP()
        for l in layers_wanted:
            df_tmp = df_site.loc[df_site['layerName']==l, ['TIME_STR', 'value']].set_index('TIME_STR')
            df_results.loc[:, l] = df_tmp['value']
        ls_df_results.append(df_results)

        s = site_name(lat, lon)
        csvname = '_'.join(['IBM_processed', s]) + '.csv'
        col_to_write = [
            'Year', 
            'Month', 
            'Day', 
            'Hour', 
            'Minute',
            # 'dswr_mean', 
            # 'dswr_95_upper_ci', 
            # 'dswr_95_lower_ci', 
            # 'global_solar_irradiance'
        ] + layers_wanted
        if write_flag:
            df_results[col_to_write].to_csv(csvname, index=False)

    IP()

def process_paris_global_solar_irradiance(dir_work, write_flag=False):
    '''
    Only process actual GHI data, since it is hourly data.
    '''

    # Sadly we cannot use the load_dir function, bummer.
    # dir_work = r'C:\Users\bxl180002\Downloads\RampSolar\IBM_May'

    ls_csv, ls_df = load_dir(dir_work)
    os.chdir(dir_work)

    df_all = pd.concat(ls_df, ignore_index=True)
    df_all = df_all[['latitude', 'longitude', 'timestamp', 'layerName', 'value']].drop_duplicates().reset_index(drop=True)
    df_all.dropna(how='any', subset=['latitude', 'longitude'], inplace=True)
    # Stupid datatime only works for INDEX, and the 'dt' accessor must be used for Series
    df_all.loc[:, 'timestamp'] = pd.to_datetime(df_all['timestamp']).dt.tz_localize('UTC') #.dt.tz_convert(timezone)
    locations = df_all[['latitude', 'longitude']].drop_duplicates().values
    df_all.loc[:, 'TIME_STR'] = df_all['timestamp'].dt.strftime('%Y-%m-%d %H')
    layers = df_all['layerName'].unique()
    layers_wanted = ['global_solar_irradiance',]
    ls_df_results = list()
    for lat, lon in locations:
        df_site = df_all.loc[
            (df_all['latitude'] == lat) & (df_all['longitude'] == lon)
            &
            (df_all['layerName'] == 'global_solar_irradiance'),
            :
        ]
        index_datetime = pd.date_range(
            df_site['timestamp'].min(), 
            df_site['timestamp'].max(), 
            freq='1H',
        )
        df_results = pd.DataFrame(
            index=index_datetime.strftime('%Y-%m-%d %H')
        )
        df_results.loc[:, 'timestamp'] = index_datetime
        df_results.loc[:, 'Year']      = index_datetime.year
        df_results.loc[:, 'Month']     = index_datetime.month
        df_results.loc[:, 'Day']       = index_datetime.day
        df_results.loc[:, 'Hour']      = index_datetime.hour
        df_results.loc[:, 'Minute']    = index_datetime.minute
        for l in layers_wanted:
            df_tmp = df_site.loc[df_site['layerName']==l, ['TIME_STR', 'value']].set_index('TIME_STR')
            df_results.loc[:, l] = df_tmp['value']
        ls_df_results.append(df_results)

        s = site_name(lat, lon)
        csvname = '_'.join(['IBM_processed', s]) + '.hourly' + '.csv'
        col_to_write = [
            'Year', 
            'Month', 
            'Day', 
            'Hour', 
            'Minute', 
        ] + layers_wanted
        if write_flag:
            df_results[col_to_write].to_csv(csvname, index=False)
    IP()

def read_paris_April():
    # Somehow read_paris1 does not work... so I rewrite the whole thing.
    # This is the one that I used to collect IBM data in June 2019
    os.chdir('./IBM/April')
    df_results = pd.DataFrame()
    # From Rui's code
    server = 'https://pairs.res.ibm.com/v2/query'
    global_pairs_auth = ('bxl180002@utdallas.edu', '?HOoper40$')

    dict_layers = {
        "GHI mean":        {"type": "vector", "id": "P274C4230"}, # GHI mean 
        "GHI upper 95%":   {"type": "vector", "id": "P275C4234"}, # GHI upper 95% CI
        "GHI lower 95%":   {"type": "vector", "id": "P276C4238"}, # GHI lower 95% CI
        "DNI mean":        {"type": "vector", "id": "P280C4254"}, # DNI mean
        "DNI lower 95%":   {"type": "vector", "id": "P283C4266"}, # DNI lower 95% CI
        "DNI upper 95%":   {"type": "vector", "id": "P282C4262"}, # DNI upper 95% CI
        "DHI mean":        {"type": "vector", "id": "P284C4270"}, # DHI mean
        "DHI lower 95%":   {"type": "vector", "id": "P286C4278"}, # DHI lower 95% CI
        "DHI upper 95%":   {"type": "vector", "id": "P285C4274"}, # DHI upper 95% CI
        "Power mean":      {"type": "vector", "id": "P277C4242"}, # Power, mean
        "Power lower 95%": {"type": "vector", "id": "P279C4250"}, # Power, 95% low
        "Power upper 95%": {"type": "vector", "id": "P278C4246"}, # Power, 95% high
        "GHI actual":      {"type": "vector", "id": "P221C3743"}, # GHI, actual. Unfortunately they do not provide actual DNI and DHI.
    }

    coord_caiso = {
        "CA_Topaz": ["35.38", "-120.18"],
        "RSAC1":    ["38.47", "-122.71"],
        "RLKC1":    ["40.25", "-123.31"],
        "SBVC1":    ["34.45", "-119.70"],
        "KNNC1":    ["40.71", "-123.92"],
        "MIAC1":    ["37.41", "-119.74"],
        "MNCC1":    ["34.31", "-117.50"],
        "STFC1":    ["34.12", "-117.94"],
        "DEMC1":    ["35.53", "-118.63"],
        "COWC1":    ["39.12", "-123.07"],
    }

    coord_miso = {
        "AMOA4": ["33.58", "-91.80"],
        "FRMI4": ["40.64", "-91.72"],
        "BNRI2": ["37.24", "-89.37"],
        "SULI3": ["39.07", "-87.35"],
        "NATL1": ["31.49", "-93.19"],
        "BDLM4": ["42.62", "-85.65"],
        "CASM5": ["47.37", "-94.61"],
        "CKWM6": ["30.52", "-88.98"],
        "TS428": ["46.89", "-103.37"],
        "RHRS2": ["43.87", "-103.44"],
    }
    intervals = [
        {"start": "2019-04-01T00:00:00Z", "end": "2019-05-01T00:00:00Z"},
    ]

    query_json = {
        "layers": [
            {"type": "vector", "id": "P274C4230"}, # GHI mean 
            {"type": "vector", "id": "P275C4234"}, # GHI upper 95% CI
            {"type": "vector", "id": "P276C4238"}, # GHI lower 95% CI
            {"type": "vector", "id": "P280C4254"}, # DNI mean
            {"type": "vector", "id": "P283C4266"}, # DNI lower 95% CI
            {"type": "vector", "id": "P282C4262"}, # DNI upper 95% CI
            {"type": "vector", "id": "P284C4270"}, # DHI mean
            {"type": "vector", "id": "P286C4278"}, # DHI lower 95% CI
            {"type": "vector", "id": "P285C4274"}, # DHI upper 95% CI
            {"type": "vector", "id": "P277C4242"}, # Power, mean
            {"type": "vector", "id": "P279C4250"}, # Power, 95% low
            {"type": "vector", "id": "P278C4246"}, # Power, 95% high
            {"type": "vector", "id": "P221C3743"}, # GHI, actual. Unfortunately they do not provide actual DNI and DHI.
        ],
        "spatial": {
            "type": "point", 
            # "coordinates": ["35.38", "-120.18"], # CAISO Station 1
            "coordinates": ["38.47", "-122.71"], # CAISO Station 2
            # "coordinates": ["40.25", "-123.31"], # CAISO Station 3
            # "coordinates": ["34.45", "-119.70"], # CAISO Station 4
            # "coordinates": ["40.71", "-123.92"], # CAISO Station 5
            # "coordinates": ["37.41", "-119.74"], # CAISO Station 6
            # "coordinates": ["34.31", "-117.50"], # CAISO Station 7
            # "coordinates": ["34.12", "-117.94"], # CAISO Station 8
            # "coordinates": ["35.53", "-118.63"], # CAISO Station 9
            # "coordinates": ["39.12", "-123.07"], # CAISO Station 10
        },
        "temporal": {
            "intervals": [{"start": "2018-09-16T00:00:00Z", "end": "2018-09-30T00:00:00Z"}],
        }
    }

    for s in coord_caiso:
        c = coord_caiso[s]
        ls_df = list()
        query_json["spatial"]["coordinates"] = c
        for layer_name in dict_layers:
            query_json["temporal"]["intervals"] = [{"start": "2019-04-01T07:00:00Z", "end": "2019-05-01T07:00:00Z"}] # Starts and ends at 0 PDT
            query_json["layers"] = [dict_layers[layer_name]]
            response = requests.post(
                url=server,
                json=query_json, 
                auth=global_pairs_auth,
            )

            ntry = 1
            while (response.status_code != 200) & (ntry<=10):
                response = requests.post(
                    url=server,
                    json=query_json, 
                    auth=global_pairs_auth,
                )
                ntry = ntry+1
                time.sleep(5) # Let's wait for 5 seconds to give a second try

            if response.status_code == 200:
                ls_df.append(pd.DataFrame(response.json()['data']))
                print s, layer_name, 'Done!', 'ntry =', ntry
            else:
                print s, layer_name, 'no'

        data = pd.concat(ls_df, axis=0)
        data.loc[:, 'timestamp'] = pd.to_datetime(data['timestamp'], unit='ms')

        csvname = 'IBM_raw_' + s + '.csv'
        data.to_csv(csvname, index=False)

    IP()

def read_paris_May():
    # Somehow read_paris1 does not work... so I rewrite the whole thing.
    # This is the one that I used to collect IBM data in June 2019
    os.chdir('./IBM/May')
    df_results = pd.DataFrame()
    # From Rui's code
    server = 'https://pairs.res.ibm.com/v2/query'
    global_pairs_auth = ('bxl180002@utdallas.edu', '?HOoper40$')

    dict_layers = {
        "GHI mean":        {"type": "vector", "id": "P274C4230"}, # GHI mean 
        "GHI upper 95%":   {"type": "vector", "id": "P275C4234"}, # GHI upper 95% CI
        "GHI lower 95%":   {"type": "vector", "id": "P276C4238"}, # GHI lower 95% CI
        "DNI mean":        {"type": "vector", "id": "P280C4254"}, # DNI mean
        "DNI lower 95%":   {"type": "vector", "id": "P283C4266"}, # DNI lower 95% CI
        "DNI upper 95%":   {"type": "vector", "id": "P282C4262"}, # DNI upper 95% CI
        "DHI mean":        {"type": "vector", "id": "P284C4270"}, # DHI mean
        "DHI lower 95%":   {"type": "vector", "id": "P286C4278"}, # DHI lower 95% CI
        "DHI upper 95%":   {"type": "vector", "id": "P285C4274"}, # DHI upper 95% CI
        "Power mean":      {"type": "vector", "id": "P277C4242"}, # Power, mean
        "Power lower 95%": {"type": "vector", "id": "P279C4250"}, # Power, 95% low
        "Power upper 95%": {"type": "vector", "id": "P278C4246"}, # Power, 95% high
        "GHI actual":      {"type": "vector", "id": "P221C3743"}, # GHI, actual. Unfortunately they do not provide actual DNI and DHI.
    }

    coord_caiso = {
        "CA_Topaz": ["35.38", "-120.18"],
        "RSAC1":    ["38.47", "-122.71"],
        "RLKC1":    ["40.25", "-123.31"],
        "SBVC1":    ["34.45", "-119.70"],
        "KNNC1":    ["40.71", "-123.92"],
        "MIAC1":    ["37.41", "-119.74"],
        "MNCC1":    ["34.31", "-117.50"],
        "STFC1":    ["34.12", "-117.94"],
        "DEMC1":    ["35.53", "-118.63"],
        "COWC1":    ["39.12", "-123.07"],
    }

    coord_miso = {
        "AMOA4": ["33.58", "-91.80"],
        "FRMI4": ["40.64", "-91.72"],
        "BNRI2": ["37.24", "-89.37"],
        "SULI3": ["39.07", "-87.35"],
        "NATL1": ["31.49", "-93.19"],
        "BDLM4": ["42.62", "-85.65"],
        "CASM5": ["47.37", "-94.61"],
        "CKWM6": ["30.52", "-88.98"],
        "TS428": ["46.89", "-103.37"],
        "RHRS2": ["43.87", "-103.44"],
    }
    # intervals = [
    #     {"start": "2019-05-01T00:00:00Z", "end": "2019-06-01T00:00:00Z"},
    # ]

    query_json = {
        "layers": [
            {"type": "vector", "id": "P274C4230"}, # GHI mean 
            {"type": "vector", "id": "P275C4234"}, # GHI upper 95% CI
            {"type": "vector", "id": "P276C4238"}, # GHI lower 95% CI
            {"type": "vector", "id": "P280C4254"}, # DNI mean
            {"type": "vector", "id": "P283C4266"}, # DNI lower 95% CI
            {"type": "vector", "id": "P282C4262"}, # DNI upper 95% CI
            {"type": "vector", "id": "P284C4270"}, # DHI mean
            {"type": "vector", "id": "P286C4278"}, # DHI lower 95% CI
            {"type": "vector", "id": "P285C4274"}, # DHI upper 95% CI
            {"type": "vector", "id": "P277C4242"}, # Power, mean
            {"type": "vector", "id": "P279C4250"}, # Power, 95% low
            {"type": "vector", "id": "P278C4246"}, # Power, 95% high
            {"type": "vector", "id": "P221C3743"}, # GHI, actual. Unfortunately they do not provide actual DNI and DHI.
        ],
        "spatial": {
            "type": "point", 
            # "coordinates": ["35.38", "-120.18"], # CAISO Station 1
            "coordinates": ["38.47", "-122.71"], # CAISO Station 2
            # "coordinates": ["40.25", "-123.31"], # CAISO Station 3
            # "coordinates": ["34.45", "-119.70"], # CAISO Station 4
            # "coordinates": ["40.71", "-123.92"], # CAISO Station 5
            # "coordinates": ["37.41", "-119.74"], # CAISO Station 6
            # "coordinates": ["34.31", "-117.50"], # CAISO Station 7
            # "coordinates": ["34.12", "-117.94"], # CAISO Station 8
            # "coordinates": ["35.53", "-118.63"], # CAISO Station 9
            # "coordinates": ["39.12", "-123.07"], # CAISO Station 10
        },
        "temporal": {
            "intervals": [{"start": "2018-09-16T00:00:00Z", "end": "2018-09-30T00:00:00Z"}],
        }
    }

    for s in coord_caiso:
        c = coord_caiso[s]
        ls_df = list()
        query_json["spatial"]["coordinates"] = c
        for layer_name in dict_layers:
            query_json["temporal"]["intervals"] = [{"start": "2019-05-01T07:00:00Z", "end": "2019-06-01T07:00:00Z"}] # Starts and ends at 0 PDT
            query_json["layers"] = [dict_layers[layer_name]]
            response = requests.post(
                url=server,
                json=query_json, 
                auth=global_pairs_auth,
            )
            ntry = 1
            while (response.status_code != 200) & (ntry<=10):
                response = requests.post(
                    url=server,
                    json=query_json, 
                    auth=global_pairs_auth,
                )
                ntry = ntry+1
                time.sleep(5) # Let's wait for 5 seconds to give a second try

            if response.status_code == 200:
                ls_df.append(pd.DataFrame(response.json()['data']))
                print s, layer_name, 'Done!', 'ntry =', ntry
            else:
                print s, layer_name, 'no'
                # IP()

        data = pd.concat(ls_df, axis=0)
        data.loc[:, 'timestamp'] = pd.to_datetime(data['timestamp'], unit='ms')

        csvname = 'IBM_raw_' + s + '.csv'
        data.to_csv(csvname, index=False)

def query_paris(query_json):
    '''
    Query PAIRS based on query_json
    '''
    server = 'https://pairs.res.ibm.com/v2/query'
    # global_pairs_auth = ('bxl180002@utdallas.edu', '?HOoper40$') # Old from Nov. 2018
    # global_pairs_auth = {"name": "binghui", "email": "binghui.li@utdallas.edu", "password": "5,T9Wk*!", "groupId":239} # New from Aug. 2019
    global_pairs_auth = ("binghui.li@utdallas.edu", "5,T9Wk*!") # New from Aug. 2019
    response = requests.post(
        url=server,
        json=query_json, 
        auth=global_pairs_auth,
    )
    return response
    # if response.status_code == 200:
    #     ls_df.append(pd.DataFrame(response.json()['data']))
    #     print s, layer_name, 'Done!'
    # else:
    #     print s, layer_name, 'no'

def read_paris_5min():
    # Somehow read_paris1 does not work... so I rewrite the whole thing.
    os.chdir('./IBM/May.5min')
    df_results = pd.DataFrame()
    # From Rui's code
    server = 'https://pairs.res.ibm.com/v2/query'
    global_pairs_auth = ('bxl180002@utdallas.edu', '?HOoper40$')

    dict_layers = {
        "5 min forecast":  {
            "type": "vector",
            "id": "P512C5638",
             "dimensions": [{"name": "horizon", "value": 2}]
        },
    }

    coord_caiso = {
        "CA_Topaz": ["35.38", "-120.18"],
        "RSAC1":    ["38.47", "-122.71"],
        "RLKC1":    ["40.25", "-123.31"],
        "SBVC1":    ["34.45", "-119.70"],
        "KNNC1":    ["40.71", "-123.92"],
        "MIAC1":    ["37.41", "-119.74"],
        "MNCC1":    ["34.31", "-117.50"],
        "STFC1":    ["34.12", "-117.94"],
        "DEMC1":    ["35.53", "-118.63"],
        "COWC1":    ["39.12", "-123.07"],
    }

    coord_miso = {
        "AMOA4": ["33.58", "-91.80"],
        "FRMI4": ["40.64", "-91.72"],
        "BNRI2": ["37.24", "-89.37"],
        "SULI3": ["39.07", "-87.35"],
        "NATL1": ["31.49", "-93.19"],
        "BDLM4": ["42.62", "-85.65"],
        "CASM5": ["47.37", "-94.61"],
        "CKWM6": ["30.52", "-88.98"],
        "TS428": ["46.89", "-103.37"],
        "RHRS2": ["43.87", "-103.44"],
    }
    intervals = [
        {"start": "2019-04-01T00:00:00Z", "end": "2019-05-01T00:00:00Z"},
    ]

    query_json = {
        "layers": [
            {"type": "vector", "id": "P274C4230"}, # GHI mean 
            {"type": "vector", "id": "P275C4234"}, # GHI upper 95% CI
            {"type": "vector", "id": "P276C4238"}, # GHI lower 95% CI
            {"type": "vector", "id": "P280C4254"}, # DNI mean
            {"type": "vector", "id": "P283C4266"}, # DNI lower 95% CI
            {"type": "vector", "id": "P282C4262"}, # DNI upper 95% CI
            {"type": "vector", "id": "P284C4270"}, # DHI mean
            {"type": "vector", "id": "P286C4278"}, # DHI lower 95% CI
            {"type": "vector", "id": "P285C4274"}, # DHI upper 95% CI
            {"type": "vector", "id": "P277C4242"}, # Power, mean
            {"type": "vector", "id": "P279C4250"}, # Power, 95% low
            {"type": "vector", "id": "P278C4246"}, # Power, 95% high
            {"type": "vector", "id": "P221C3743"}, # GHI, actual. Unfortunately they do not provide actual DNI and DHI.
        ],
        "spatial": {
            "type": "point", 
            # "coordinates": ["35.38", "-120.18"], # CAISO Station 1
            "coordinates": ["38.47", "-122.71"], # CAISO Station 2
            # "coordinates": ["40.25", "-123.31"], # CAISO Station 3
            # "coordinates": ["34.45", "-119.70"], # CAISO Station 4
            # "coordinates": ["40.71", "-123.92"], # CAISO Station 5
            # "coordinates": ["37.41", "-119.74"], # CAISO Station 6
            # "coordinates": ["34.31", "-117.50"], # CAISO Station 7
            # "coordinates": ["34.12", "-117.94"], # CAISO Station 8
            # "coordinates": ["35.53", "-118.63"], # CAISO Station 9
            # "coordinates": ["39.12", "-123.07"], # CAISO Station 10
        },
        "temporal": {
            "intervals": [{"start": "2018-09-16T00:00:00Z", "end": "2018-09-30T00:00:00Z"}],
        }
    }

    for s in coord_caiso:
        c = coord_caiso[s]
        ls_df = list()
        query_json["spatial"]["coordinates"] = c
        for layer_name in dict_layers:
            query_json["temporal"]["intervals"] = [{"start": "2019-05-01T07:00:00Z", "end": "2019-06-01T07:00:00Z"}] # Start and end at 0 PDT this time
            query_json["layers"] = [dict_layers[layer_name]]
            response = requests.post(
                url=server,
                json=query_json, 
                auth=global_pairs_auth,
            )
            if response.status_code == 200:
                ls_df.append(pd.DataFrame(response.json()['data']))
                print s, layer_name, 'Done!'
            else:
                print s, layer_name, 'no'

        data = pd.concat(ls_df, axis=0)
        data.loc[:, 'timestamp'] = pd.to_datetime(data['timestamp'], unit='ms')

        csvname = 'IBM_raw_' + s + '.csv'
        data.to_csv(csvname, index=False)

def process_paris_more_quantiles(dir_work, write_flag=False):
    '''
    Process IBM's raw data, extract all layers, fill missing data with NA
    This is for the data with more quantiles
    '''

    # Sadly we cannot use the load_dir function, bummer.

    list_colname =  ['p5',   'p25',  'p50', 'p75',  'p95',  'Mean']
    list_pstr =     ['0.05', '0.25', '0.5', '0.75', '0.95', '-1']

    ls_csv, ls_df = load_dir(dir_work)
    os.chdir(dir_work)

    df_all = pd.concat(ls_df, ignore_index=True)
    df_all = df_all[['latitude', 'longitude', 'timestamp', 'layerName', 'value', 'property']].drop_duplicates().reset_index(drop=True)
    df_all.dropna(how='any', subset=['latitude', 'longitude'], inplace=True)
    # Stupid datatime only works for INDEX, and the 'dt' accessor must be used for Series
    df_all.loc[:, 'timestamp'] = pd.to_datetime(df_all['timestamp']).dt.tz_localize('UTC') #.dt.tz_convert(timezone)
    locations = df_all[['latitude', 'longitude']].drop_duplicates().values
    df_all.loc[:, 'TIME_STR'] = df_all['timestamp'].dt.strftime('%Y-%m-%d %H:%M')
    ls_df_results = list()
    for lat, lon in locations:
        df_site = df_all.loc[
            (df_all['latitude'] == lat)
            &
            (df_all['longitude'] == lon)
            ,
            :
        ] #.set_index('TIME_STR')

        # Find the property name of each percentile in this site, it changes with site... This should be a one-to-one mapping.
        property_all = df_site['property'].unique().tolist()
        map_property2pstr = dict()
        for pstr in list_pstr:
            for p in property_all:
                if pstr in p:
                    map_property2pstr[p] = pstr

        index_datetime = pd.date_range(
            df_site['timestamp'].min(), 
            df_site['timestamp'].max(), 
            freq='5min',
        ) # .tz_localize('UTC').tz_convert(timezone)
        df_results = pd.DataFrame(
            index=index_datetime.strftime('%Y-%m-%d %H:%M')
        )
        df_results.loc[:, 'timestamp'] = index_datetime
        df_results.loc[:, 'Year']      = index_datetime.year
        df_results.loc[:, 'Month']     = index_datetime.month
        df_results.loc[:, 'Day']       = index_datetime.day
        df_results.loc[:, 'Hour']      = index_datetime.hour
        df_results.loc[:, 'Minute']    = index_datetime.minute
        # for l in layers:
        #     df_tmp = df_site.loc[df_site['layerName']==l, :].groupby(['TIME_STR']).mean()
        #     df_results.loc[:, l] = df_tmp['value']
        # IP()
        for l in property_all:
            pstr = map_property2pstr[l]
            df_tmp = df_site.loc[df_site['property']==l, ['TIME_STR', 'value']].set_index('TIME_STR')
            colname = list_colname[list_pstr.index(pstr)]
            df_results.loc[:, colname] = df_tmp['value']
        ls_df_results.append(df_results)

        s = site_name(lat, lon)
        csvname = '_'.join(['IBM_processed', s]) + '.csv'
        col_to_write = [
            'Year', 
            'Month', 
            'Day', 
            'Hour', 
            'Minute',
            # 'dswr_mean', 
            # 'dswr_95_upper_ci', 
            # 'dswr_95_lower_ci', 
            # 'global_solar_irradiance'
        ] + list_colname
        if write_flag:
            df_results[col_to_write].to_csv(csvname, index=False)

        ls_vio = [0]*4
        ls_vio[0] = df_results.loc[(df_results['p5']  >= df_results['p25'])&(df_results['p25']>0), :].shape[0]
        ls_vio[1] = df_results.loc[(df_results['p25'] >= df_results['p50'])&(df_results['p50']>0), :].shape[0]
        ls_vio[2] = df_results.loc[(df_results['p50'] >= df_results['p75'])&(df_results['p75']>0), :].shape[0]
        ls_vio[3] = df_results.loc[(df_results['p75'] >= df_results['p95'])&(df_results['p95']>0), :].shape[0]
        print "Violations:", s, ls_vio

###########################################################
# This is a bit weird, but this function is to generate a time series for the 
# AGC data in the wind ramp project by using the ACE data from CAISO

def generate_agc_for_118():
    from scipy.stats import norm
    # Load ACE data
    df_ace = pd.read_csv('for_reg_ACE.csv')
    df_ace['T'] = pd.to_datetime(df_ace['T'])

    # Load total load
    df_load = pd.read_csv('for_flexiramp_summary.csv')
    
    # Select load in 2018 and calculate ACE in percentage
    i_selected = (pd.to_datetime(df_load['OPR_DT']).dt.year==2018)
    load = df_load.loc[i_selected, 'LOAD_B_RTD'].to_numpy().repeat(5) # Each interval has 5 minutes
    ace  = df_ace['ACE'].to_numpy()
    ace_percentage = ace/load

    # Calculate mean and standard deviation
    mu = ace_percentage[~np.isnan(ace_percentage)].mean()
    std = ace_percentage[~np.isnan(ace_percentage)].std()
    print 'MEAN and STD: ', mu, std

    # Plot it out
    plt.hist(ace_percentage[~np.isnan(ace_percentage)], bins=50, density=True, alpha=0.6, color='g')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)

    # Now, load the 118 system data
    t_agc = 6 # 6 seconds
    n_agc = 300/t_agc # Number of agc intervals per ed interval
    df_load_ed = pd.read_csv(r'C:\Users\bxl180002\git\FlexibleRampSCUC\118bus\ed_loads.csv')
    df_gen_ed  = pd.read_csv(r'C:\Users\bxl180002\git\FlexibleRampSCUC\118bus\ed_generator.csv')
    
    # Produce load AGC data
    dict_load_agc = dict()
    ls_bus = [i for i in df_load_ed.columns if i.startswith('Bus')]
    for b in ls_bus:
        load_ed = df_load_ed[b].to_numpy()
        len_load = load_ed.shape[0]
        load_agc = np.interp(np.linspace(0, 288, 288*n_agc), range(0, 289), load_ed.tolist()+ [load_ed[-1]])
        rand_scaler = np.ones(288*n_agc) + np.random.randn(288*n_agc)*std
        dict_load_agc[b] = load_agc*rand_scaler
    df_load_agc = pd.DataFrame(dict_load_agc)
    df_load_agc['Slot'] = df_load_agc.index+1
    plt.figure()
    plt.plot(
        range(0, 288), 
        df_load_ed[ls_bus].sum(axis=1), 
        label='ED',
    )
    plt.plot(
        np.linspace(0, 288, 288*n_agc), 
        df_load_agc[ls_bus].sum(axis=1),
        label='AGC',
    )
    plt.legend()
    plt.title('LOAD')

    # Produce renewable AGC data, note: only for non-thermal gens
    dict_gen_agc = dict()
    ls_gen = [
        i for i in df_gen_ed.columns 
        if i.startswith('Wind') 
        or i.startswith('Solar')
        or i.startswith('Hydro')
    ]
    for g in ls_gen:
        power_ed = df_gen_ed[g].to_numpy()
        len_power = power_ed.shape[0]
        power_agc = np.interp(np.linspace(0, 288, 288*n_agc), range(0, 289), power_ed.tolist()+ [power_ed[-1]])
        rand_scaler = np.ones(288*n_agc) + np.random.randn(288*n_agc)*std
        dict_gen_agc[g] = power_agc*rand_scaler
    df_gen_agc = pd.DataFrame(dict_gen_agc)
    df_gen_agc['Slot'] = df_gen_agc.index+1
    plt.figure()
    plt.plot(
        range(0, 288), 
        df_gen_ed[ls_gen].sum(axis=1),
        label='ED',
    )
    plt.plot(
        np.linspace(0, 288, 288*n_agc), 
        df_gen_agc[ls_gen].sum(axis=1),
        label='AGC',
    )
    plt.legend()
    plt.title('Renewables total')
    plt.show()

    IP()


###########################################################
# CAISO OASIS data collection starts here
def collect_caiso_api(dir_work, d_begin, d_end, urlparam, days_interval=29):
    dir_home = os.getcwd()
    os.chdir(dir_work)
    url = 'http://oasis.caiso.com/oasisapi/SingleZip'

    n_interval = (d_end - d_begin).days/days_interval
    interval = datetime.timedelta(days=days_interval)
    ls_df = list()
    query_err = False
    for i in range(0, n_interval+1):
        csvname = str(i) + '.csv'
        if os.path.isfile(csvname):
            continue
        if query_err:
            pass
        else:
            d1 = d_begin + interval*i
            d2 = d_begin + interval*(i+1)
            if d2 > d_end:
                d2 = d_end
            if d1 > d_end:
                break
        print "Time: {:>s}, {:>s} to {:>s}: Reading...".format(
            datetime.datetime.today().strftime('%H:%M'),
            d1.strftime('%Y%m%d %H:%M'),
            d2.strftime('%Y%m%d %H:%M'),
        )
        urlparam['startdatetime'] = d1.strftime('%Y%m%dT%H:00-0000')  # UTC time
        urlparam['enddatetime']   = d2.strftime('%Y%m%dT%H:00-0000')  # UTC time
        fullurl = url + '?' + urllib.urlencode(urlparam)
        response = urllib.urlopen(fullurl)
        s = response.read()
        str_file = StringIO.StringIO(s)

        # try:
        #     z = ZipFile(str_file)
        # except BadZipfile:
        #     time.sleep(5)
        #     query_err = True
        #     continue
        z = ZipFile(str_file)

        m = z.namelist()[0]  # Single zip, so only one member
        with z.open(m, 'r') as f:
            df_csv = pd.read_csv(f)
            ls_df.append(df_csv)

        # Sadly CAISO only allows for one query every 5 seconds...
        time.sleep(5)

        # # Write csv files
        df_csv.to_csv(csvname, index=False)

    os.chdir(dir_home)
    return ls_df

def collect_AS_REQ_or_RESULTS(whichone='req'):
    # Collect as requirements (whichone='req') and results (whichone='res').
    # Note that PST is UTC-8, PDT is UTC-7
    if platform.system() == 'Windows':
        dir_work = (
            'C:\\Users\\bxl180002\\Downloads\\RampSolar\\CAISO\\AS_REQ' 
            if whichone == 'req' 
            else 'C:\\Users\\bxl180002\\Downloads\\RampSolar\\CAISO\\AS_RESULTS'
        )
    elif platform.system() == 'Linux':
        pass  # To be added if running on the hpc
    else:
        pass

    queryname = 'AS_REQ' if whichone == 'req' else 'AS_RESULTS'
    urlparam = {
        'queryname':     queryname,  # Mandatory for singlezip
        'startdatetime': None,       # Mandatory for singlezip
        'enddatetime':   None,       # Mandatory for singlezip
        'market_run_id': 'ALL',      # Mandatory for singlezip
        'version':       1,          # Mandatory for singlezip
        'as_type':       'ALL',
        'as_region':     'ALL',
        'resultformat':  6,  # Let's just use csv format, 6 for csv, 1 for xml
    }

    d_begin = datetime.datetime(2017, 1, 1, 8, 0)  # This is PST 2017-01-01 00:00
    d_end   = datetime.datetime(2019, 6, 1, 7, 0)  # This is PDT 2019-05-31 24:00
    if not os.path.isdir(dir_work):
        os.mkdir(dir_work)
    ls_df = collect_caiso_api(dir_work, d_begin, d_end, urlparam)

    IP()

def collect_SLD_FCST_rtd():
    # Collect demand forecast for RTD market
    # Note that PST is UTC-8, PDT is UTC-7

    if platform.system() == 'Windows':
        dir_work = 'C:\\Users\\bxl180002\\Downloads\\RampSolar\\CAISO\\SLD_FCST_rtd'
    elif platform.system() == 'Linux':
        pass  # To be added if running on the hpc
    else:
        dir_work = '/Users/bli/Downloads/caiso'

    urlparam = {
        'queryname':      'SLD_FCST',  # Mandatory for singlezip
        'startdatetime':   None,  # Mandatory for singlezip
        'enddatetime':     None,  # Mandatory for singlezip
        'market_run_id':  'ALL',  # Mandatory for singlezip, all market types
        'execution_type': 'RTD',
        'version':        1,  # Mandatory for singlezip
        'resultformat':   6,  # Let's just use csv format, 6 for csv, 1 for xml
    }
    d_begin = datetime.datetime(2014, 1, 1, 8, 0)  # This is PST 2014-01-01 00:00
    d_end   = datetime.datetime(2019, 6, 1, 7, 0)  # This is PDT 2019-05-31 24:00
    if not os.path.isdir(dir_work):
        os.mkdir(dir_work)
    ls_df = collect_caiso_api(dir_work, d_begin, d_end, urlparam)

def collect_SLD_FCST_rtpd():
    # Collect demand forecast for RTPD
    # Note that PST is UTC-8, PDT is UTC-7
    # Some how 'ALL' does not includes RTPD forecast in collect_SLD_FCST, hence 
    # this function...

    if platform.system() == 'Windows':
        dir_work = 'C:\\Users\\bxl180002\\Downloads\\RampSolar\\CAISO\\SLD_FCST_rtpd'
    elif platform.system() == 'Linux':
        pass  # To be added if running on the hpc
    else:
        dir_work = '/Users/bli/Downloads/caiso'

    urlparam = {
        'queryname':      'SLD_FCST',  # Mandatory for singlezip
        'startdatetime':   None,  # Mandatory for singlezip
        'enddatetime':     None,  # Mandatory for singlezip
        'market_run_id':  'RTM',  # Mandatory for singlezip, all market types
        'execution_type': 'RTPD',
        'version':        1,  # Mandatory for singlezip
        'resultformat':   6,  # Let's just use csv format, 6 for csv, 1 for xml
    }
    d_begin = datetime.datetime(2015, 1, 1, 8, 0)  # This is PST 2015-01-01 00:00
    d_end   = datetime.datetime(2019, 6, 1, 7, 0)  # This is PDT 2019-05-31 24:00
    if not os.path.isdir(dir_work):
        os.mkdir(dir_work)
    ls_df = collect_caiso_api(dir_work, d_begin, d_end, urlparam, days_interval=29)

def collect_SLD_REN_FCST_rtd():
    # Collect solar and wind forecast
    # Note that PST is UTC-8, PDT is UTC-7

    if platform.system() == 'Windows':
        dir_work = 'C:\\Users\\bxl180002\\Downloads\\RampSolar\\CAISO\\SLD_REN_FCST_rtd'
    elif platform.system() == 'Linux':
        pass  # To be added if running on the hpc
    else:
        dir_work = '/Users/bli/Downloads/caiso'

    urlparam = {
        'queryname':      'SLD_REN_FCST',  # Mandatory for singlezip
        'startdatetime':   None,  # Mandatory for singlezip
        'enddatetime':     None,  # Mandatory for singlezip
        # Mandatory for singlezip, out of some unknown reason, if ALL market 
        # runs is selected, the response url automatically change to groupzip, 
        # then only the first day's data is downloaded, and we'll get an error 
        # after a random number of iterations. That's why only the RTD data is 
        # downloaded 
        'market_run_id':  'RTD',
        # 'execution_type': 'ALL',
        'version':        1,  # Mandatory for singlezip
        'resultformat':   6,  # Let's just use csv format, 6 for csv, 1 for xml
    }
    # When I start download this data, somehow 2014 does not work, so I have to use 2015...
    d_begin = datetime.datetime(2015, 1, 1, 8, 0)  # This is PST 2015-01-01 00:00
    d_end   = datetime.datetime(2019, 6, 1, 7, 0)  # This is PDT 2019-05-31 24:00
    if not os.path.isdir(dir_work):
        os.mkdir(dir_work)
    ls_df = collect_caiso_api(dir_work, d_begin, d_end, urlparam)

def collect_SLD_REN_FCST_rtpd():
    # Collect solar and wind forecast
    # Note that PST is UTC-8, PDT is UTC-7

    if platform.system() == 'Windows':
        dir_work = 'C:\\Users\\bxl180002\\Downloads\\RampSolar\\CAISO\\SLD_REN_FCST_rtpd'
    elif platform.system() == 'Linux':
        pass  # To be added if running on the hpc
    else:
        dir_work = '/Users/bli/Downloads/caiso'

    urlparam = {
        'queryname':      'SLD_REN_FCST',  # Mandatory for singlezip
        'startdatetime':   None,  # Mandatory for singlezip
        'enddatetime':     None,  # Mandatory for singlezip
        # Mandatory for singlezip, out of some unknown reason, if ALL market 
        # runs is selected, the response url automatically change to groupzip, 
        # then only the first day's data is downloaded, and we'll get an error 
        # after a random number of iterations. That's why only the RTD data is 
        # downloaded 
        'market_run_id':  'RTPD',
        # 'execution_type': 'ALL',
        'version':        1,  # Mandatory for singlezip
        'resultformat':   6,  # Let's just use csv format, 6 for csv, 1 for xml
    }
    # When I start download this data, somehow 2014 does not work, so I have to use 2015...
    d_begin = datetime.datetime(2015, 1, 1, 8, 0)  # This is PST 2017-01-01 00:00
    d_end   = datetime.datetime(2019, 6, 1, 7, 0)  # This is PDT 2019-05-31 24:00
    if not os.path.isdir(dir_work):
        os.mkdir(dir_work)
    ls_df = collect_caiso_api(dir_work, d_begin, d_end, urlparam)

def collect_SLD_ADV_FCST():
    # Collect advisory demand forecast for RTD and RTPD
    # Note that PST is UTC-8, PDT is UTC-7
    if platform.system() == 'Windows':
        dir_work = 'C:\\Users\\bxl180002\\Downloads\\RampSolar\\CAISO\\SLD_ADV_FCST'
    elif platform.system() == 'Linux':
        pass  # To be added if running on the hpc
    else:
        dir_work = '/Users/bli/Downloads/caiso'

    urlparam = {
        'queryname':      'SLD_ADV_FCST',  # Mandatory for singlezip
        'startdatetime':   None,  # Mandatory for singlezip
        'enddatetime':     None,  # Mandatory for singlezip
        'version':        4,  # Mandatory for singlezip
        'resultformat':   6,  # Let's just use csv format, 6 for csv, 1 for xml
    }
    # Sadly CAISO only provides this data up to 2017...
    d_begin = datetime.datetime(2017, 1, 1, 8, 0)  # This is PST 2017-01-01 00:00
    d_end   = datetime.datetime(2019, 6, 1, 7, 0)  # This is PDT 2019-05-31 24:00
    if not os.path.isdir(dir_work):
        os.mkdir(dir_work)
    ls_df = collect_caiso_api(dir_work, d_begin, d_end, urlparam)

def collect_ENE_FLEX_RAMP_REQT():
    # Collect flexible ramp reserve requirements for RTD and RTPD
    # Note that PST is UTC-8, PDT is UTC-7
    # Flexible ramp products example URL
    url_example = (
        'http://oasis.caiso.com/oasisapi/SingleZip?queryname=ENE_FLEX_RAMP_REQT&'
        'market_run_id=ALL&baa_grp_id=ALL&startdatetime=20180301T07:00-0000&'
        'enddatetime=20180302T07:00-0000&version=4&resultformat=6'
    ) # This is from CAISO's OASIS manual

    if platform.system() == 'Windows':
        dir_work = 'C:\\Users\\bxl180002\\Downloads\\RampSolar\\CAISO\\ENE_FLEX_RAMP_REQT'
    elif platform.system() == 'Linux':
        pass  # To be added if running on the hpc
    else:
        dir_work = '/Users/bli/Downloads/caiso'

    urlparam = {
        'queryname':      'ENE_FLEX_RAMP_REQT',  # Mandatory for singlezip
        'startdatetime':   None,  # Mandatory for singlezip
        'enddatetime':     None,  # Mandatory for singlezip
        'version':        4,  # Mandatory for singlezip
        'resultformat':   6,  # Let's just use csv format, 6 for csv, 1 for xml
        'baa_grp_id':     'ALL',
        'market_run_id':  'ALL',
    }
    # Sadly CAISO only provides this data up to 2017...
    d_begin = datetime.datetime(2017, 1, 1, 8, 0)  # This is PST 2017-01-01 00:00
    d_end   = datetime.datetime(2019, 6, 1, 7, 0)  # This is PDT 2019-05-31 24:00
    if not os.path.isdir(dir_work):
        os.mkdir(dir_work)
    ls_df = collect_caiso_api(dir_work, d_begin, d_end, urlparam)

###########################################################
# CAISO OASIS and ACE raw data process starts here
def load_dir(dir_work):
    # Read all csv files in dir_work and return a list of dataframes
    dir_work = os.path.abspath(dir_work)
    ls_csv  = list()
    ls_df   = list()
    if os.path.isdir(dir_work): # Is a directory with many csv files
        dir_home = os.getcwd()
        os.chdir(dir_work)
        allfiles = os.listdir(dir_work)
        for fname in allfiles:
            if fname.endswith('csv'):
                df_csv = pd.read_csv(fname)
                ls_csv.append(fname)
                ls_df.append(df_csv)
        os.chdir(dir_home)
        return ls_csv, ls_df
    else: # Is a zipfile
        z = ZipFile(dir_work)
        for m in z.namelist():
            if m.endswith('csv'):
                with z.open(m, 'r') as f:
                    df_csv = pd.read_csv(f)
                    ls_csv.append(m)
                    ls_df.append(df_csv)
        return ls_csv, ls_df

def process_raw_for_flexiramp():
    '''
    This function reads the raw data from CAISO's OASIS, extracts load 
    forecasts (binding and advisory), renewalbe forecasts (B), flexible ramp 
    reserve requirements and return two pandas dataframes: RTD and RTPD.
    '''
    # Load OASIS data, renewable forecast
    _, ls_df = load_dir('./CAISO_OASIS/SLD_REN_FCST_rtd.zip')
    df_ren_f = pd.concat(ls_df, ignore_index=True)
    _, ls_df = load_dir('./CAISO_OASIS/SLD_REN_FCST_rtpd.zip')
    ls_df.append(df_ren_f)
    df_ren_f = pd.concat(ls_df, ignore_index=True)

    # Load OASIS data, binding load forecast
    _, ls_df = load_dir('./CAISO_OASIS/SLD_FCST_rtd.zip')
    df_load_bf = pd.concat(ls_df, ignore_index=True)
    _, ls_df = load_dir('./CAISO_OASIS/SLD_FCST_rtpd.zip')
    ls_df.append(df_load_bf)
    df_load_bf = pd.concat(ls_df, ignore_index=True)

    # Load OASIS data, advisory load forecast
    _, ls_df = load_dir('./CAISO_OASIS/SLD_ADV_FCST.zip')
    df_load_af = pd.concat(ls_df, ignore_index=True)
    df_load_af = df_load_af.rename(columns={'VALUE': 'MW'})

    # Load OASIS data, flexible ramp reserve requirements
    _, ls_df = load_dir('./CAISO_OASIS/ENE_FLEX_RAMP_REQT.zip')
    df_frrq = pd.concat(ls_df, ignore_index=True)

    print('Data loaded!')

    # Infer timestamp for all data, note CAISO's original data is in UTC time
    df_ren_f.loc[:, 'TIMESTAMP'] = pd.to_datetime(df_ren_f['INTERVALENDTIME_GMT'], infer_datetime_format=True)
    df_load_bf.loc[:, 'TIMESTAMP'] = pd.to_datetime(df_load_bf['INTERVALENDTIME_GMT'], infer_datetime_format=True)
    df_load_af.loc[:, 'TIMESTAMP'] = pd.to_datetime(df_load_af['INTERVAL_END_GMT'], infer_datetime_format=True)
    df_frrq.loc[:, 'TIMESTAMP'] = pd.to_datetime(df_frrq['INTERVALENDTIME_GMT'], infer_datetime_format=True)

    # Now, processing RTPD data
    ############################################################################
    t_start_rtpd = min([
        df_frrq.loc[df_frrq['MARKET_RUN_ID']=='RTPD', 'TIMESTAMP'].min(),
        df_load_af.loc[df_load_af['MKT_TYPE']=='RTPD', 'TIMESTAMP'].min(),
        df_load_bf.loc[df_load_bf['EXECUTION_TYPE']=='RTPD', 'TIMESTAMP'].min(),
        df_ren_f.loc[df_ren_f['MARKET_RUN_ID']=='RTPD', 'TIMESTAMP'].min(),
    ])
    t_end_rtpd = max([
        df_frrq.loc[df_frrq['MARKET_RUN_ID']=='RTPD', 'TIMESTAMP'].max(),
        df_load_af.loc[df_load_af['MKT_TYPE']=='RTPD', 'TIMESTAMP'].max(),
        df_load_bf.loc[df_load_bf['EXECUTION_TYPE']=='RTPD', 'TIMESTAMP'].max(),
        df_ren_f.loc[df_ren_f['MARKET_RUN_ID']=='RTPD', 'TIMESTAMP'].max(),
    ])
    df_results_rtpd = pd.DataFrame(
        index=pd.date_range(start=t_start_rtpd, end=t_end_rtpd, freq='15min')
    )

    df_tmp = df_load_bf.loc[
        (df_load_bf['TAC_AREA_NAME']=='CA ISO-TAC')&(df_load_bf['EXECUTION_TYPE']=='RTPD'), 
        ['MW', 'TIMESTAMP']
    ].set_index('TIMESTAMP')
    df_results_rtpd.loc[:, 'LOAD_B_RTPD'] = df_tmp['MW']

    df_tmp = df_load_af.loc[
        (df_load_af['BAA_GRP_ID']=='CA ISO-TAC')&(df_load_af['MKT_TYPE']=='RTPD'), 
        ['MW', 'TIMESTAMP']
    ].set_index('TIMESTAMP')
    df_results_rtpd.loc[:, 'LOAD_A_RTPD'] = df_tmp['MW']

    df_ren_f.loc[:, 'R_H'] = df_ren_f['RENEWABLE_TYPE'].astype('str') + '_' + df_ren_f['TRADING_HUB'].astype('str')
    for rh in df_ren_f['R_H'].unique():
        c = '_'.join([rh, 'RTPD'])
        df_tmp = df_ren_f.loc[
            (df_ren_f['R_H']==rh)&(df_ren_f['MARKET_RUN_ID']=='RTPD'),
            ['TIMESTAMP', 'MW']
        ].set_index('TIMESTAMP')
        df_results_rtpd.loc[:, c] = df_tmp['MW']

    for ud in df_frrq['RAMP_TYPE'].unique():
        c = '_'.join([ud, 'RTPD'])
        df_tmp = df_frrq.loc[
            (df_frrq['RAMP_TYPE']==ud)&(df_frrq['MARKET_RUN_ID']=='RTPD')&(df_frrq['BAA_GRP_ID']=='CISO'),
            ['TIMESTAMP', 'UNCERTAINITY']
        ].set_index('TIMESTAMP')
        df_results_rtpd.loc[:, c] = df_tmp['UNCERTAINITY']

    # Now, processing RTD data
    ############################################################################
    t_start_rtd = min([
        df_frrq.loc[df_frrq['MARKET_RUN_ID']=='RTD', 'TIMESTAMP'].min(),
        df_load_af.loc[df_load_af['MKT_TYPE']=='RTD', 'TIMESTAMP'].min(),
        df_load_bf.loc[df_load_bf['EXECUTION_TYPE']=='RTD', 'TIMESTAMP'].min(),
        df_ren_f.loc[df_ren_f['MARKET_RUN_ID']=='RTD', 'TIMESTAMP'].min(),
    ])
    t_end_rtd = max([
        df_frrq.loc[df_frrq['MARKET_RUN_ID']=='RTD', 'TIMESTAMP'].max(),
        df_load_af.loc[df_load_af['MKT_TYPE']=='RTD', 'TIMESTAMP'].max(),
        df_load_bf.loc[df_load_bf['EXECUTION_TYPE']=='RTD', 'TIMESTAMP'].max(),
        df_ren_f.loc[df_ren_f['MARKET_RUN_ID']=='RTD', 'TIMESTAMP'].max(),
    ])
    df_results_rtd = pd.DataFrame(index=pd.date_range(start=t_start_rtd, end=t_end_rtd, freq='5min'))

    df_tmp = df_load_bf.loc[
        (df_load_bf['TAC_AREA_NAME']=='CA ISO-TAC')&(df_load_bf['EXECUTION_TYPE']=='RTD'), 
        ['MW', 'TIMESTAMP']
    ].set_index('TIMESTAMP')
    df_results_rtd.loc[:, 'LOAD_B_RTD'] = df_tmp['MW']

    df_tmp = df_load_af.loc[
        (df_load_af['BAA_GRP_ID']=='CA ISO-TAC')&(df_load_af['MKT_TYPE']=='RTD'), 
        ['MW', 'TIMESTAMP']
    ].set_index('TIMESTAMP')
    df_results_rtd.loc[:, 'LOAD_A_RTD'] = df_tmp['MW']

    for rh in df_ren_f['R_H'].unique():
        c = '_'.join([rh, 'RTD'])
        df_tmp = df_ren_f.loc[
            (df_ren_f['R_H']==rh)&(df_ren_f['MARKET_RUN_ID']=='RTD'),
            ['TIMESTAMP', 'MW']
        ].set_index('TIMESTAMP')
        df_results_rtd.loc[:, c] = df_tmp['MW']

    for ud in df_frrq['RAMP_TYPE'].unique():
        c = '_'.join([ud, 'RTD'])
        df_tmp = df_frrq.loc[
            (df_frrq['RAMP_TYPE']==ud)&(df_frrq['MARKET_RUN_ID']=='RTD')&(df_frrq['BAA_GRP_ID']=='CISO'),
            ['TIMESTAMP', 'UNCERTAINITY']
        ].set_index('TIMESTAMP')
        df_results_rtd.loc[:, c] = df_tmp['UNCERTAINITY']

    return df_results_rtd, df_results_rtpd

def process_actual(write_flag=False):
    '''
    This function is from previous process_raw_for_flexiramp function and only 
    keeps the part that deals with actual data, which is NOT from OASIS, but from
    CAISO's Today's Outlook. Not that I found some inconsistency with the OASIS
    data so this data is only for reference. This function is kept just for backup.
    '''
    # Load actual net load data
    ls_csv, ls_df = load_dir('./CAISO_OASIS/NetDemand.zip')
    for i in range(0, len(ls_csv)):
        fname = ls_csv[i]
        df    = ls_df[i]
        thisday = datetime.datetime.strptime(fname, 'CAISO-netdemand-%Y%m%d.csv')
        df = df.T.iloc[1:-1, :] # The first row is the header
        df.columns = ['HA FRCST', 'LOAD', 'NET LOAD']
        df['OPR_DT'] = thisday.strftime('%Y-%m-%d')
        df['YEAR']   = thisday.year
        df['MONTH']  = thisday.month
        df['DAY']    = thisday.day
        df['HOUR']     = pd.to_datetime(df.index).hour + 1 # CAISO OPR_HR starts from 1
        df['INTERVAL'] = pd.to_datetime(df.index).minute/5 + 1 # CAISO OPR_INTERVAL starts from 1
        ls_df[i] = df
    df_netload = pd.concat(ls_df, ignore_index=True)

    # Load actual solar and wind data, data availability: 20180410 to 20190531
    ls_csv, ls_df = load_dir('./CAISO_OASIS/Renewables.zip')
    for i in range(0, len(ls_csv)):
        fname = ls_csv[i]
        df    = ls_df[i]
        thisday = datetime.datetime.strptime(fname, 'CAISO-renewables-%Y%m%d.csv')
        df = df.T.iloc[1:-1, :] # The first row is the header
        df.columns = ['Solar', 'Wind', 'Geothermal', 'Biomass', 'Biogass', 'Small hydro', 'Batteries']
        df['OPR_DT'] = thisday.strftime('%Y-%m-%d')
        df['YEAR']   = thisday.year
        df['MONTH']  = thisday.month
        df['DAY']    = thisday.day
        df['HOUR']     = pd.to_datetime(df.index).hour + 1 # CAISO OPR_HR starts from 1
        df['INTERVAL'] = pd.to_datetime(df.index).minute/5 + 1 # CAISO OPR_INTERVAL starts from 1
        ls_df[i] = df
    df_actualver = pd.concat(ls_df, ignore_index=True)

    # Determine the starting time and ending time
    ls_dates = df_netload['OPR_DT'].unique()
    ls_dates.sort()

    # Prepare results container
    interval_perday = range(1, 13)*24
    hour_perday     = [i for i in range(1, 25) for j in range(1, 13)]
    dict_results = {'OPR_DT': [], 'OPR_HR': [], 'OPR_INTERVAL': []}
    for d in ls_dates:
        dict_results['OPR_DT'] += [d]*288
        dict_results['OPR_HR'] += hour_perday
        dict_results['OPR_INTERVAL'] += interval_perday
    df_results = pd.DataFrame(dict_results)
    df_results['TIME_STR'] = df_results['OPR_DT'] + '-' + df_results['OPR_HR'].astype('str') + '-' + df_results['OPR_INTERVAL'].astype('str')
    df_results = df_results.set_index('TIME_STR')

    # Add actual net load
    df_netload['TIME_STR'] = df_netload['OPR_DT'] + '-' +  df_netload['HOUR'].astype('str') + '-' + df_netload['INTERVAL'].astype(str)
    df_netload = df_netload.set_index('TIME_STR')
    df_results['NET_LOAD_ACTUAL'] = df_netload['NET LOAD']
    df_results['LOAD_ACTUAL']     = df_netload['LOAD']

    # Add actual renewables, we only care about solar and wind, note solar may 
    # also includes solar thermal.
    df_actualver['TIME_STR'] = df_actualver['OPR_DT'] + '-' +  df_actualver['HOUR'].astype('str') + '-' + df_actualver['INTERVAL'].astype(str)
    df_actualver = df_actualver.set_index('TIME_STR')
    df_results['SOLAR_ACTUAL'] = df_actualver['Solar']
    df_results['WIND_ACTUAL']  = df_actualver['Wind']

    if write_flag:
        df_results.to_csv('for_flexiramp_summary.csv', index=False)
    else:
        IP()

def process_raw_reg_req(write_flag=False):
    ''' 
    Process raw regulation requirements data from CAISO, DAM and RTM, output is 
    for_reg_reg_requirement.csv by default. 
    Note that DAM requirements are for every 1 hour, while RTM requirements are 
    for every 15 min.
    '''
    _, ls_df = load_dir('./CAISO_OASIS/AS_REQ.zip')
    df = pd.concat(ls_df, ignore_index=True)

    print 'Data loaded!'

    # Determine the starting time and ending time
    ls_dates = df['OPR_DT'].unique()
    ls_dates.sort()

    # Prepare results container
    interval_perday = range(1, 5)*24
    hour_perday     = [i for i in range(1, 25) for j in range(1, 5)]
    dict_results = {'OPR_DT': [], 'OPR_HR': [], 'OPR_INTERVAL': []}
    for d in ls_dates:
        dict_results['OPR_DT'] += [d]*24*4
        dict_results['OPR_HR'] += hour_perday
        dict_results['OPR_INTERVAL'] += interval_perday
    df_results = pd.DataFrame(dict_results)
    df_results.loc[:, 'TIME_STR'] = df_results['OPR_DT'] + '-' + df_results['OPR_HR'].astype('str') + '-' + df_results['OPR_INTERVAL'].astype('str')
    df_results = df_results.set_index('TIME_STR')

    # Assign regulation
    dict_btype = { # Bound type
        'MAX':   'AS Regional AS Requirement Maximum',
        'MIN':   'AS Regional AS Requirement Minimum',
    }
    for rtype in ['RU', 'RD']:
        for btype in ['MAX', 'MIN']:
            for mtype in ['RTM', 'DAM']: # Market type
                df_this = df.loc[
                    (df['MARKET_RUN_ID']==mtype) 
                    & (df['ANC_TYPE']==rtype) 
                    & (df['ANC_REGION']=='AS_CAISO_EXP') 
                    & (df['LABEL']==dict_btype[btype]), 
                    :
                ]
                df_this.loc[:, 'TIME_STR'] = df_this['OPR_DT'] +'-'+ df_this['OPR_HR'].astype('str') +'-'+ df_this['OPR_INTERVAL'].astype('str')
                df_this = df_this.set_index('TIME_STR')
                if mtype is 'DAM':
                    ls_tmp = list()
                    for i in range(1, 5):
                        df_tmp = df_this.copy()
                        df_tmp.loc[:, 'TIME_STR'] = df_tmp['OPR_DT'] +'-'+ df_tmp['OPR_HR'].astype('str') +'-'+ str(i)
                        df_tmp = df_tmp.set_index('TIME_STR')
                        ls_tmp.append(df_tmp)
                    df_this = pd.concat(ls_tmp)

                k = '_'.join([rtype, btype, mtype])
                df_results[k] = df_this['MW']

    if write_flag:
        df_results.to_csv('for_reg_reg_requirement.csv', index=False)
    else:
        IP()

    # df_results['RU_MAX_diff'] = df_results['RU_MAX_RTM'] - df_results['RU_MAX_DAM']
    # df_results['RU_MIN_diff'] = df_results['RU_MIN_RTM'] - df_results['RU_MIN_DAM']
    # df_results['RD_MAX_diff'] = df_results['RD_MAX_RTM'] - df_results['RD_MAX_DAM']
    # df_results['RD_MIN_diff'] = df_results['RD_MIN_RTM'] - df_results['RD_MIN_DAM']
    # df_results[['RU_MAX_diff', 'RU_MIN_diff', 'RD_MAX_diff', 'RD_MIN_diff']].plot()
    # plt.show()

def process_raw_reg_results(write_flag=False):
    ''' 
    Process raw regulation results data from CAISO, DAM and RTM, output is ...
    Note that DAM results are for every 1 hour, while RTM results are for every 
    15 min.
    '''
    _, ls_df = load_dir('./CAISO_OASIS/AS_RESULTS.zip')
    df = pd.concat(ls_df, ignore_index=True)

    print 'Data loaded!'

    # Determine the starting time and ending time
    ls_dates = df['OPR_DT'].unique()
    ls_dates.sort()

    # Prepare results container
    interval_perday = range(1, 5)*24 # The finiest interval is every 15 min (RTM)
    hour_perday     = [i for i in range(1, 25) for j in range(1, 5)]
    dict_results = {'OPR_DT': [], 'OPR_HR': [], 'OPR_INTERVAL': []}
    for d in ls_dates:
        dict_results['OPR_DT'] += [d]*24*4
        dict_results['OPR_HR'] += hour_perday
        dict_results['OPR_INTERVAL'] += interval_perday
    df_results = pd.DataFrame(dict_results)
    df_results.loc[:, 'TIME_STR'] = df_results['OPR_DT'] + '-' + df_results['OPR_HR'].astype('str') + '-' + df_results['OPR_INTERVAL'].astype('str')
    df_results = df_results.set_index('TIME_STR')

    # Assign regulation
    ls_as_region = ['AS_CAISO', 'AS_CAISO_EXP'] # Let's take a look at both regions
    for rtype in ['RU', 'RD']:
        for region in ls_as_region:
            for mtype in ['RTM', 'DAM']: # Market type
                df_this = df.loc[
                    (df['MARKET_RUN_ID']==mtype) 
                    & (df['ANC_TYPE']==rtype) 
                    & (df['ANC_REGION']==region)
                    & (df['RESULT_TYPE']=='AS_MW'), # This is the total procurement
                    :
                ]
                df_this.loc[:, 'TIME_STR'] = df_this['OPR_DT'] +'-'+ df_this['OPR_HR'].astype('str') +'-'+ df_this['OPR_INTERVAL'].astype('str')
                df_this = df_this.set_index('TIME_STR')
                if mtype is 'DAM':
                    # Because OPR_INTERVAL for all DAM products are 0, we have to add them manually
                    ls_tmp = list()
                    for i in range(1, 5):
                        df_tmp = df_this.copy()
                        df_tmp.loc[:, 'TIME_STR'] = df_tmp['OPR_DT'] +'-'+ df_tmp['OPR_HR'].astype('str') +'-'+ str(i)
                        df_tmp = df_tmp.set_index('TIME_STR')
                        ls_tmp.append(df_tmp)
                    df_this = pd.concat(ls_tmp)

                k = '_'.join([rtype, mtype, region])
                df_results[k] = df_this['MW']

    if write_flag:
        df_results.to_csv('for_reg_reg_results.csv', index=False)
    else:
        IP()

def process_raw_ace_error(write_flag=False, plot_flag=False):
    ''' 
    Process raw ACE data from CAISO, output file is for_reg_ACE.csv by default
    '''
    xls = r'C:\Users\bxl180002\Downloads\RampSolar\CAISO\HistoricalACE\2018 CAISO ACE Data Pull\2018 CAISO Area Control Error 1 minute average.xlsx'

    pd_xlsx = pd.ExcelFile(xls)
    ls_df = []
    for sheet in pd_xlsx.sheet_names:
        df = pd_xlsx.parse(sheet, header=None)
        df = df.drop([0, 1], axis=1) # The first two columns are always empty.
        df.columns = ['T', 'ACE']
        ls_df.append(df)
    print 'Data read!'

    # Remove invalid data
    df_ace = pd.concat(ls_df, ignore_index=True)
    df_ace = df_ace.replace(r'\s+', np.nan, regex=True) # Replace empty string with nan
    df_ace = df_ace.dropna() # Then, drop all rows with nan
    df_ace.loc[:, 'T'] = pd.to_datetime(df_ace['T']) # Why do we get an error if we do not include this line?

    # Daylight saving time handling, remove duplicates
    # The raw data is in local time, which uses daylight saving time.
    # The daylight saving time converts PST 2018-03-11 02:00 - 2018-11-04 00:59
    # into PDT 2018-03-11 03:00 - 2018-11-04 01:59

    # This part has no duplicates, so, easy to handle
    i_pdt = (
        ( df_ace['T'] >= datetime.datetime(2018, 3, 11, 3, 0) )
        & ( df_ace['T'] <= datetime.datetime(2018, 11, 4, 0, 59) )
    )
    df_ace.loc[i_pdt, 'T'] = df_ace.loc[i_pdt, 'T'] - pd.Timedelta(hours=1) 

    # This part has duplicates, for each pair of duplicates indexed by the 
    # same time stamp, the first one is in PDT time and the second one is in 
    # PST time.
    i_pdt = (
        ( df_ace['T'] >= datetime.datetime(2018, 11, 4, 1, 0) )
        & ( df_ace['T'] <= datetime.datetime(2018, 11, 4, 1, 59) )
    )
    m_updated = set()
    for i, _ in df_ace[i_pdt].iterrows():
        m = df_ace.loc[i, 'T'].minute
        if m not in m_updated:
            df_ace.loc[i, 'T'] = df_ace.loc[i, 'T'] - datetime.timedelta(hours=1)
            m_updated.add(m)

    # Add index identifier
    df_ace = df_ace.sort_values(by=['T'])
    df_ace['OPR_HR'] = df_ace['T'].dt.hour+1
    df_ace['OPR_INTERVAL'] = df_ace['T'].dt.minute+1
    df_ace['OPR_DT'] = df_ace['T'].dt.strftime('%Y-%m-%d')
    df_ace['TIME_STR'] = df_ace['OPR_DT'] + '-' + df_ace['OPR_HR'].astype('str') + '-' + df_ace['OPR_INTERVAL'].astype('str')

    # What's remaining as duplicated index correspond to identical ACE values, 
    # so can remove either one safely
    df_ace = df_ace.set_index('TIME_STR')
    df_ace = df_ace[~df_ace.index.duplicated(keep='first')]

    # We will only include days in 2018
    ls_d = df_ace.loc[df_ace['T'].dt.year==2018, 'T'].dt.strftime('%Y-%m-%d').unique().tolist()
    ls_d.sort()

    # Prepare result container
    interval_perday = range(1, 61)*24
    hour_perday     = [i for i in range(1, 25) for j in range(1, 61)]
    dict_results = {'OPR_DT': [], 'OPR_HR': [], 'OPR_INTERVAL': []}
    for d in ls_d:
        dict_results['OPR_DT'] += [d]*60*24
        dict_results['OPR_HR'] += hour_perday
        dict_results['OPR_INTERVAL'] += interval_perday
    df_results = pd.DataFrame(dict_results)
    df_results['TIME_STR'] = df_results['OPR_DT'] + '-' + df_results['OPR_HR'].astype('str') + '-' + df_results['OPR_INTERVAL'].astype('str')
    df_results = df_results.set_index('TIME_STR')

    # Now we can assign values to the result container.
    df_results['ACE'] = df_ace['ACE']
    df_results['T'] = df_ace['T']

    if write_flag:
        df_results.to_csv('for_reg_ACE.csv', index=False)
    else:
        IP()

    if plot_flag:
        # Plot time sequence for all months
        # for i in range(0, len(ls_df_new)):
        #     df = ls_df_new[i]
        #     t = pd_xlsx.sheet_names[i]
        #     plt.figure()
        #     plt.plot(df['ACE'])
        #     plt.title(t)

        # Plot distributions of every hour for each month
        for df in ls_df_new:
            for m in df['T'].dt.month.unique():
                if len(df.loc[df['T'].dt.month==m, 'ACE']) < 10: # Too few samples
                    continue
            plt.figure()
            for h in df['T'].dt.hour.unique():
                ace_counts, ace_edge = np.histogram(df.loc[df['T'].dt.hour==h, 'ACE'], bins= 50)
                # Numpy's "-1" is different from Python
                ace_center = (ace_edge[0: -1] + ace_edge[1:])/2 
                bin_width = np.mean(ace_edge[1:] - ace_edge[0:-1])
                plt.plot(
                    ace_center, 
                    ace_counts.astype('float64')/sum(ace_counts)/bin_width, 
                    'b', 
                    linewidth=0.5
                )
                plt.xlim([-1000, 1000])
                plt.ylabel('PDF')
                plt.xlabel('ACE (MW)')
                plt.title(m)

        # Plot time sequence area graph and 5% and 95% confidence levels.
        for df in ls_df_new:
            for m in df['T'].dt.month.unique():
                if len(df.loc[df['T'].dt.month==m, 'ACE']) < 10: # Too few samples
                    continue
                else:
                    plt.figure()
                    ax = plt.subplot(1, 1, 1)
                    for d in df['T'].dt.day.unique():
                        i_selected = (df['T'].dt.month==m)&(df['T'].dt.day==d)
                        x = -df.loc[i_selected, 'ACE'] # Note: Negative ACE
                        ax.fill_between(range(0, len(x)), 0, x, alpha=0.1, color='b')
                        # for j in range(0, 24):
                        #     plt.axvline(x=59+j*60, alpha=0.1, color='k', linewidth=0.5)
                    x = list()
                    y05 = list()
                    y95 = list()
                    for h in df['T'].dt.hour.unique():
                        i_selected = (df['T'].dt.month==m)&(df['T'].dt.hour==h)
                        ecdf = ECDF(df.loc[i_selected, 'ACE'])
                        p95 = np.interp(0.95, ecdf.y, ecdf.x)
                        p05 = np.interp(0.05, ecdf.y, ecdf.x)
                        x += range(h*60, h*60+60)
                        y95 += [-p95 for i in range(h*60, h*60+60)]
                        y05 += [-p05 for i in range(h*60, h*60+60)]
                    ax.plot(x, y95, 'r', linewidth=0.5)
                    ax.plot(x, y05, 'r', linewidth=0.5)
                    ax.set_ylim([-1000, 1000])
                    ax.set_xticks([60*i for i in range(0, 24)])
                    ax.set_xticklabels(range(0, 24))
                    plt.title(m)

        plt.show()

def prepare_data_for_mucun(write_flag=False):
    # RTD B and RTD actual
    df = load_and_process_netload_data(use_persistence=True)
    df.loc[:, 'TIME_STAMP'] = pd.to_datetime(
        df['OPR_DT'].astype(str) 
        + 
        ' ' 
        + 
        (df['OPR_HR']-1).map('{:02g}'.format) 
        + 
        ':'
        + 
        ((df['OPR_INTERVAL']-1)*5).map('{:02g}'.format), 
        format='%Y-%m-%d %H:%M'
    )
    df.loc[:, 'Year']   = df['TIME_STAMP'].dt.year
    df.loc[:, 'Month']  = df['TIME_STAMP'].dt.month
    df.loc[:, 'Day']    = df['TIME_STAMP'].dt.day
    df.loc[:, 'Hour']   = df['TIME_STAMP'].dt.hour
    df.loc[:, 'Minute'] = df['TIME_STAMP'].dt.minute
    if write_flag:
        df[['Year', 'Month', 'Day', 'Hour', 'Minute', 'NET_LOAD_ACTUAL', 'NET_LOAD_B_RTD']].to_csv('RTD_B.csv', index=False)
        df[['Year', 'Month', 'Day', 'Hour', 'Minute', 'NET_LOAD_ACTUAL', 'NET_LOAD_A_RTD']].to_csv('RTD_A.csv', index=False)

    # RTPD, just actual data, not forecast
    df = load_and_process_netload_data(use_persistence=False)
    df.loc[:, 'RTPD_INTERVAL'] = (df['OPR_INTERVAL']-1) - (df['OPR_INTERVAL']-1)%3
    df.loc[:, 'RTPD_INTERVAL'] = (df.loc[:, 'RTPD_INTERVAL']/3).astype(int) + 1
    df.loc[:, 'TIME_STR_RTPD'] = df['OPR_DT'].astype(str) + ' ' + (df['OPR_HR']-1).map('{:02g}'.format) + ':'+ ((df['RTPD_INTERVAL']-1)*15).map('{:02g}'.format)
    df_rtpd = df.groupby(by='TIME_STR_RTPD').mean()
    df_rtpd.loc[:, 'TIME_STAMP'] = pd.to_datetime(
        df_rtpd.index, 
        format='%Y-%m-%d %H:%M'
    )
    df_rtpd.loc[:, 'Year']   = df_rtpd['TIME_STAMP'].dt.year
    df_rtpd.loc[:, 'Month']  = df_rtpd['TIME_STAMP'].dt.month
    df_rtpd.loc[:, 'Day']    = df_rtpd['TIME_STAMP'].dt.day
    df_rtpd.loc[:, 'Hour']   = df_rtpd['TIME_STAMP'].dt.hour
    df_rtpd.loc[:, 'Minute'] = df_rtpd['TIME_STAMP'].dt.minute
    if write_flag:
        df_rtpd[['Year', 'Month', 'Day', 'Hour', 'Minute', 'NET_LOAD_ACTUAL', 'NET_LOAD_B_RTPD']].to_csv('RTPD_B.csv', index=False)


###########################################################
# Flexible ramp reserve requirement analysis starts here

def find_percentiles(df_data, ls_p):
    ''' Return percentiles per hour '''
    ls_strdays = df_data.columns.difference(['OPR_HR', 'OPR_INTERVAL'])
    NH = 24 # Number of hours in a day
    NI = len(df_data)/NH # Number of intervals in a hour, RTD: 12, RTPD: 4

    dict_percentiles = {p: list() for p in ls_p}
    for h in range(1, NH+1):
        samples = df_data.loc[df_data['OPR_HR']==h, ls_strdays].to_numpy().flatten()
        ecdf = ECDF(samples[~np.isnan(samples)])
        for p in dict_percentiles.iterkeys():
            dict_percentiles[p].append(np.interp(p, ecdf.y, ecdf.x))
    return dict_percentiles

def clear_sky_mean_pv(index_datetime, lat, lon, deltat):
    '''
    Calculate clear-sky mean PV output, index_datetime must be a pd.DatetimeIndex
    object, deltat is the number of minutes between two adjacent timestamp in 
    index_datetime
    '''

    # Expand the time sequence for mean
    deltatime = [np.timedelta64(-i, 'm') for i in range(deltat, 0, -1)]
    index_datetime_formean2d = np.repeat(np.expand_dims(index_datetime.values, axis=1), deltat, axis=1)
    deltatime_formean2d      = np.repeat(np.expand_dims(deltatime, axis = 0), index_datetime.shape[0], axis=0)
    index_datetime_formean = (index_datetime_formean2d - deltatime_formean2d).flatten()
    index_datetime_formean = pd.DatetimeIndex(index_datetime_formean)

    # get the module and inverter specifications from SAM
    sandia_modules = pvlib.pvsystem.retrieve_sam('SandiaMod')
    sapm_inverters = pvlib.pvsystem.retrieve_sam('cecinverter')

    # module = sandia_modules['Canadian_Solar_CS5P_220M___2009_']
    # inverter = sapm_inverters['ABB__MICRO_0_25_I_OUTD_US_208_208V__CEC_2014_']
    # array_ns = 1 # 1 serial units
    # array_np = 1 # 1 parallel strings

    # The following is Elina's module parameters
    module = sandia_modules['First_Solar_FS_272___2009_']
    inverter = sapm_inverters['Power_One__ULTRA_1500_TL_OUTD_1_US_690_x_y_z_690V__CEC_2013_']
    array_ns = 11 # 11 serial units
    array_np = 2360 # 2360 parallel strings

    # specify constant ambient air temp and wind for simplicity
    temp_air = 20
    wind_speed = 0

    system = {
        'module': module,
        'inverter': inverter,
        'surface_azimuth': 180,
        'surface_tilt': 0, # Rui's using 30, we use 0 follow Elina's matlab script
    }

    timezone = 'US/Pacific'
    altitude = 200 # For simplicity
    solpos = pvlib.solarposition.get_solarposition(index_datetime_formean, lat, lon)
    dni_extra = pvlib.irradiance.get_extra_radiation(index_datetime_formean)
    airmass = pvlib.atmosphere.get_relative_airmass(solpos['apparent_zenith'])
    pressure = pvlib.atmosphere.alt2pres(altitude)
    am_abs = pvlib.atmosphere.get_absolute_airmass(airmass, pressure)

    # Use clear-sky GHI for clear-sky power calculation
    loc = pvlib.location.Location(lat, lon, 'US/Pacific', altitude)
    ghi_cs = loc.get_clearsky(index_datetime_formean, model='ineichen')

    # Calculate AC power output
    tracker_data = pvlib.tracking.singleaxis(
        solpos['apparent_zenith'],
        solpos['azimuth'],
        axis_tilt=0, 
        axis_azimuth=180, 
        max_angle=90,
        backtrack=True, 
        gcr=2.0/7.0,
    )
    aoi = tracker_data['aoi']
    # aoi = pvlib.irradiance.aoi(
    #     system['surface_tilt'], 
    #     system['surface_azimuth'],
    #     solpos['apparent_zenith'].values, 
    #     solpos['azimuth'].values,
    # )
    total_irrad = pvlib.irradiance.get_total_irradiance(
        # system['surface_tilt'],
        # system['surface_azimuth'],
        tracker_data['surface_tilt'],
        tracker_data['surface_azimuth'],
        solpos['apparent_zenith'].values,
        solpos['azimuth'].values,
        ghi_cs['dni'].values, 
        ghi_cs['ghi'].values, 
        ghi_cs['dhi'].values,
        dni_extra=dni_extra.values,
        model='haydavies',
    )
    temps = pvlib.pvsystem.sapm_celltemp(
        total_irrad['poa_global'],
        wind_speed, 
        temp_air,
    )
    effective_irradiance = pvlib.pvsystem.sapm_effective_irradiance(
        total_irrad['poa_direct'], 
        total_irrad['poa_diffuse'],
        am_abs.values, 
        aoi, 
        module
    )
    dc = pvlib.pvsystem.sapm(effective_irradiance, temps['temp_cell'].values, module)
    ac = pvlib.pvsystem.snlinverter(dc['v_mp']*array_ns, dc['p_mp']*array_ns*array_np, inverter)

    ac[np.isnan(ac)] = 0
    power_cs_formean = ac
    power_cs = power_cs_formean.reshape(index_datetime.shape[0], deltat).mean(axis=1)
    power_cs[power_cs<=0] = 0
    return power_cs

def load_and_process_netload_data(use_persistence=False):
    '''
    This function load for_flexiramp_summary.csv and calculates net load and net load 
    forecast errors
    '''
    import pvlib
    df = pd.read_csv('for_flexiramp_summary.csv')

    # We have all forecast for binding intervals, so just use it
    df['Solar_NP15_B_RTPD'] = df['Solar_NP15_RTPD']
    df['Solar_SP15_B_RTPD'] = df['Solar_SP15_RTPD']
    df['Solar_ZP26_B_RTPD'] = df['Solar_ZP26_RTPD']
    df['Wind_NP15_B_RTPD']  = df['Wind_NP15_RTPD']
    df['Wind_SP15_B_RTPD']  = df['Wind_SP15_RTPD']

    df['Solar_NP15_B_RTD'] = df['Solar_NP15_RTD']
    df['Solar_SP15_B_RTD'] = df['Solar_SP15_RTD']
    df['Solar_ZP26_B_RTD'] = df['Solar_ZP26_RTD']
    df['Wind_NP15_B_RTD']  = df['Wind_NP15_RTD']
    df['Wind_SP15_B_RTD']  = df['Wind_SP15_RTD']

    if not use_persistence:
        df['Solar_NP15_A_RTD'] = df['Solar_NP15_RTD']
        df['Solar_SP15_A_RTD'] = df['Solar_SP15_RTD']
        df['Solar_ZP26_A_RTD'] = df['Solar_ZP26_RTD']
        df['Wind_NP15_A_RTD']  = df['Wind_NP15_RTD']
        df['Wind_SP15_A_RTD']  = df['Wind_SP15_RTD']

        df['Solar_NP15_A_RTPD'] = df['Solar_NP15_RTPD']
        df['Solar_SP15_A_RTPD'] = df['Solar_SP15_RTPD']
        df['Solar_ZP26_A_RTPD'] = df['Solar_ZP26_RTPD']
        df['Wind_NP15_A_RTPD']  = df['Wind_NP15_RTPD']
        df['Wind_SP15_A_RTPD']  = df['Wind_SP15_RTPD']
    else:
        # RTD solar
        ########################################################################
        # We use persistence forecast, wind is shifted by 1 interval, and 
        # solar from Cong's forecast
        df.loc[:, 'timestamp'] = pd.to_datetime(
            df['OPR_DT'].astype(str) 
            + 
            ' ' 
            + 
            (df['OPR_HR']-1).map('{:02g}'.format) 
            + 
            ':'
            + 
            ((df['OPR_INTERVAL']-1)*5).map('{:02g}'.format), 
            format='%Y-%m-%d %H:%M'
        ).dt.tz_localize('US/Pacific', nonexistent='NaT', ambiguous='NaT')

        loc = dict()
        loc['NP15'] = pvlib.location.Location(38.548, -121.531, 'US/Pacific',4)
        loc['ZP26'] = pvlib.location.Location(35.853, -120.212, 'US/Pacific',446)
        loc['SP15'] = pvlib.location.Location(33.907, -116.257, 'US/Pacific',753)

        # # Use clear-sky GHI as full output
        # dict_df_csghi = dict()
        # for k in loc.keys():
        #     df_csghi = loc[k].get_clearsky(pd.DatetimeIndex(df['timestamp']), model='ineichen')
        #     # Due to daylight saving time conversion, we'll get NaT at the 
        #     # transitioning time, however, we know for sure GHI = 0 at that time
        #     df_csghi.loc['NaT', :] = 0
        #     dict_df_csghi[k] = df_csghi

        # # Use persistence clear-sky index to get advisory solar power forecast
        # for k in loc.keys():
        #     cs_ghi = dict_df_csghi[k]['ghi'].values
        #     power  = df['Solar_'+k+'_B_RTD'].values
        #     csindex = np.zeros(power.shape)
        #     csindex[~(cs_ghi==0)] = power[~(cs_ghi==0)]/cs_ghi[~(cs_ghi==0)]
        #     csindex[np.isnan(csindex)] = 0
        #     csindex[np.isinf(csindex)] = 0
        #     frcst_adv = np.concatenate(
        #         [[np.nan], dict_df_csghi[k]['ghi'].values[1:]*csindex[0:-1]]
        #     )
        #     df.loc[:, 'Solar_'+k+'_A_RTD'] = frcst_adv

        #     df.loc[:, 'AB_diff_'+k] = np.abs(df['Solar_'+k+'_A_RTD'] - df['Solar_'+k+'_B_RTD'])/df['Solar_'+k+'_B_RTD']

        # Use clear-sky power as full output
        # Use persistence clear-sky index to get advisory solar power forecast
        index_datetime = pd.DatetimeIndex(df['timestamp'])
        for k in loc.keys():
            power_cs = clear_sky_mean_pv(index_datetime, loc[k].latitude, loc[k].longitude, 5)
            power  = df['Solar_'+k+'_B_RTD'].values
            csindex = np.zeros(power.shape)
            csindex[~(power_cs==0)] = power[~(power_cs==0)]/power_cs[~(power_cs==0)]
            csindex[np.isnan(csindex)] = 0
            csindex[np.isinf(csindex)] = 0
            frcst_adv = np.concatenate(
                [[np.nan], power_cs[1:]*csindex[0:-1]]
            )
            df.loc[:, 'Solar_'+k+'_A_RTD'] = frcst_adv

            df.loc[:, 'AB_diff_'+k] = np.abs(df['Solar_'+k+'_A_RTD'] - df['Solar_'+k+'_B_RTD'])/df['Solar_'+k+'_B_RTD']

            # Fix numeric error when CS GHI is small and adv forecast is too large.
            # Basically, whenever B-A forecast percentage error > 100%, use binding 
            # forecast.
            df.loc[df['AB_diff_'+k]>1, 'Solar_'+k+'_A_RTD'] = df.loc[df['AB_diff_'+k]>1, 'Solar_'+k+'_B_RTD']

            # # Uncomment the following part to visualize
            # fig = plt.figure()
            # ax = plt.subplot(111)
            # ax1 = ax.twinx()
            # df[['timestamp', 'Solar_'+k+'_A_RTD', 'Solar_'+k+'_B_RTD']].plot(ax=ax, x='timestamp')
            # df[['timestamp', 'AB_diff_'+k]].plot(ax=ax1, x='timestamp', style='r-')

            plt.figure()
            ax1 = plt.subplot(111)
            ax2 = ax1.twinx()
            df[['timestamp', 'Solar_'+k+'_B_RTD', 'Solar_'+k+'_A_RTD']].plot(ax = ax1, x='timestamp')
            # dict_df_csghi[k]['ghi'].plot(ax=ax2,x='timestamp', color='green')
            ax2.plot(df['timestamp'], power_cs, color='green')
            ax1.set_ylabel('Power (MW)')
            ax2.set_ylabel('Clear-sky GHI (W/m^2)')

        # RTD wind
        ########################################################################
        df['Wind_NP15_A_RTD']  = np.concatenate([ [np.nan], df['Wind_NP15_RTD'].values[0:-1] ])
        df['Wind_SP15_A_RTD']  = np.concatenate([ [np.nan], df['Wind_SP15_RTD'].values[0:-1] ])

        # RTPD solar and wind, just use binding forecast as advisory forecat 
        # because the binding RTPD forecasts was never used in flexiramp calculation
        ########################################################################
        df['Solar_NP15_A_RTPD'] = df['Solar_NP15_RTPD']
        df['Solar_SP15_A_RTPD'] = df['Solar_SP15_RTPD']
        df['Solar_ZP26_A_RTPD'] = df['Solar_ZP26_RTPD']
        df['Wind_NP15_A_RTPD']  = df['Wind_NP15_RTPD']
        df['Wind_SP15_A_RTPD']  = df['Wind_SP15_RTPD']

        # df_rtd_a  = pd.read_csv('C:/Users/bxl180002/Downloads/RampSolar/CAISO/RTD_Forecasts.csv')
        # df_rtpd_a = pd.read_csv('C:/Users/bxl180002/Downloads/RampSolar/CAISO/RTPD_Forecasts.csv')
        # df_rtd_a  = pd.read_csv('~/Downloads/Tmp/RTD_Forecasts.csv')
        # df_rtpd_a = pd.read_csv('~/Downloads/Tmp/RTPD_Forecasts.csv')
        # df_rtd_a['Wind_A_NP15_RTD']   = np.concatenate([ [np.nan], df['Wind_NP15_RTD'].values[0:-1] ])
        # df_rtd_a['Wind_A_SP15_RTD']   = np.concatenate([ [np.nan], df['Wind_SP15_RTD'].values[0:-1] ])
        # df_rtpd_a['Wind_A_NP15_RTPD'] = np.concatenate([ [np.nan], df.loc[df['OPR_INTERVAL']%3 == 0, 'Wind_NP15_RTPD'].values[0:-1] ])
        # df_rtpd_a['Wind_A_SP15_RTPD'] = np.concatenate([ [np.nan], df.loc[df['OPR_INTERVAL']%3 == 0, 'Wind_SP15_RTPD'].values[0:-1] ])

        # df['Solar_NP15_A_RTD'] = df_rtd_a['Forc_Solar_NP15_RTD'].values
        # df['Solar_SP15_A_RTD'] = df_rtd_a['Forc_Solar_SP15_RTD'].values
        # df['Solar_ZP26_A_RTD'] = df_rtd_a['Forc_Solar_ZP26_RTD'].values
        # df['Wind_NP15_A_RTD']  = df_rtd_a['Wind_A_NP15_RTD'].values
        # df['Wind_SP15_A_RTD']  = df_rtd_a['Wind_A_SP15_RTD'].values

        # df['Solar_NP15_A_RTPD'] = df_rtpd_a['Forc_Solar_NP15_RTPD'].values.repeat(3)
        # df['Solar_SP15_A_RTPD'] = df_rtpd_a['Forc_Solar_SP15_RTPD'].values.repeat(3)
        # df['Solar_ZP26_A_RTPD'] = df_rtpd_a['Forc_Solar_ZP26_RTPD'].values.repeat(3)
        # df['Wind_NP15_A_RTPD']  = df_rtpd_a['Wind_A_NP15_RTPD'].values.repeat(3)
        # df['Wind_SP15_A_RTPD']  = df_rtpd_a['Wind_A_SP15_RTPD'].values.repeat(3)

        # # Fix spikes due to persistence forecast.
        # threshold = 600.0/100 # Percentage threshold: 100%
        # columns_with_spikes = [
        #     ('Solar_NP15_B_RTD', 'Solar_NP15_A_RTD'),
        #     ('Solar_SP15_B_RTD', 'Solar_SP15_A_RTD'),
        #     ('Solar_ZP26_B_RTD', 'Solar_ZP26_A_RTD'),
        #     ('Solar_NP15_B_RTPD', 'Solar_NP15_A_RTPD'),
        #     ('Solar_SP15_B_RTPD', 'Solar_SP15_A_RTPD'),
        #     ('Solar_ZP26_B_RTPD', 'Solar_ZP26_A_RTPD'),
        # ]
        # print 'Summary of persistence forecast of solar power, thredhold = {:>.2f}%'.format(threshold*100)
        # for c_b, c_a in columns_with_spikes:
        #     # Column of binding forecast, column of advisory forecast
        #     tmp = (df[c_a] - df[c_b])/df[c_b] # Relative errors
        #     index_spikes = tmp[tmp>threshold].index
        #     # df.loc[index_spikes, c_a] = (df.loc[index_spikes-1, c_a].values + df.loc[index_spikes+1, c_a].values)/2
        #     df.loc[index_spikes, c_a] = df.loc[index_spikes, c_b]
        #     print '{:<20s} has {:>4g} spikes out of {:>7g} points. Fraction: {:>4.3}%'.format(
        #         c_a,
        #         index_spikes.size,
        #         df.index.size,
        #         float(index_spikes.size)/df.index.size*100,
        #     )
        # IP()

    # Now we can calculate the net load for both the advisory and binding intervals
    df['NET_LOAD_B_RTD'] = (
        df['LOAD_B_RTD'] - 
        df['Wind_NP15_B_RTD'] -
        df['Wind_SP15_B_RTD'] - 
        df['Solar_NP15_B_RTD'] - 
        df['Solar_SP15_B_RTD'] - 
        df['Solar_ZP26_B_RTD']
    )
    df['NET_LOAD_B_RTPD'] = (
        df['LOAD_B_RTPD'] -
        df['Wind_NP15_B_RTPD'] -
        df['Wind_SP15_B_RTPD'] -
        df['Solar_NP15_B_RTPD'] -
        df['Solar_SP15_B_RTPD'] -
        df['Solar_ZP26_B_RTPD']
    )

    df['NET_LOAD_A_RTD'] = (
        # df['LOAD_A_RTD'] - 
        df['LOAD_B_RTD'] - # To remove gross load ramp
        df['Wind_NP15_A_RTD'] -
        df['Wind_SP15_A_RTD'] - 
        df['Solar_NP15_A_RTD'] - 
        df['Solar_SP15_A_RTD'] - 
        df['Solar_ZP26_A_RTD']
    )
    df['NET_LOAD_A_RTPD'] = (
        df['LOAD_A_RTPD'] -
        df['Wind_NP15_A_RTPD'] -
        df['Wind_SP15_A_RTPD'] -
        df['Solar_NP15_A_RTPD'] -
        df['Solar_SP15_A_RTPD'] -
        df['Solar_ZP26_A_RTPD']
    )


    # Now we calculate forecast errors
    df['FORECAST_ERROR_act_Brtd'] = df['NET_LOAD_ACTUAL'] - df['NET_LOAD_B_RTD']
    df['FORECAST_ERROR_Brtd_Artd'] = df['NET_LOAD_B_RTD'] - df['NET_LOAD_A_RTD']
    df['FORECAST_ERROR_Brtd_Artpd'] = df['NET_LOAD_B_RTD'] - df['NET_LOAD_A_RTPD']

    # dict_legend = {
    #     'FORECAST_ERROR_act_Brtd': 'Actual - RTD B',
    #     'FORECAST_ERROR_Brtd_Artd': 'RTD B - RTD A', 
    #     'FORECAST_ERROR_Brtd_Artpd': 'RTD B - RTPD A',
    # }

    # A simple distributions
    # ls = list()
    # for c in df.columns:
    #     if c.startswith('FORECAST_ERROR'):
    #         ls.append(c)
    #         sample = df.loc[~pd.isna(df[c]), c].to_numpy()
    #         counts, edge = np.histogram(sample, bins=100)
    #         # Numpy's "-1" is different from Python
    #         center = (edge[0: -1] + edge[1:])/2 
    #         bin_width = np.mean(edge[1:] - edge[0:-1])
    #         plt.plot(
    #             center, 
    #             counts.astype('float64')/sum(counts)/bin_width,
    #         )
    #         plt.legend([dict_legend[i] for i in ls])

    return df.drop(
        [
            'Solar_NP15_RTD',
            'Solar_SP15_RTD',
            'Solar_ZP26_RTD',
            'Wind_NP15_RTD',
            'Wind_SP15_RTD',
            'Solar_NP15_RTPD',
            'Solar_SP15_RTPD',
            'Solar_ZP26_RTPD',
            'Wind_NP15_RTPD',
            'Wind_SP15_RTPD',
        ],
        axis=1,
    )

def tmp_plot_forecast():
    # This function is to give an overview to Cong's persistence forecast
    df = load_and_process_netload_data(use_persistence=True)
    df['TIME_STAMP'] = pd.to_datetime(
        df['OPR_DT']+'-'+(df['OPR_HR']-1).map('{:02g}'.format)+'-'+(5*df['OPR_INTERVAL']-5).map('{:02g}'.format),
        format='%Y-%m-%d-%H-%M'
    )
    df[['TIME_STAMP', 'Solar_NP15_B_RTD', 'Solar_NP15_A_RTD']].plot(title='solar NP15 RTD', x='TIME_STAMP')
    df[['TIME_STAMP', 'Solar_SP15_B_RTD', 'Solar_SP15_A_RTD']].plot(title='solar SP15 RTD', x='TIME_STAMP')
    df[['TIME_STAMP', 'Solar_ZP26_B_RTD', 'Solar_ZP26_A_RTD']].plot(title='solar ZP15 RTD', x='TIME_STAMP')
    df[['TIME_STAMP', 'Solar_NP15_B_RTPD', 'Solar_NP15_A_RTPD']].plot(title='solar NP15 RTPD', x='TIME_STAMP')
    df[['TIME_STAMP', 'Solar_SP15_B_RTPD', 'Solar_SP15_A_RTPD']].plot(title='solar SP15 RTPD', x='TIME_STAMP')
    df[['TIME_STAMP', 'Solar_ZP26_B_RTPD', 'Solar_ZP26_A_RTPD']].plot(title='solar ZP15 RTPD', x='TIME_STAMP')
    df[['TIME_STAMP', 'Wind_NP15_B_RTD', 'Wind_NP15_A_RTD']].plot(title='wind NP15 RTD', x='TIME_STAMP')
    df[['TIME_STAMP', 'Wind_SP15_B_RTD', 'Wind_SP15_A_RTD']].plot(title='wind SP15 RTD', x='TIME_STAMP')
    df[['TIME_STAMP', 'Wind_NP15_B_RTPD', 'Wind_NP15_A_RTPD']].plot(title='wind NP15 RTPD', x='TIME_STAMP')
    df[['TIME_STAMP', 'Wind_SP15_B_RTPD', 'Wind_SP15_A_RTPD']].plot(title='wind SP15 RTPD', x='TIME_STAMP')
    plt.show()

def baseline_flexiramp_for_day(YYYY, MM, DD, use_persistence=False):
    '''
    Make graphs for day YYYY/MM/DD
    '''

    def plot_ramp_stats(ax, df_data, direction=None):
        ls_strdays = df_data.columns.difference(['OPR_HR', 'OPR_INTERVAL'])
        NH = 24 # Number of hours in a day
        NI = len(df_data)/NH # Number of intervals in a hour, RTD: 12, RTPD: 4

        # Superimposed bar plots
        # for strd in ls_strdays:
        #     ax.bar(df_data.index, df_data[strd].tolist(), alpha=0.2, color='b')

        # Uncertainty components
        dict_percentiles = find_percentiles(df_data, [0.025, 0.975])
        if not direction: # Both directions, RTD
            for k in dict_percentiles:
                ax.plot(
                    df_data.index, 
                    np.array(dict_percentiles[k]).repeat(NI), 
                    'r', 
                    linewidth=1
                )
        elif direction is 'U': # RTPD FRU, max RTD B - RTPD A is used
            ax.plot(
                df_data.index, 
                np.array(dict_percentiles[0.975]).repeat(NI), 
                'r', 
                linewidth=1
            )
        elif direction is 'D': # RTPD FRD, min RTD B - RTPD A is used
            ax.plot(
                df_data.index, 
                np.array(dict_percentiles[0.025]).repeat(NI), 
                'r', 
                linewidth=1
            )


    ############################################################################
    # Example starts here.
    ############################################################################
    df = load_and_process_netload_data(use_persistence=use_persistence)

    # Let's try to get the uncertainty requirements for a weekday and a weekend
    D_TARGET = datetime.datetime(YYYY, MM, DD)
    ls_hist_days = list()

    # Collect days of 30 weekdays and weekends before the target day
    d = D_TARGET
    if D_TARGET.weekday() < 5: # This is a weekday
        while len(ls_hist_days) < 30:
            d = d - datetime.timedelta(days=1)
            if d.weekday() < 5: # Python starts from 0, so 5 is Saturday
                ls_hist_days.append(d)
    else: # This is a weekend
        while len(ls_hist_days) < 30:
            d = d - datetime.timedelta(days=1)
            if d.weekday() >= 5: # Python starts from 0, so 5 is Saturday
                ls_hist_days.append(d)

    # Collect the RTD net load forecast errors for the 30 weekdays and weekends
    df_rtd = pd.DataFrame(index = range(288), columns=['OPR_HR', 'OPR_INTERVAL'])
    for d in ls_hist_days:
        str_d = d.strftime('%Y-%m-%d')
        df_rtd.loc[:, str_d] = df.loc[df['OPR_DT']==str_d, 'FORECAST_ERROR_Brtd_Artd'].tolist()
    df_rtd['OPR_HR'] = df.loc[df['OPR_DT']==str_d, 'OPR_HR'].tolist()
    df_rtd['OPR_INTERVAL'] = df.loc[df['OPR_DT']==str_d, 'OPR_INTERVAL'].tolist()

    # Collect the RTPD net load forecast errors for the 30 weekdays and weekends
    df_max = pd.DataFrame(index = range(96), columns=['OPR_HR', 'OPR_INTERVAL'])
    df_min = pd.DataFrame(index = range(96), columns=['OPR_HR', 'OPR_INTERVAL'])
    for d in ls_hist_days:
        str_d = d.strftime('%Y-%m-%d')
        ls_tmp = list()
        for i in range(3):
            ls_tmp.append(
                df.loc[
                    (df['OPR_DT']==str_d)&(df['OPR_INTERVAL']%3==i), 
                    'FORECAST_ERROR_Brtd_Artpd'
                ].tolist()
            )
        df_max.loc[:, str_d] = np.array(ls_tmp).T.max(axis=1)
        df_min.loc[:, str_d] = np.array(ls_tmp).T.min(axis=1)
    df_max['OPR_HR'] = df.loc[(df['OPR_DT']==str_d)&(df['OPR_INTERVAL']%3==0), 'OPR_HR'].tolist()
    df_min['OPR_HR'] = df_max['OPR_HR'].tolist()
    df_max['OPR_INTERVAL'] = (df.loc[(df['OPR_DT']==str_d)&(df['OPR_INTERVAL']%3==0), 'OPR_INTERVAL']/3).astype('int').tolist()
    df_min['OPR_INTERVAL'] = df_max['OPR_INTERVAL'].tolist()

    # Try to make a bar plot for RTD weekdays
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    plot_ramp_stats(ax, df_rtd)
    ax.set_title('RTD, up and down, ' + D_TARGET.strftime('%m/%d/%y, %a'))
    ax.set_ylabel('MW')
    ax.set_xlabel('Hour')
    df_actual = df.loc[
        (df['OPR_DT']==D_TARGET.strftime('%Y-%m-%d'))&(df['OPR_INTERVAL']==1), 
        ['UP_RTD', 'DOWN_RTD']
    ].reset_index(drop=True) # Actual up and down ramp reserve requirements
    ax.plot(
        range(0, 288),
        df_actual['UP_RTD'].repeat(12), 
        'b',
        linewidth=1
    )
    ax.plot(
        range(0, 288),
        -df_actual['DOWN_RTD'].repeat(12), 
        'b',
        linewidth=1
    )
    ax.set_xticks(range(0, 288,12))
    ax.set_xticklabels([str(i/12) for i in (range(0, 288, 12))])
    h_b  = mlines.Line2D([], [], color='r', linestyle = '-') # Baseline
    h_o = mlines.Line2D([], [], color='b', linestyle = '-') # Oasis
    ax.legend([h_b, h_o], ['Baseline', 'OASIS'])

    # Try to make a bar plot for RTPD market
    df_actual = df.loc[
        (df['OPR_DT']==D_TARGET.strftime('%Y-%m-%d'))&(df['OPR_INTERVAL']==1), 
        ['UP_RTPD', 'DOWN_RTPD']
    ].reset_index(drop=True) # Actual up and down ramp reserve requirements
    plt.figure()
    gs = gridspec.GridSpec(
        2, 1,
        height_ratios=[1]*2,
        width_ratios=[1]*1,
        hspace=0.5,
        wspace=0,
    )
    # FRU
    ax = plt.subplot( gs[0, 0])
    plot_ramp_stats(ax, df_max, 'U') # FRU
    ax.set_title('RTPD, up, ' + D_TARGET.strftime('%m/%d/%y, %a'))
    # ax.set_title('2/22/19, Fri.')
    # ax.annotate(
    #     'RTPD, up',
    #     xy=(0.1, 0.9), 
    #     xycoords='axes fraction',
    #     xytext=(.1, .9), 
    #     textcoords='axes fraction',
    #     # fontsize=12,
    # )
    ax.set_ylabel('MW')
    ax.set_xlabel('Hour')
    ax.plot(
        range(0, 96),
        df_actual['UP_RTPD'].repeat(4), 
        'b',
        linewidth=1
    )
    ax.set_xticks(range(0, 96,4))
    ax.set_xticklabels([str(i/4) for i in (range(0, 96,4))])
    h_b  = mlines.Line2D([], [], color='r', linestyle = '-') # Baseline
    h_o = mlines.Line2D([], [], color='b', linestyle = '-') # Oasis
    ax.legend([h_b, h_o], ['Baseline', 'OASIS'])

    # FRD
    ax = plt.subplot( gs[1, 0])
    plot_ramp_stats(ax, df_min, 'D') # FRD
    ax.set_title('RTPD, down, ' + D_TARGET.strftime('%m/%d/%y, %a'))
    # ax.annotate(
    #     'RTPD, down',
    #     xy=(0.1, 0.9), 
    #     xycoords='axes fraction',
    #     xytext=(.1, .9), 
    #     textcoords='axes fraction',
    #     # fontsize=12,
    # )
    ax.set_ylabel('MW')
    ax.set_xlabel('Hour')
    ax.plot(
        range(0, 96),
        -df_actual['DOWN_RTPD'].repeat(4), 
        'b',
        linewidth=1
    )
    ax.set_xticks(range(0, 96,4))
    ax.set_xticklabels([str(i/4) for i in (range(0, 96,4))])
    h_b  = mlines.Line2D([], [], color='r', linestyle = '-') # Baseline
    h_o = mlines.Line2D([], [], color='b', linestyle = '-') # Oasis
    ax.legend([h_b, h_o], ['Baseline', 'OASIS'])

    # Load and net load of 30 days, use actual net load and forecasted RTD B load
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    for d in ls_hist_days:
        h_nl = plt.plot(
            df_rtd.index,
            df.loc[df['OPR_DT']==d.strftime('%Y-%m-%d'), 'NET_LOAD_ACTUAL'].tolist(),
            color='b', 
            linewidth=0.5,
        )
        h_l  = plt.plot(
            df_rtd.index,
            df.loc[df['OPR_DT']==d.strftime('%Y-%m-%d'), 'LOAD_B_RTD'].tolist(),
            color='k',
            linewidth=0.5
        )
    ax.set_title(
        'CAISO load and net load, 30 ' +
        ('weekdays ' if D_TARGET.weekday() < 5 else 'weekends ') + 
        'before ' + 
        D_TARGET.strftime('%m/%d/%y, %a')
    )
    ax.set_ylabel('MW')
    ax.set_xlabel('Hour')
    ax.set_xticks(range(0, 288,12))
    ax.set_xticklabels([str(i/12) for i in (range(0, 288, 12))])
    h_l  = mlines.Line2D([], [], color='k', linestyle = '-')
    h_nl = mlines.Line2D([], [], color='b', linestyle = '-')
    ax.legend([h_nl, h_l], ['Net load', 'Load'])

    # Plot distributions of selected hours
    ls_h = [1, 6, 12, 18]
    ls_color = ['b', 'r', 'g', 'k']
    plt.figure()
    for h in ls_h:
        df_tmp = df_rtd
        ls_strdays = df_tmp.columns.difference(['OPR_HR', 'OPR_INTERVAL'])
        samples = df_tmp.loc[df_tmp['OPR_HR']==h, ls_strdays].to_numpy().flatten()
        counts, edge = np.histogram(samples[~np.isnan(samples)], bins=25)
        # Numpy's "-1" is different from Python
        center = (edge[0: -1] + edge[1:])/2 
        bin_width = np.mean(edge[1:] - edge[0:-1])
        plt.plot(
            center, 
            counts.astype('float64')/sum(counts)/bin_width,
            color=ls_color[ls_h.index(h)],
            linewidth = 0.5,
        )
    plt.legend([str(i) for i in ls_h])

    plt.figure()
    for h in ls_h:
        df_tmp = df_max
        ls_strdays = df_tmp.columns.difference(['OPR_HR', 'OPR_INTERVAL'])
        samples = df_tmp.loc[df_tmp['OPR_HR']==h, ls_strdays].to_numpy().flatten()
        counts, edge = np.histogram(samples[~np.isnan(samples)], bins=25)
        # Numpy's "-1" is different from Python
        center = (edge[0: -1] + edge[1:])/2 
        bin_width = np.mean(edge[1:] - edge[0:-1])
        plt.plot(
            center, 
            counts.astype('float64')/sum(counts)/bin_width,
            color=ls_color[ls_h.index(h)],
            linewidth = 0.5,
        )
    plt.legend([str(i) for i in ls_h])

    plt.figure()
    for h in ls_h:
        df_tmp = df_min
        ls_strdays = df_tmp.columns.difference(['OPR_HR', 'OPR_INTERVAL'])
        samples = df_tmp.loc[df_tmp['OPR_HR']==h, ls_strdays].to_numpy().flatten()
        counts, edge = np.histogram(samples[~np.isnan(samples)], bins=25)
        # Numpy's "-1" is different from Python
        center = (edge[0: -1] + edge[1:])/2 
        bin_width = np.mean(edge[1:] - edge[0:-1])
        plt.plot(
            center, 
            counts.astype('float64')/sum(counts)/bin_width,
            color=ls_color[ls_h.index(h)],
            linewidth = 0.5,
        )
    plt.legend([str(i) for i in ls_h])

    # plt.figure()
    # for h in range(1, 25):
    #     samples = df_wkend.loc[df_wkday['OPR_HR']==h, df_wkend.columns.difference(['OPR_HR', 'OPR_INTERVAL'])].to_numpy().flatten()
    #     counts, edge = np.histogram(samples[~np.isnan(samples)], bins=100)
    #     # Numpy's "-1" is different from Python
    #     center = (edge[0: -1] + edge[1:])/2 
    #     bin_width = np.mean(edge[1:] - edge[0:-1])
    #     plt.plot(
    #         center, 
    #         counts.astype('float64')/sum(counts)/bin_width,
    #         color='b',
    #         linewidth = 0.5,
    #         alpha=0.5
    #     )

    # IP()
    plt.show()

def baseline_flexiramp():
    df = load_and_process_netload_data(use_persistence=True)
    ############################################################################
    # Derive FRU and FRD for all weekdays and weekends
    ############################################################################
    df['WEEKDAY'] = pd.to_datetime(df['OPR_DT']).dt.weekday
    dict_ls_days  = dict()
    dict_ls_days['WEEKDAY'] = df.loc[df['WEEKDAY']<=4, 'OPR_DT'].unique().tolist()
    dict_ls_days['WEEKEND'] = df.loc[df['WEEKDAY']>4,  'OPR_DT'].unique().tolist()
    str_dend = df.loc[:, 'OPR_DT'].unique()[-1]

    # Collect net load forecast errors for RTD markets
    dict_df_rtd = dict()
    for w in dict_ls_days.iterkeys():
        dict_df_rtd[w] = pd.DataFrame(index = range(288), columns=['OPR_HR', 'OPR_INTERVAL'])
        for d in dict_ls_days[w]:
            dict_df_rtd[w].loc[:, d] = df.loc[df['OPR_DT']==d, 'FORECAST_ERROR_Brtd_Artd'].tolist()
        dict_df_rtd[w]['OPR_HR'] = df.loc[df['OPR_DT']==d, 'OPR_HR'].tolist()
        dict_df_rtd[w]['OPR_INTERVAL'] = df.loc[df['OPR_DT']==d, 'OPR_INTERVAL'].tolist()

    # Collect net load forecast errors for RTPD markets: max and min 
    dict_df_rtpd_max = dict()
    dict_df_rtpd_min = dict()
    for w in dict_ls_days.keys():
        dict_df_rtpd_max[w] = pd.DataFrame(index = range(96), columns=['OPR_HR', 'OPR_INTERVAL'])
        dict_df_rtpd_min[w] = pd.DataFrame(index = range(96), columns=['OPR_HR', 'OPR_INTERVAL'])
        for d in dict_ls_days[w]:
            ls_tmp = list()
            for i in range(3):
                ls_tmp.append(
                    df.loc[
                        (df['OPR_DT']==d)&(df['OPR_INTERVAL']%3==i), 
                        'FORECAST_ERROR_Brtd_Artpd'
                    ].tolist()
                )
            dict_df_rtpd_max[w].loc[:, d] = np.array(ls_tmp).T.max(axis=1)
            dict_df_rtpd_min[w].loc[:, d] = np.array(ls_tmp).T.min(axis=1)
        dict_df_rtpd_max[w]['OPR_HR'] = df.loc[(df['OPR_DT']==d)&(df['OPR_INTERVAL']%3==0), 'OPR_HR'].tolist()
        dict_df_rtpd_min[w]['OPR_HR'] = dict_df_rtpd_max[w]['OPR_HR'].tolist()
        dict_df_rtpd_max[w]['OPR_INTERVAL'] = (
            df.loc[(df['OPR_DT']==d)&(df['OPR_INTERVAL']%3==0), 'OPR_INTERVAL']/3
        ).astype('int').tolist()
        dict_df_rtpd_min[w]['OPR_INTERVAL'] = dict_df_rtpd_max[w]['OPR_INTERVAL'].tolist()

    # Calculate uncertainties

    # First, add d + 1 to the studied period
    # d_plus1 = datetime.datetime.strptime(str_dend, '%Y-%m-%d') + datetime.timedelta(days=1)
    # if d_plus1.weekday() <= 4:
    #     dict_ls_days['WEEKDAY'].append(d_plus1.strftime('%Y-%m-%d'))
    # else:
    #     dict_ls_days['WEEKEND'].append(d_plus1.strftime('%Y-%m-%d'))

    ############################################################################
    # First, calculate 5 and 95 percentiles of NL forecast errors, note: it's hourly
    ############################################################################
    dict_df_uc = dict() # Uncertainty
    for w in dict_ls_days:
        dict_uc = {
            'OPR_DT':       list(),
            'OPR_HR':       list(),
            'P95_RTD':      list(),
            'P05_RTD':      list(),
            'P95_MAX_RTPD': list(),
            'P05_MIN_RTPD': list(),
        }
        for i in range(30, len(dict_ls_days[w])):
            d = dict_ls_days[w][i]
            selected_columns = dict_ls_days[w][i-30:i] + ['OPR_HR', 'OPR_INTERVAL']
            dict_p_rtd    = find_percentiles(dict_df_rtd[w].loc[:, selected_columns], [0.025, 0.975])
            dict_p95_rtpd = find_percentiles(dict_df_rtpd_max[w].loc[:, selected_columns], [0.975])
            dict_p05_rtpd = find_percentiles(dict_df_rtpd_min[w].loc[:, selected_columns], [0.025])
            dict_uc['P95_RTD'] += dict_p_rtd[0.975]
            dict_uc['P05_RTD'] += dict_p_rtd[0.025]
            dict_uc['P95_MAX_RTPD'] += dict_p95_rtpd[0.975]
            dict_uc['P05_MIN_RTPD'] += dict_p05_rtpd[0.025]
            dict_uc['OPR_HR'] += range(1, len(dict_p_rtd[0.975]) + 1)
            dict_uc['OPR_DT'] += [d for j in dict_p_rtd[0.975]]
        dict_df_uc[w] = pd.DataFrame(dict_uc)
    df_uc = pd.concat([dict_df_uc[k] for k in dict_df_uc.iterkeys()], ignore_index=True)
    df_uc = df_uc.sort_values(by=['OPR_DT', 'OPR_HR'])
    df_uc = df_uc.reset_index(drop=True)

    ############################################################################
    # Now we calculate the FRU/FRD requirements
    ############################################################################
    # RTD FRU and FRD
    i_df_rtd = df[df['OPR_DT']>='2017-04-16'].index
    i_uc_rtd = df_uc[df_uc['OPR_DT']>='2017-04-16'].index
    dict_fr_rtd = {
        'OPR_DT':          df.loc[i_df_rtd[0:-1], 'OPR_DT'].tolist(),
        'OPR_HR':          df.loc[i_df_rtd[0:-1], 'OPR_HR'].tolist(),
        'OPR_INTERVAL':    df.loc[i_df_rtd[0:-1], 'OPR_INTERVAL'].tolist(),
        'NET_LOAD_CHANGE': df.loc[i_df_rtd[1:], 'NET_LOAD_B_RTD'].to_numpy() - df.loc[i_df_rtd[0:-1], 'NET_LOAD_A_RTD'].to_numpy(),
        'P95_RTD':         df_uc.loc[i_uc_rtd, 'P95_RTD'].repeat(12).to_numpy()[0:-1], # -1 indicates the second last element in numpy
        'P05_RTD':         df_uc.loc[i_uc_rtd, 'P05_RTD'].repeat(12).to_numpy()[0:-1],
    }

    # RTPD FRU and FRD
    i_df_rtpd = df[(df['OPR_INTERVAL']%3==0)&(df['OPR_DT']>='2017-04-16')].index
    i_uc_rtpd = df_uc[df_uc['OPR_DT']>='2017-04-16'].index
    dict_fr_rtpd = {
        'OPR_DT':          df.loc[i_df_rtpd[0:-1], 'OPR_DT'].tolist(),
        'OPR_HR':          df.loc[i_df_rtpd[0:-1], 'OPR_HR'].tolist(),
        'OPR_INTERVAL':    (df.loc[i_df_rtpd[0:-1], 'OPR_INTERVAL']/3).astype('int').tolist(),
        'NET_LOAD_CHANGE': df.loc[i_df_rtpd[1:], 'NET_LOAD_B_RTD'].to_numpy() - df.loc[i_df_rtpd[0:-1], 'NET_LOAD_A_RTD'].to_numpy(),
        'P95_RTPD':        df_uc.loc[i_uc_rtpd, 'P95_RTD'].repeat(4).to_numpy()[0:-1],
        'P05_RTPD':        df_uc.loc[i_uc_rtpd, 'P05_RTD'].repeat(4).to_numpy()[0:-1],
    }

    df_fr_rtd = pd.DataFrame(dict_fr_rtd)
    df_fr_rtd['FRU'] = df_fr_rtd['NET_LOAD_CHANGE']+df_fr_rtd['P95_RTD']
    df_fr_rtd['FRD'] = df_fr_rtd['NET_LOAD_CHANGE']+df_fr_rtd['P05_RTD']
    df_fr_rtd.loc[df_fr_rtd['FRU']<0, 'FRU']=0
    df_fr_rtd.loc[df_fr_rtd['FRD']>0, 'FRD']=0

    df_fr_rtpd = pd.DataFrame(dict_fr_rtpd)
    df_fr_rtpd['FRU'] = df_fr_rtpd['NET_LOAD_CHANGE']+df_fr_rtpd['P95_RTPD']
    df_fr_rtpd['FRD'] = df_fr_rtpd['NET_LOAD_CHANGE']+df_fr_rtpd['P05_RTPD']
    df_fr_rtpd.loc[df_fr_rtpd['FRU']<0, 'FRU']=0
    df_fr_rtpd.loc[df_fr_rtpd['FRD']>0, 'FRD']=0

    # Plot it out
    plt.figure()
    plt.fill_between(df_fr_rtd.index.tolist(), df_fr_rtd['FRD'].tolist(), df_fr_rtd['FRU'].tolist(), facecolor='b')
    plt.title('RTD')
    plt.xlabel('# of 5-min interval')
    plt.ylabel('MW')

    plt.figure()
    plt.fill_between(df_fr_rtpd.index.tolist(), df_fr_rtpd['FRD'].tolist(), df_fr_rtpd['FRU'].tolist(), facecolor='b')
    plt.title('RTPD')
    plt.xlabel('# of 15-min interval')
    plt.ylabel('MW')

    plt.show()

    IP()

def estimate_validation():
    # Compare my estimates of uncertainties with CAISO OASIS data
    df_oasis    = pd.read_csv('for_flexiramp_summary.csv') # for_flexiramp_summary.csv can be created by function process_raw_for_flexiramp()
    df_oasis = df_oasis.loc[df_oasis['OPR_INTERVAL']==1, :] # Hourly data is used
    df_oasis['TIME_STAMP'] = pd.to_datetime(df_oasis['OPR_DT']+'-'+(df_oasis['OPR_HR']-1).map('{:02g}'.format), format='%Y-%m-%d-%H')

    df_est = pd.read_csv('df_uc.csv') # df_uc.csv is from pandas dataframe df_uc in function baseline_flexiramp()
    df_est['TIME_STAMP'] = pd.to_datetime(df_est['OPR_DT']+'-'+(df_est['OPR_HR']-1).map('{:02g}'.format), format='%Y-%m-%d-%H')

    # Starting index because of lack of several weekends in df_est
    i_oasis_s = df_oasis.loc[df_oasis['OPR_DT']=='2017-04-16', :].index.min()
    i_oasis_e = df_oasis.loc[df_oasis['OPR_DT']=='2019-02-28', :].index.max()
    i_est_s = df_est.loc[df_est['OPR_DT']=='2017-04-16', :].index.min()
    i_est_e = df_est.loc[df_est['OPR_DT']=='2019-02-28', :].index.max()

    df_fru_rtpd = pd.DataFrame(
        {
            'Time':  df_oasis.loc[i_oasis_s: i_oasis_e, 'TIME_STAMP'].tolist(),
            'OASIS': df_oasis.loc[i_oasis_s: i_oasis_e, 'UP_RTPD'].tolist(),
            'Estimate': df_est.loc[i_est_s: i_est_e, 'P95_MAX_RTPD'].tolist()
        }
    )
    df_frd_rtpd = pd.DataFrame(
        {
            'Time':  df_oasis.loc[i_oasis_s: i_oasis_e, 'TIME_STAMP'].tolist(),
            'OASIS': df_oasis.loc[i_oasis_s: i_oasis_e, 'DOWN_RTPD'].tolist(),
            'Estimate': (-df_est.loc[i_est_s: i_est_e, 'P05_MIN_RTPD']).tolist()
        }
    )

    df_fru_rtd = pd.DataFrame(
        {
            'Time':  df_oasis.loc[i_oasis_s: i_oasis_e, 'TIME_STAMP'].tolist(),
            'OASIS': df_oasis.loc[i_oasis_s: i_oasis_e, 'UP_RTD'].tolist(),
            'Estimate': df_est.loc[i_est_s: i_est_e, 'P95_RTD'].tolist()
        }
    )
    df_frd_rtd = pd.DataFrame(
        {
            'Time':  df_oasis.loc[i_oasis_s: i_oasis_e, 'TIME_STAMP'].tolist(),
            'OASIS': df_oasis.loc[i_oasis_s: i_oasis_e, 'DOWN_RTD'].tolist(),
            'Estimate': (-df_est.loc[i_est_s: i_est_e, 'P05_RTD']).tolist()
        }
    )

    df_fru_rtpd.plot(title='RTPD UP uncertainties', x='Time')
    df_frd_rtpd.plot(title='RTPD DN uncertainties', x='Time')
    df_fru_rtd.plot(title='RTD UP uncertainties', x='Time')
    df_frd_rtd.plot(title='RTD DN uncertainties', x='Time')
    plt.show()

###########################################################
# Regulation requirement analysis starts here

def baseline_reg_for_day(YYYY, MM, DD, oasis='DA'):
    '''
    Make graphs comparing REG with CAISO RU/RD history for day YYYY/MM/DD.
    We hope to compare REG* with CAISO's historical RU/RD requirements, however, 
    to calculate REG* we need the actual deployed RU/RD for each ACE interval, 
    which is one minute. Now we don't have the one-minute actually deployed 
    regulating reserve. So now the graph displayed REG instead of REG*.
    '''
    # Load ACE error data
    df_ace = pd.read_csv('for_reg_ACE.csv')
    df_ace['T'] = pd.to_datetime(df_ace['T'])

    # Load regulation requirements from CAISO OASIS
    df_reg_req = pd.read_csv('for_reg_reg_requirement.csv')

    df_sample = pd.DataFrame(
        {
            'OPR_HR': [i for i in range(1, 25) for j in range(60)],
            'OPR_INTERVAL': range(1, 61)*24,
        }
    )
    NX = len(df_sample) # Number of data point in a day

    # Plot time sequence area graph and 5% and 95% confidence levels.
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    for d in df_ace.loc[df_ace['T'].dt.month==MM, 'T'].dt.day.unique():
        i_selected = (df_ace['T'].dt.month==MM)&(df_ace['T'].dt.day==d)
        df_tmp = df_ace[i_selected]
        k = df_tmp['OPR_DT'].unique()[0]
        df_sample[k] = -df_tmp['ACE'].to_numpy() # Negative because REG = -ACE
        ax.fill_between(
            range(0, NX),
            0,
            -df_tmp['ACE'], # Negative because we compare -ACE and REG on the graph
            alpha=0.1,
            facecolor='b',
        )

    # Find confidence intervals
    dict_percentile = find_percentiles(df_sample, [0.05, 0.95])
    y05 = dict_percentile[0.05]
    y95 = dict_percentile[0.95]
    y95 = np.array(y95).repeat(60)
    y05 = np.array(y05).repeat(60)
    ax.plot(range(0, len(y95)), y95, 'r', linewidth=0.5, label='Baseline')
    ax.plot(range(0, len(y05)), y05, 'r', linewidth=0.5)
    ax.set_ylim([-1000, 1000])
    ax.set_xticks([60*i for i in range(0, 24)])
    ax.set_xticklabels(range(0, 24))
    # plt.title(MM)

    # Now, add CAISO requirements for that day
    i_select = df_reg_req['OPR_DT'] == '{:4g}-{:02g}-{:02g}'.format(YYYY, MM, DD)
    if oasis=='DA':
        ax.fill_between(
            range(0, NX), 
            df_reg_req.loc[i_select, 'RU_MIN_DAM'].repeat(15), 
            df_reg_req.loc[i_select, 'RU_MAX_DAM'].repeat(15), 
            alpha=0.2, 
            facecolor='k',
            label='OASIS',
        )
        ax.fill_between(
            range(0, NX), 
            -df_reg_req.loc[i_select, 'RD_MIN_DAM'].repeat(15), 
            -df_reg_req.loc[i_select, 'RD_MAX_DAM'].repeat(15), 
            alpha=0.2, 
            facecolor='k',
        )
    elif oasis=='RT':
        ax.fill_between(
            range(0, NX), 
            df_reg_req.loc[i_select, 'RU_MIN_RTM'].repeat(15), 
            df_reg_req.loc[i_select, 'RU_MAX_RTM'].repeat(15), 
            alpha=0.2, 
            color='k',
            label='OASIS',
        )
        ax.fill_between(
            range(0, NX), 
            -df_reg_req.loc[i_select, 'RD_MIN_RTM'].repeat(15), 
            -df_reg_req.loc[i_select, 'RD_MAX_RTM'].repeat(15), 
            alpha=0.2, 
            facecolor='k',
        )

    ax.set_title(datetime.datetime(YYYY, MM, DD).strftime('%m/%d/%y, %a'))
    ax.set_ylabel('MW')
    ax.set_xlabel('Hour')
    ax.legend()
    plt.show()

def lab_reg():
    from pandas.plotting import register_matplotlib_converters
    register_matplotlib_converters()

    # Load ACE error data
    df_ace = pd.read_csv('for_reg_ACE.csv')
    df_ace['T'] = pd.to_datetime(df_ace['T'])

    # Load regulation requirements from CAISO OASIS
    df_reg_req = pd.read_csv('for_reg_reg_requirement.csv')

    # Load regulation results from CAISO OASIS
    df_reg_res = pd.read_csv('for_reg_reg_results.csv')
    df_reg_res.loc[:, 'RU_TOT_AS_CAISO_EXP'] = df_reg_res[['RU_RTM_AS_CAISO_EXP', 'RU_DAM_AS_CAISO_EXP']].sum(axis=1)
    df_reg_res.loc[:, 'RD_TOT_AS_CAISO_EXP'] = df_reg_res[['RD_RTM_AS_CAISO_EXP', 'RD_DAM_AS_CAISO_EXP']].sum(axis=1)

    # Load net load probabilistic forecast
    df_prob = pd.read_csv(r'C:\Users\bxl180002\Downloads\RampSolar\Mucun\tobinghui\binghui_10min.csv')
    df_prob.loc[:, 'timestamp'] = pd.to_datetime(df_prob[['Year', 'Month', 'Day', 'Hour', 'Minute']])
    df_prob.set_index('timestamp', inplace=True)

    # Calculate 95 and 5 prediction interval
    df_prob.loc[:, 'Q97-D'] = df_prob['Q97'] - df_prob['Forecast']
    df_prob.loc[:, 'D-Q3']  = df_prob['Forecast'] - df_prob['Q3']

    df = pd.concat(
        [
            df_reg_req,
            df_reg_res[['RU_TOT_AS_CAISO_EXP', 'RD_TOT_AS_CAISO_EXP']],
        ],
        axis=1
    )

    df.loc[:, 'timestamp'] = pd.to_datetime(
        df['OPR_DT'].astype(str) 
        + 
        ' ' 
        + 
        (df['OPR_HR']-1).map('{:02g}'.format) 
        + 
        ':'
        + 
        ((df['OPR_INTERVAL']-1)*15).map('{:02g}'.format), 
        format='%Y-%m-%d %H:%M'
    )
    df.set_index('timestamp', inplace=True)

    year_selected  = 2019
    month_selected = 5

    for d in df.index[(df.index.year==year_selected) & (df.index.month==month_selected)].day.unique():
        fig = plt.figure()
        for year_selected in [2018, 2019]:
            ax = plt.subplot(1, 2, year_selected-2017)

            irows = ((df.index.year == year_selected) & (df.index.month == month_selected) & (df.index.day == d))

            ax.fill_between(
                df.index[irows],
                df_reg_req.loc[irows, 'RU_MIN_RTM'],
                df_reg_req.loc[irows, 'RU_MAX_RTM'],
                alpha=0.2,
                color='k',
            )
            ax.fill_between(
                df.index[irows],
                -df_reg_req.loc[irows, 'RD_MIN_RTM'],
                -df_reg_req.loc[irows, 'RD_MAX_RTM'],
                alpha=0.2,
                facecolor='k',
            )
            ax.plot(df.index[irows], df.loc[irows, 'RU_TOT_AS_CAISO_EXP'], 'r', linewidth=1)
            ax.plot(df.index[irows], -df.loc[irows, 'RD_TOT_AS_CAISO_EXP'], 'r', linewidth=1)

            if year_selected == 2018:
                ax.plot(
                    df_ace.loc[(df_ace['T'].dt.month==5) & (df_ace['T'].dt.day==d), 'T'],
                    df_ace.loc[(df_ace['T'].dt.month==5) & (df_ace['T'].dt.day==d), 'ACE'],
                )

            ax.set_title(str(year_selected) + '/' + str(month_selected) + '/' + str(d))

        # fig.savefig('./figs/' +str(d)+'.png')

        # Mucun's t-10 probabilistic prediction interval
        # ax.plot(df_prob.index[df_prob.index.day==d], df_prob.loc[df_prob.index.day==d, 'Q97-D'], 'b', linewidth=1)
        # ax.plot(df_prob.index[df_prob.index.day==d], -df_prob.loc[df_prob.index.day==d, 'D-Q3'], 'b', linewidth=1)
    plt.show()

    # Compare prob requirements with posted requirements
    for d in df.index[(df.index.year==year_selected) & (df.index.month==month_selected)].day.unique():
        fig = plt.figure()
        for year_selected in [2019]:
            ax = plt.subplot(1, 1, 1)

            irows = ((df.index.year == year_selected) & (df.index.month == month_selected) & (df.index.day == d))

            # ax.fill_between(
            #     df.index[irows],
            #     df_reg_req.loc[irows, 'RU_MIN_RTM'],
            #     df_reg_req.loc[irows, 'RU_MAX_RTM'],
            #     alpha=0.2,
            #     color='k',
            # )
            # ax.fill_between(
            #     df.index[irows],
            #     -df_reg_req.loc[irows, 'RD_MIN_RTM'],
            #     -df_reg_req.loc[irows, 'RD_MAX_RTM'],
            #     alpha=0.2,
            #     facecolor='k',
            # )
            h1 = df.loc[irows, 'RU_TOT_AS_CAISO_EXP'].plot(ax=ax, color='r')
            (-df.loc[irows, 'RD_TOT_AS_CAISO_EXP']).plot(ax=ax, color='r')

            if year_selected == 2018:
                ax.plot(
                    df_ace.loc[(df_ace['T'].dt.month==5) & (df_ace['T'].dt.day==d), 'T'],
                    df_ace.loc[(df_ace['T'].dt.month==5) & (df_ace['T'].dt.day==d), 'ACE'],
                )

            ax.set_title(str(year_selected) + '/' + str(month_selected) + '/' + str(d))


        # Mucun's t-10 probabilistic prediction interval
        df_prob.loc[df_prob.index.day==d, 'Q97-D'].plot(ax=ax, color='b')
        (-df_prob.loc[df_prob.index.day==d, 'D-Q3']).plot(ax=ax, color='b')
        # ax.plot(df_prob.index[df_prob.index.day==d], df_prob.loc[df_prob.index.day==d, 'Q97-D'], 'b', linewidth=1)
        # ax.plot(df_prob.index[df_prob.index.day==d], -df_prob.loc[df_prob.index.day==d, 'D-Q3'], 'b', linewidth=1)
        h2 = ax.set_xlabel('Time')
        ax.set_ylabel('MW')
        fig.savefig('./figs/reg_may' +str(d)+'.png')

        # plt.legend([h1, h2], ['OASIS', 'Prob'])

    plt.show()

    plt.figure()
    ax = plt.subplot(1, 1, 1)

    # ax.fill_between(
    #     df_prob.index,
    #     df_prob.loc[:, 'Q5'],
    #     df_prob.loc[:, 'Q95'],
    #     alpha=0.2,
    #     facecolor='k',
    # )
    # ax.fill_between(
    #     df_prob.index,
    #     df_prob.loc[:, 'Q25'],
    #     df_prob.loc[:, 'Q75'],
    #     alpha=0.4,
    #     facecolor='k',
    # )
    # ax.fill_between(
    #     df_prob.index,
    #     df_prob.loc[:, 'Q55'],
    #     df_prob.loc[:, 'Q65'],
    #     alpha=0.6,
    #     facecolor='k',
    # )

    # ax.plot(df_prob.index, df_prob['Deterministic'], 'b', linewidth=1)
    # ax.plot(df_prob.index, df_prob['Actual'], 'k', linewidth=1)

    # (df_prob['D-Q3']).plot(ax=ax) # Down uncertainty of net load, so we need reg-up
    # (-df_prob['Q97-D']).plot(ax=ax) # Up uncertainty of net load, so we need reg-dn
    # df.loc[df.index.month>=5, 'RU_TOT_AS_CAISO_EXP'].plot(ax=ax)
    # (-df.loc[df.index.month>=5, 'RD_TOT_AS_CAISO_EXP']).plot(ax=ax)
    # plt.show()

    # # Synthetic required reg up/down terms (REG*)
    # array_reg_up_syn = df.loc[df.index.year==2018, 'RU_TOT_AS_CAISO_EXP'].values.repeat(15, axis=0) - df_ace['ACE'].values
    # array_reg_dn_syn = -df.loc[df.index.year==2018, 'RD_TOT_AS_CAISO_EXP'].values.repeat(15, axis=0) - df_ace['ACE'].values
    # df_syn_reg_star = pd.DataFrame(
    #     np.column_stack([array_reg_up_syn, array_reg_dn_syn]),
    #     index=df_ace['T'],
    #     columns=['RU_STAR_SYN', 'RD_STAR_SYN'],
    # )

    # MM = 5
    # for selected_day in df.index[(df.index.year==2018) & (df.index.month==MM)].day.unique():

    #     df_sample_u = pd.DataFrame(
    #         {
    #             'OPR_HR': [i for i in range(1, 25) for j in range(60)],
    #             'OPR_INTERVAL': range(1, 61)*24,
    #         }
    #     )
    #     df_sample_d = pd.DataFrame(
    #         {
    #             'OPR_HR': [i for i in range(1, 25) for j in range(60)],
    #             'OPR_INTERVAL': range(1, 61)*24,
    #         }
    #     )
    #     plt.figure()
    #     ax = plt.subplot(1, 1, 1)
    #     for d in df_syn_reg_star.index[df_syn_reg_star.index.month==MM].day.unique():
    #         i_selected = (df_syn_reg_star.index.month==MM)&(df_syn_reg_star.index.day==d)
    #         df_tmp = df_syn_reg_star[i_selected]
    #         k = df_tmp.index.day.unique()[0]
    #         df_sample_u[k] = df_tmp['RU_STAR_SYN'].values
    #         df_sample_d[k] = df_tmp['RD_STAR_SYN'].values

    #     # Find confidence intervals
    #     dict_percentile_u = find_percentiles(df_sample_u, [0.05, 0.95])
    #     dict_percentile_d = find_percentiles(df_sample_d, [0.05, 0.95])
    #     reg_up = dict_percentile_u[0.95]
    #     reg_dn = dict_percentile_d[0.95]
    #     reg_up = np.array(reg_up).repeat(60)
    #     reg_dn = np.array(reg_dn).repeat(60)

    #     # df_reg_est = pd.DataFrame(
    #     #     np.column_stack([reg_up, reg_dn]),
    #     #     columns=['RU_EST', 'RD_EST'],
    #     #     index=pd.date_range(start='5/1/2018', end='5/2/2018', freq='min')[0:-1],
    #     # )

    #     ax.plot(range(0, len(reg_up)), reg_up, 'k', linewidth=0.5, label='Baseline')
    #     ax.plot(range(0, len(reg_dn)), reg_dn, 'k', linewidth=0.5)
    #     ax.set_ylim([-1000, 1000])
    #     ax.set_xticks([60*i for i in range(0, 24)])
    #     ax.set_xticklabels(range(0, 24))

    #     i_select = df_reg_req['OPR_DT'] == '{:4g}-{:02g}-{:02g}'.format(2019, MM, selected_day)
    #     ax.fill_between(
    #         range(0, i_selected.sum()), 
    #         df_reg_req.loc[i_select, 'RU_MIN_RTM'].repeat(15), 
    #         df_reg_req.loc[i_select, 'RU_MAX_RTM'].repeat(15), 
    #         alpha=0.2, 
    #         color='k',
    #         label='OASIS',
    #     )
    #     ax.fill_between(
    #         range(0, i_selected.sum()), 
    #         -df_reg_req.loc[i_select, 'RD_MIN_RTM'].repeat(15), 
    #         -df_reg_req.loc[i_select, 'RD_MAX_RTM'].repeat(15), 
    #         alpha=0.2, 
    #         facecolor='k',
    #     )
    
    # plt.show()

    IP()


if __name__ == '__main__':

    # dist_tx2kb_spdis()
    # read_nsrdb()
    # process_nsrdb()

    # Wind ramp AGC data generator
    ############################################################################
    # generate_agc_for_118()

    # IBM PARIS data processing
    ############################################################################
    # process_paris('./IBM/April', write_flag=False)
    # process_paris_more_quantiles('C:\\Users\\bxl180002\\git\\SF2\\IBM\\May.more_quantiles.5min', write_flag=False)
    # read_paris_April()
    # read_paris_May()
    # process_paris_global_solar_irradiance('./IBM/April', write_flag=False)
    # read_paris_5min()

    # CAISO OASIS data collection and save
    ############################################################################
    # collect_SLD_FCST_rtd()
    # collect_SLD_FCST_rtpd()
    # collect_AS_REQ_or_RESULTS()
    # collect_AS_REQ_or_RESULTS('res')
    # collect_SLD_REN_FCST_rtd()
    # collect_SLD_REN_FCST_rtpd()
    # collect_SLD_ADV_FCST()
    # collect_ENE_FLEX_RAMP_REQT()

    # CAISO flexiramp reserve analysis
    ############################################################################
    df = process_raw_for_flexiramp()
    # process_actual()
    # tmp_plot_forecast()
    # baseline_flexiramp_for_day(2019, 5, 31, use_persistence=True)
    # baseline_flexiramp()
    # estimate_validation()

    # CAISO regulation analysis
    ############################################################################
    # process_raw_reg_req()
    # process_raw_reg_results()
    # process_raw_ace_error()
    # baseline_reg_for_day(2019, 2, 22, oasis='DA')
    # lab_reg()

    # Others
    ############################################################################
    # prepare_data_for_mucun(write_flag=False)

    # IP()

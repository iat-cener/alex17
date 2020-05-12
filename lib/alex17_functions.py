import pandas as pd
import numpy as np
import datetime
import netCDF4
import utm
import xarray as xr
import rasterio as rio
import matplotlib.pyplot as plt 

from lib.variables_dictionary.variables import Variables
from lib.variables_dictionary.variables import nc_global_attributes_from_yaml

def read_sim(filename):
    M = xr.open_dataset(filename)
    U = M.eastward_wind
    V = M.northward_wind
    W = M.upward_air_velocity
    Th = M.air_potential_temperature
    TKE = M.specific_turbulent_kinetic_energy
    S = (U**2 + V**2)**0.5
    WD = (270-np.rad2deg(np.arctan2(V,U)))%360
    WD.attrs['units'] = 'degrees from North'
    WD.attrs['long_name'] = 'wind direction'
    S.attrs['units'] = 'm s-1'
    S.attrs['long_name'] = 'horizontal wind speed'
    return U, V, W, Th, TKE, S, WD

def read_obs(filename):
    M = xr.open_dataset(filename)
    S = M.wind_speed
    WD = M.wind_from_direction
    Sstd = M.wind_speed_std
    return S, Sstd, WD


def basemap_plot(src, masts, Z_transect, ref, ax):
    # Add overviews to raster to plot faster at lower resolution (https://rasterio.readthedocs.io/en/latest/topics/overviews.html)
    #from rasterio.enums import Resampling
    #factors = [2, 4, 8, 16]
    #dst = rio.open('./inputs/DTM_Alaiz_2m.tif', 'r+')
    #dst.build_overviews(factors, Resampling.average)
    #dst.update_tags(ns='rio_overview', resampling='average')
    #dst.close()
    A_ind = Z_transect['Name'].str.contains('A') # Tajonar ridge scan
    B_ind = Z_transect['Name'].str.contains('B') # Elortz valley scan
    C_ind = Z_transect['Name'].str.contains('C') # Alaiz ridge scan
    oview = src.overviews(1)[2] # choose overview (0 is largest, -1 is the smallest)
    topo = src.read(1, out_shape=(1, int(src.height // oview), int(src.width // oview)))
    spatial_extent = [src.bounds.left - ref[0], src.bounds.right - ref[0], src.bounds.bottom - ref[1], src.bounds.top - ref[1]]
    topo_ma = np.ma.masked_where(topo == 0 , topo, copy=True) 
    topo_ma.shape
    h_topo = ax.imshow(topo_ma, cmap='terrain', extent=spatial_extent, vmin=300, vmax=1200) #[left, right, bottom, top]
    h_masts = ax.scatter(masts['x'], masts['y'], s = 10, marker='s', c='k', label = 'Masts')
    h_A = ax.scatter(Z_transect[A_ind]['x'], Z_transect[A_ind]['y'], s = 2, marker='.', c='blue', label = 'A-transect')
    h_B = ax.scatter(Z_transect[B_ind]['x'], Z_transect[B_ind]['y'], s = 2, marker='.', c='black', label = 'B-transect')
    h_C = ax.scatter(Z_transect[C_ind]['x'], Z_transect[C_ind]['y'], s = 2, marker='.', c='red', label = 'C-transect')
    for i, txt in enumerate(masts['Name']):
        ax.annotate(txt, (masts['x'][i]+50, masts['y'][i]+50))   
    ax.set_title('ALEX17 sites')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.legend(handles = [h_masts,h_A,h_B,h_C])
    plt.colorbar(h_topo, ax = ax)
    return [h_masts,h_A,h_B,h_C]

    
def Ztransect_plot(masts, Z_transect, ax):
    A_ind = Z_transect['Name'].str.contains('A') # Tajonar ridge scan
    B_ind = Z_transect['Name'].str.contains('B') # Elortz valley scan
    C_ind = Z_transect['Name'].str.contains('C') # Alaiz ridge scan
    h_topoZ = (Z_transect['z']-125.).plot.area(stacked=False, color = 'lightgrey', alpha = 1, ax = ax)
    h_topoA = Z_transect[A_ind]['z'].plot.line(style = '.', color = 'blue', ms = 3, ax = ax)
    h_topoB = Z_transect[B_ind]['z'].plot.line(style = '.', color = 'black', ms = 3, ax = ax)
    h_topoC = Z_transect[C_ind]['z'].plot.line(style = '.', color = 'red', ms = 3, ax = ax)
    masts_inZ = [] # index of Z_transect position nearest to each mast
    for i, row in masts.iterrows():
        d = np.sqrt((Z_transect['x'] - masts['x'][i])**2 + (Z_transect['y'] - masts['y'][i])**2)
        masts_inZ.append(d[d == d.min()].index[0])
    for x in masts_inZ:
        ax.axvline(x, color = 'silver', linestyle = '--', zorder = 0)
    ax.set_title('Z-transect profile at 125 m above ground level')
    ax.set_xlabel('Z-transect position')
    ax.set_ylabel('z [m]')
    ax.set_axisbelow(True)
    ax.yaxis.grid(zorder = 0)
    ax.set_ylim([0,1000])
    ax.set_xticks(masts_inZ)
    ax.set_xticklabels(masts['Name'])
    
    return [h_topoZ, h_topoA, h_topoB, h_topoC]
    
def get_height_index(f_get_column):
    # height<1000m index
    max_height = 1000
    column = f_get_column(42.695, -1.558)
    all_heights = column.height.values
    last_height_id = np.searchsorted(all_heights, max_height) + 1
    height_index = all_heights[:last_height_id]
    return last_height_id, height_index


def create_nc_file(output_dataset, variables_dimensions, variables_indexes, var_dict, variables_to_write):
    for i in range(len(variables_dimensions)):
        _ = var_dict.nc_create_dimension(output_dataset, variables_dimensions[i], variables_indexes[i])
    variables = {}
    for var_name in variables_to_write:
        variables[var_name[0]] = var_dict.nc_create(output_dataset, var_name[0], variables_dimensions)
    return variables

def time_from_str(time_str):
    ###  str: 2018-09-30 00:00,
    ###       012345678901234567890
    year = int(time_str[0:4])
    month = int(time_str[5:7])
    day = int(time_str[8:10])
    hours = int(time_str[11:13])
    minutes = int(time_str[14:16])
    seconds = 0
    m_sec = 0
    # print(year, month, day, hours, minutes,seconds, m_sec)
    dt = datetime.datetime(year, month, day, hours, minutes, seconds, m_sec)

    seconds_since = (dt - datetime.datetime(1970, 1, 1)).total_seconds()
    return seconds_since


def create_alex17_file_1(file_1_path, f_get_column, variables_to_write):
    # initialise the variables dictionary class
    var_dict = Variables()

    # ---------------------------- create indexes
    # load the met mast locations
    met_mast_locations = pd.read_csv('./inputs/validation_profiles_XYZ.csv')
    utm_zone_number = 30
    utm_zone_letter = 'T'

    # height
    last_height_id, height_index = get_height_index(f_get_column)

    # time
    date_from = time_from_str("2018-09-30 00:00")
    date_to = time_from_str("2018-10-04 00:00")
    #date_to = time_from_str("2018-10-01 00:00")
    time_delta = 10 * 60
    #time_delta = 30 * 60
    time_index = list(range(int(date_from), int(date_to), int(time_delta)))

    print('1. simID_zt.nc:')
    # define the output file structure
    output_dataset = netCDF4.Dataset(file_1_path, 'w', format="NETCDF4")

    variables_indexes = [met_mast_locations['Name'].values, height_index, time_index]
    variables_dimensions = ('id', 'height', 'time')
    variables = create_nc_file(output_dataset, variables_dimensions, variables_indexes, var_dict, variables_to_write)
    nc_global_attributes_from_yaml(output_dataset, './config/Marinet.yaml')

    # fill the output file with results
    for index, row in met_mast_locations.iterrows():
        (lat, lon) = utm.to_latlon(float(row['easting[m]']), float(row['northing[m]']), utm_zone_number,
                                   utm_zone_letter)
        column = f_get_column(lat, lon)  # (time, height, variables)
        column = column.interp(time=time_index)  # interpolate into the desired time resolution
        for v_i in range(len(variables_to_write)):
            # (variables_name v_i)(id index, height :, time:) = (time :, height :, variables v_i).T
            variables[variables_to_write[v_i][0]][index, :, :] = column[:, :last_height_id, v_i].T
    output_dataset.close()


def create_alex17_file_2(file_path, f_get_point, variables_to_write):
    print('2. alex17_00a_box.nc:')

    utm_zone_number = 30
    utm_zone_letter = 'T'

    # initialise the variables dictionary class
    var_dict = Variables()

    # ---------------------------- create indexes
    # x - easting[m]
    easting_index = list(range(612000, 622000, 100))
    # y - northing[m]
    northing_index = list(range(4726000, 4736000, 100))
    # height
    height_index = [40, 125]
    #height_index = [40,]
    # time
    date_from = time_from_str("2018-09-30 00:00")
    date_to = time_from_str("2018-10-04 00:00")
    #date_to = time_from_str("2018-10-01 00:00")
    time_delta = 10 * 60
    #time_delta = 30 * 60
    time_index = list(range(int(date_from), int(date_to), int(time_delta)))

    # define the output file structure
    output_dataset = netCDF4.Dataset(file_path, 'w', format="NETCDF4")
    variables_indexes = [time_index, height_index, easting_index, northing_index]
    variables_dimensions = ('time', 'height', 'easting', 'northing')
    variables = create_nc_file(output_dataset, variables_dimensions, variables_indexes, var_dict, variables_to_write)
    nc_global_attributes_from_yaml(output_dataset, './config/Marinet.yaml')

    # fill the output file with results
    for i_x in range(len(easting_index)):
        x = easting_index[i_x]
        for i_y in range(len(northing_index)):
            y = northing_index[i_y]
            (lat, lon) = utm.to_latlon(x, y, utm_zone_number, utm_zone_letter)
            for i_h in range(len(height_index)):
                height = height_index[i_h]
                point = f_get_point(lat, lon, height)  # (time, variables)
                point = point.interp(time=time_index)  # interpolate into the desired time resolution
                for v_i in range(len(variables_to_write)):
                    # (variables_name v_i)('time' :, 'height' i_h, 'easting' i_x, 'northing' i_y) = (time :, variables v_i)
                    variables[variables_to_write[v_i][0]][:, i_h, i_x, i_y] = point[:, v_i]
    output_dataset.close()


def create_alex17_file_3(file_path, f_get_column, variables_to_write):
    print('3. simID_yzt.nc:')
    utm_zone_number = 30
    utm_zone_letter = 'T'

    # initialise the variables dictionary class
    var_dict = Variables()

    # ---------------------------- create indexes
    # y - northing[m]
    sampling_locations = pd.read_csv('./inputs/validation_ZTransect_XYZ.csv')

    # time
    date_from = time_from_str("2018-09-30 00:00")
    date_to = time_from_str("2018-10-04 00:00")
    #date_to = time_from_str("2018-10-01 00:00")
    time_delta = 10 * 60
    #time_delta = 30 * 60
    time_index = list(range(int(date_from), int(date_to), int(time_delta)))

    # height
    last_height_id, height_index = get_height_index(f_get_column)

    # define the output file structure
    output_dataset = netCDF4.Dataset(file_path, 'w', format="NETCDF4")
    variables_indexes = [sampling_locations['Name'].values, height_index, time_index]
    variables_dimensions = ('id', 'height', 'time')
    variables = create_nc_file(output_dataset, variables_dimensions, variables_indexes, var_dict, variables_to_write)
    nc_global_attributes_from_yaml(output_dataset, './config/Marinet.yaml')

    # fill the output file with results
    for index, row in sampling_locations.iterrows():
        (lat, lon) = utm.to_latlon(float(row['easting[m]']), float(row['northing[m]']), utm_zone_number,
                                   utm_zone_letter)
        column = f_get_column(lat, lon)  # (time, height, variables)
        column = column.interp(time=time_index)  # interpolate into the desired time resolution
        for v_i in range(len(variables_to_write)):
            # (variables_name v_i)(id index, height :, time:) = (time :, height :, variables v_i).T
            variables[variables_to_write[v_i][0]][index, :, :] = column[:, :last_height_id, v_i].T
    output_dataset.close()
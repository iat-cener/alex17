from netCDF4 import Dataset
import xarray as xr
from lib.variables_dictionary.variables import Variables
from lib.variables_dictionary.variables import nc_global_attributes_from_yaml
import pandas as pd
import numpy as np
import datetime
import netCDF4
import utm


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


def create_nc_file(output_dataset, variables_dimensions, variables_indexes, var_dict, variables_to_write):
    for i in range(len(variables_dimensions)):
        _ = var_dict.nc_create_dimension(output_dataset, variables_dimensions[i], variables_indexes[i])
    variables = {}
    for var_name in variables_to_write:
        variables[var_name] = var_dict.nc_create(output_dataset, var_name, variables_dimensions)
    return variables


# initialise the variables dictionary class
var_dict = Variables()

# ---------------------------- create indexes

# time
date_from = time_from_str("2018-09-30 00:00")
date_to = time_from_str("2018-10-04 00:00")
# date_to = time_from_str("2018-10-01 00:00")
time_delta = 10 * 60
# time_delta = 30 * 60
time_index = list(range(int(date_from), int(date_to), int(time_delta)))

cup_data = xr.open_dataset('./inputs/Nwinds_cup_masts.nc')
mp5_data = xr.open_dataset('./inputs/Nwinds_MP5_reference.nc')
sonic_data = xr.open_dataset('./inputs/Nwinds_sonic_masts.nc')

# height
height_index = np.array([10., 20., 40., 60., 80., 100., 120.])
met_masts = np.array(['MP5', 'M1', 'M2', 'M3', 'M5', 'M6', 'M7'])

variables_to_read = ['wind_speed', 'wind_direction', 'wind_speed_std']
variables_to_write = ['wind_speed', 'wind_from_direction', 'wind_speed_std']


output_dataset = netCDF4.Dataset('observations/masts_obs.nc', 'w', format="NETCDF4")

variables_indexes = [met_masts, height_index, time_index]
variables_dimensions = ('id', 'height', 'time')
variables = create_nc_file(output_dataset, variables_dimensions, variables_indexes, var_dict, variables_to_write)

nc_global_attributes_from_yaml(output_dataset, './observations/observations.yaml')

# fill the output file with results
for i_met_mast in range(len(met_masts)):
    met_mast = met_masts[i_met_mast]
    if met_mast == 'MP5':
        column = mp5_data  # (time, height, variables)
        heights = mp5_data['Height'].values
    elif met_mast == 'M1':
        column = cup_data.sel(Mast_ID=0)  # (time, height, variables)
        heights = cup_data['Height'].values
    elif met_mast == 'M5':
        column = cup_data.sel(Mast_ID=1)  # (time, height, variables)
        heights = cup_data['Height'].values
    elif met_mast == 'M2':
        column = sonic_data.sel(Mast_ID=0)  # (time, height, variables)
        heights = sonic_data['Height'].values
    elif met_mast == 'M3':
        column = sonic_data.sel(Mast_ID=1)  # (time, height, variables)
        heights = sonic_data['Height'].values
    elif met_mast == 'M6':
        column = sonic_data.sel(Mast_ID=2)  # (time, height, variables)
        heights = sonic_data['Height'].values
    elif met_mast == 'M7':
        column = sonic_data.sel(Mast_ID=3)  # (time, height, variables)
        heights = sonic_data['Height'].values

    # column = column.interp(Time=time_index)  # interpolate into the desired time resolution
    for i_var in range(len(variables_to_write)):
        var_name = variables_to_write[i_var]
        for height in heights:
            # (variables_name v_i)(id index, height :, time:) = (time :, height :, variables v_i).T
            variables[var_name][i_met_mast, np.where(height_index == height)[0][0], :] = \
                column[variables_to_read[i_var]].loc[height, :].values
output_dataset.close()

print()


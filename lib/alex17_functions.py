import pandas as pd
import numpy as np
import datetime
import netCDF4
import utm
import xarray as xr

from lib.variables_dictionary.variables import Variables
from lib.WrfReader import WrfReader


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
        _ = create_dimension(output_dataset, variables_dimensions[i], var_dict, variables_indexes[i])
    variables = {}
    for var_name in variables_to_write:
        variables[var_name[0]] = create_variable(output_dataset, var_name[0], variables_dimensions, var_dict)
    return variables


def create_dimension(output_dataset, var_name, var_dict, data):
    output_dataset.createDimension(var_name, len(data))
    dim = var_dict.nc_create(output_dataset, var_name, (var_name,))
    dim[:] = data
    return dim


def create_variable(output_dataset, variable_name, dimensions, var_dict):
    # create the variable
    var = var_dict.nc_create(
        output_dataset,
        variable_name,
        dimensions)
    return var


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
    met_mast_locations = pd.read_csv('./ALEX17_inputs/validation_profiles_XYZ.csv')
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
    print('2. simID_xyt.nc:')

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
    sampling_locations = pd.read_csv('../ALEX17_inputs/validation_ZTransect_XYZ.csv')
    sampling_locations.sort_values('northing[m]', inplace=True)
    northing_index = sampling_locations['northing[m]'].values
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
    variables_indexes = [time_index, height_index, northing_index]
    variables_dimensions = ('time', 'height', 'northing')
    variables = create_nc_file(output_dataset, variables_dimensions, variables_indexes, var_dict, variables_to_write)

    # fill the output file with results
    for i_y in range(len(northing_index)):
        y = northing_index[i_y]
        x = sampling_locations['easting[m]'][i_y]
        (lat, lon) = utm.to_latlon(x, y, utm_zone_number, utm_zone_letter)
        column = f_get_column(lat, lon)  # (time, height, variables)
        column = column.interp(time=time_index)  # interpolate into the desired time resolution
        for v_i in range(len(variables_to_write)):
            # (variables_name v_i)(time :, height :, northing i_y) = (time :, height :, variables v_i).T
            variables[variables_to_write[v_i][0]][:, :, i_y] = column[:, :last_height_id, v_i]
    output_dataset.close()
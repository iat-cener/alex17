#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 17:05:28 2020
Tool for handling WRF files
@authors
    : Pawel Gancarski (pgancarski@cener.com)
    : Javier Sanz Rodrigo (jsrodrigo@cener.com)
    : Roberto Chavez (roberto.chavez@ul.com)
"""
import utm, glob, netCDF4, datetime
import numpy as np
import pandas as pd
import xarray as xr
from functools import lru_cache
from numba import jit
from wrf import extract_times, extract_dim, ALL_TIMES, ll_to_xy
from lib.nc_read_functions import get_zagl, readAllvars, get_index_of_subset_domain

# %%

@jit(nopython=True)
def array_multiply_numba(a, x):
    return a * x


@jit(nopython=True)
def array_add_numba(a, b):
    return a + b


class WrfReader:
    '''
    Class to read the wrf output files and store a subset of the data varirables.
    The subset is used to make it RAM efficient as

	Parameters
	----------
	wrf_files_path : str
		path of wrf output files

	wrf_domain_number : str
		number of the wrf domain for the data to get extracted.

	dateLimits : list
		2-element list with the datefrom - dates to which the data is to be extracted

	variables_to_extract : list
		Name of the WRF variables to extract from the wrfout files.
        You should write the names as they are stored in the ncfile except for the following variable:
            - Th = potential temperature which is calculated as:
                T (perturbation potential temperature) + T0(base state temperature)
            - TENDENCIES=advection and pressure gradient terms from the momentum budget,
              and potential temp. budget which are assumed to be stored under the names:
                'RU_TEND_ADV','RV_TEND_ADV','RU_TEND_PGF','RV_TEND_PGF','T_TEND_ADV'

	subset_of_wrfDomain -optional: list
		3-element list with the [lat,lon, L=optinal] associated with the central coordinates and extension
        of a subset domain to make to whole WRF reading more memory efficient
    '''

    def __init__(self, wrf_files_path, wrf_domain_number, dateLimits, variables_to_extract, subset_of_wrfDomain=None):

        datefrom = np.datetime64(pd.to_datetime(dateLimits[0], format='%Y-%m-%d %H:%M:%S', errors='ignore'))
        dateto = np.datetime64(pd.to_datetime(dateLimits[1], format='%Y-%m-%d %H:%M:%S', errors='ignore'))

        # So far, spinup discarded hours are hardcoded
        spinup = 0

        '''
        List with all files to read. This part of code may difer depending on how you name and store
        the wrf output files. The nc_files list sould contain all the files over which the code
        is going to loop and concatenate variables in "Time" dimension axis.
        '''
        nc_files_list = []
        folders = [sorted(glob.glob(wrf_files_path + '*'))]
        # notice that it is assumed that wrf files are named with a convention of wrfout_d0xx
        fileN = 'wrfout_d0' + str(wrf_domain_number)
        for ii in folders[0]:
            # wrfoutfiles.append(sorted(glob.glob(ii+'/WRF/'+fileN+'*'))[0])
            if ii.split('/')[-1].startswith(fileN):
                nc_files_list.append(ii)

        Nfiles = len(nc_files_list)

        print(str(Nfiles) + ' files found for wrfout_d0' + str(wrf_domain_number) + ' outputs')
        print("NOTE. It is assumed that every ncfile contains an initialized Run for which " +
              str(spinup) + " h of spinup are discarded")

        print('The following variables are to be extracted ')
        print(variables_to_extract)

        # loop over all files found
        for ii, file in enumerate(nc_files_list):

            print("Reading netCDF data from file:", file)

            # open the netcdf file
            f2 = netCDF4.Dataset(file)

            # get the netCDF time data
            ncTimes = extract_times(f2, timeidx=ALL_TIMES, squeeze=True)

            spinupDate = ncTimes[0] + np.timedelta64(spinup, 'h')

            # protects from empty dates
            if sum(np.logical_and(ncTimes >= datefrom, ncTimes <= dateto)) != 0:
                if ii == 0:
                    # if no subset_of_wrfDomain is chosen, it extracts all domain grid points
                    if subset_of_wrfDomain is None:
                        iBT = np.arange(0, extract_dim(f2, 'bottom_top'))
                        iSN = np.arange(0, extract_dim(f2, 'south_north'))
                        iWE = np.arange(0, extract_dim(f2, 'west_east'))
                    else:
                        lat_s, lon_s, L = subset_of_wrfDomain[0], subset_of_wrfDomain[1], subset_of_wrfDomain[2]
                        # get the indexes of a subset domain for memory efficient purposes
                        iBT, iSN, iWE = get_index_of_subset_domain(f2, lat_s, lon_s, L)

                    # makes sure spinup dates are discarded while filtering only period between desired dates
                    iTimes = np.logical_and.reduce((ncTimes >= spinupDate, ncTimes >= datefrom, ncTimes <= dateto))

                    # set the valid times
                    self.times = ncTimes[iTimes]
                    self.min_i_lat = min(iSN)
                    self.min_i_lon = min(iWE)

                    self.lat_lon = np.concatenate((f2.variables.get('XLAT')[0:1, iSN, iWE],
                                                   f2.variables.get('XLONG')[0:1, iSN, iWE]), axis=0)

                    # Build z above ground
                    zagl = get_zagl(f2, iTimes, iBT, iSN, iWE)

                    self.variables_data = readAllvars(f2, variables_to_extract, iTimes, iBT, iSN, iWE)
                    self.variables_names = self.variables_data.keys()
                    self.input_file = f2

                else:

                    iTimes = np.logical_and.reduce((ncTimes > spinupDate, ncTimes >= datefrom, ncTimes <= dateto))

                    self.times = np.concatenate((self.times, ncTimes[iTimes]))

                    # concatenate z above ground
                    zaglTmp = get_zagl(f2, iTimes, iBT, iSN, iWE)
                    zagl = np.concatenate((zagl, zaglTmp), axis=0)

                    # concatenate every selected variable
                    varsTmp = readAllvars(f2, variables_to_extract, iTimes, iBT, iSN, iWE)

                    for iVar in varsTmp.keys():
                        self.variables_data[iVar] = np.concatenate((self.variables_data[iVar], varsTmp[iVar]), axis=0)

                    f2.close()
            else:
                print('File : ' + file + 'does not contain any valid dates')

        # average-out the height as the variability of the actual values is low
        self.heights = np.nanmean(zagl, axis=(0, 2, 3))

        # Reference Time for the netcdf files
        self.referenceDate = np.datetime64(datetime.datetime(1970, 1, 1, 0, 0, 0))

        # convert time stamp to "seconds since <referenceDate>"
        self.seconds = pd.Series(self.times - self.referenceDate).dt.total_seconds().values

        # number of grid points in the subset wrf box
        self.nt, self.nz, self.ny, self.nx = zagl.shape

    def __del__(self):
        if hasattr(self, 'input_file'):
            self.input_file.close()

    @jit(forceobj=True)
    def get_indexes(self, lat, lon, height):
        i_lon, i_lat = ll_to_xy(self.input_file, lat, lon)
        i_h = np.searchsorted(self.heights, height) - 1
        i_lat = i_lat - self.min_i_lat
        i_lon = i_lon - self.min_i_lon
        return i_lat, i_lon, i_h

    @jit(forceobj=True)
    def get_nearest_points(self, lat, lon, height):
        i_lat, i_lon, i_h = self.get_indexes(lat, lon, height)
        x, y, _, _ = utm.from_latlon(lat, lon)
        z = height

        x_near, y_near = self.utm_from_indexes((i_lat, i_lon))

        if x_near > x:  # East from sampling point
            i_lon = i_lon - 1

        if y_near > y:  # North from sampling point
            i_lat = i_lat - 1

        sw_point = (i_lat, i_lon, i_h)
        se_point = (i_lat, i_lon + 1, i_h)
        nw_point = (i_lat + 1, i_lon, i_h)
        ne_point = (i_lat + 1, i_lon + 1, i_h)

        top_sw_point = (i_lat, i_lon, i_h + 1)
        top_se_point = (i_lat, i_lon + 1, i_h + 1)
        top_nw_point = (i_lat + 1, i_lon, i_h + 1)
        top_ne_point = (i_lat + 1, i_lon + 1, i_h + 1)

        x_f, y_f = self.utm_from_indexes(ne_point)
        x_l, y_l = self.utm_from_indexes(sw_point)
        z_l, z_f = (self.heights[i_h], self.heights[i_h + 1])
        '''
        print("-------------")
        print(x, y)
        print(x_f, y_f)
        print(x_l, y_l)
        '''

        alpha_x = (x - x_l) / (x_f - x_l)
        alpha_y = (y - y_l) / (y_f - y_l)
        alpha_z = (z - z_l) / (z_f - z_l)

        points_indexes = np.array([
            sw_point, se_point,
            nw_point, ne_point,
            top_sw_point, top_se_point,
            top_nw_point, top_ne_point
        ])

        weights = np.array([
            (1 - alpha_x) * (1 - alpha_y) * (1 - alpha_z), alpha_x * (1 - alpha_y) * (1 - alpha_z),
            (1 - alpha_x) * alpha_y * (1 - alpha_z), alpha_x * alpha_y * (1 - alpha_z),
            (1 - alpha_x) * (1 - alpha_y) * alpha_z, alpha_x * (1 - alpha_y) * alpha_z,
            (1 - alpha_x) * alpha_y * alpha_z, alpha_x * alpha_y * alpha_z
        ])
        # print(weights.sum())
        return points_indexes, weights

    def lat_lon_from_indexes(self, p):
        [lat, lon] = self.lat_lon[:, p[0], p[1]]
        # print(utm.from_latlon(lat, lon))
        return lat, lon

    def utm_from_indexes(self, p):
        lat, lon = self.lat_lon_from_indexes(p)
        x, y, _, _ = utm.from_latlon(lat, lon)
        return x, y

    def utm_h_from_indexes(self, p):
        x, y = self.utm_from_indexes(p)
        z = self.heights[p[2]]
        return x, y, z

    @lru_cache(maxsize=None)
    def get_all_vals(self, i_lat, i_lon, i_h):
        out_dict = {}
        for var_name in self.variables_names:
            if len(self.variables_data[var_name].shape) == 4:
                out_dict[var_name] = self.variables_data[var_name][:, i_h, i_lat, i_lon]
            else:
                out_dict[var_name] = self.variables_data[var_name][:, i_lat, i_lon]
        return pd.DataFrame(out_dict, index=self.seconds)

    def interp_3d(self, lat, lon, height):
        x, y, _, _ = utm.from_latlon(lat, lon)
        nearest_points_indexes, weights = self.get_nearest_points(lat, lon, height)
        [i_lat, i_lon, i_h] = nearest_points_indexes[0]
        vals = self.get_all_vals(i_lat, i_lon, i_h)
        out = array_multiply_numba(vals.to_numpy(), weights[0])
        for i in range(1, len(nearest_points_indexes)):
            [i_lat, i_lon, i_h] = nearest_points_indexes[i]
            out = array_add_numba(
                out,
                array_multiply_numba(self.get_all_vals(i_lat, i_lon, i_h).to_numpy(),
                                     weights[i])
            )

        return pd.DataFrame(out, index=self.seconds)

    def get_time(self, ):
        return self.seconds

    def get_height(self, ):
        return self.heights

    def get_data(self, varName):
        return self.variables_data[varName]

    def get_point(self, lat, lon, height):

        out = self.interp_3d(lat, lon, height)

        coords = {'time': self.seconds,
                  'variables': [x for x in self.variables_names]}
        dims = ('time', 'variables')
        out = xr.DataArray(out, coords=coords, dims=dims)
        return out

    def get_column(self, lat, lon):
        shape = (len(self.seconds),
                 len(self.heights),
                 len(self.variables_names))
        coords = {'time': self.seconds,
                  'height': self.heights,
                  'variables': [x for x in self.variables_names]}
        dims = ('time', 'height', 'variables')
        out = xr.DataArray(np.empty(shape), coords=coords, dims=dims)

        for i_h in range(len(self.heights) - 1):
            out[:, i_h, :] = self.get_point(lat, lon, self.heights[i_h])
        return out

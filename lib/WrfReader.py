import os
import utm
import glob
import numpy as np
import pandas as pd
import netCDF4
import datetime
import xarray as xr
import scipy.spatial.distance as scipy_dist

from wrf import extract_times, extract_dim, ALL_TIMES, ll_to_xy, xy_to_ll
from lib.aux_functions import *


def build_lat_lon_table(nc_file, i_we, i_sn):
    out = np.empty([len(i_we), len(i_sn), 2])

    for i in range(len(i_we)):
        for j in range(len(i_sn)):
            out[i, j] = xy_to_ll(nc_file, i_we[i], i_sn[j])
    return out


class WrfReader:
    """
    Tool for handling WRF files


    @authors
        : Pawel Gancarski (pgancarski@cener.com)
        : Javier Sanz Rodrigo (jsrodrigo@cener.com)
        : Pedro Correia (pmfernandez@cener.com)
        : Roberto Chavez (rchavez.arroyo@gmail.com)
    """

    def __init__(self, variables_to_write):
        # %% --------------- INPUT INFORMATION AND CASE DEFINITIONS ------------------------------------
        '''
        INPUTS:

        dom = 2			     # domain to extract time series from wrfout_d0[dom]*
        #dxmeso = 9000.0	 # horizontal resolution of the mesoscale grid
        #L = 45000.0         # Length of microscale domain Lx = Ly for spatial avaraging
                             # L = 0            linear interpolation to site coordinates
                             # 0 < L < dxmeso   nearest grid point
                             # L > dxmeso       spatial average over Nav points
        '''
        # Station coordinates
        # identifier. only serves for naming files (or folders) to output the data
        siteID = 'Alaiz'
        # degrees N
        lat_s = 42.695
        # degrees E
        lon_s = -1.558
        # Length of the box to stract data (horizontally)
        L = 45000
        # domain number of the WRF
        dom = 3

        if dom == 1:
            dxmeso = 27000
        elif dom == 2:
            dxmeso = 9000.0
        elif dom == 3:
            dxmeso = 3000

        print(
            'processing ' + siteID + ' lon:' + str(lon_s) + ' lat:' + str(lat_s) + ' L:' + str(L) + ' dom:' + str(dom))

        # ----------------------------------------------------------------------
        # Load simulation data
        # ----------------------------------------------------------------------
        # Simulation filename
        wrfOutput = './ALEX17_inputs/'
        outputFolder = './out/' + siteID + '/'

        # Define the dates from and to in str format for which you want to limit
        # the extraction of wrf data
        datefromSys = "2018-09-30 00:00"
        datetoSys = "2018-10-04 00:00"

        # helps to define the name of the output
        simID = ''
        pbl = 'MYNN2.5'

        # Spinup discarded hours
        spinup = 24

        Nav = int(L / dxmeso) + 1  # number of points to include in the spatial averaging box
        simdes = ': WRF3.8.1, NEWA setup, ERA5, ' + pbl + ' 27km > 9km >3km. Data from d0' + str(dom) + ' at ' + str(
            int(dxmeso / 1000)) + ' km res, Lav = ' + str(int(L / 1000)) + ' km'

        '''
        Name of the WRF variables to extract from the wrfout files.
        You should write the names as they are stored in the ncfile except for the following variables:
        - POT = potential temperature which is calculated from T (perturbation potential temperature) + T0(base state temperature)
        - TENDENCIES=advection and pressure gradient terms from the momentum budget, and potential temp. budget which
                    are assumed to be stored under the names: 'RU_TEND_ADV','RV_TEND_ADV','RU_TEND_PGF','RV_TEND_PGF','T_TEND_ADV'
        '''
        varsWRF = ['T2', 'TSK', 'UST', 'PSFC', 'HFX', 'LH', 'RMOL', 'U', 'V', 'W', 'POT', 'TENDENCIES', 'TKE_PBL']
        self.variables_to_write = variables_to_write

        # Coriolis calculation
        fc = coriolis(lat_s)

        # latlon to utm
        utmX_s, utmY_s, utm_zonenumber_s, utm_zoneletter_s = utm.from_latlon(lat_s, lon_s)

        '''
        List with all files to read. Tthis part of code may difer depending on 
        how you name and store the wrf output files the main idea is that the wrfoutfiles 
        list, contains all the files over which the code is going to loop and concatenate
        in "Time" dimension the variables. e.g. wrfout is stored by weekly simulations
        and you need several days contained between the end of one and begining of the second, 
        then both files should be listed in the wrfoutfiles list.  Or maybe you want several weeks
        of data concatenated, then all the wrf with those data sould be listed here.
        '''
        wrfoutfiles = []
        folders = [sorted(glob.glob(wrfOutput + '*'))]
        fileN = 'wrfout_d0' + str(dom)
        for ii in folders[0]:
            # wrfoutfiles.append(sorted(glob.glob(ii+'/WRF/'+fileN+'*'))[0])
            if ii.split('/')[-1].startswith(fileN):
                wrfoutfiles.append(ii)

        Nfiles = len(wrfoutfiles)

        print(str(Nfiles) + ' ' + " files found for L=" + str(L) + ' and d0' + str(dom))
        print("NOTE. It is assumed that every ncfile contains an initialized Run for which " + str(
            spinup) + " h of spinup are discarded")

        # Create output directory
        if not os.path.exists(outputFolder):
            os.makedirs(outputFolder)

        datefrom = np.datetime64(pd.to_datetime(datefromSys, format='%Y-%m-%d %H:%M'))
        dateto = np.datetime64(pd.to_datetime(datetoSys, format='%Y-%m-%d %H:%M'))

        # loop over all files found
        for ii, file in enumerate(wrfoutfiles):

            # ii=0; file = wrfoutfiles[ii]
            print("Processing file :", file)

            f2 = netCDF4.Dataset(file)

            # get the netCDF time data
            ncTimes = extract_times(f2, timeidx=ALL_TIMES, squeeze=True)

            spinupDate = ncTimes[0] + np.timedelta64(spinup, 'h')

            # protects from empty dates
            if sum(np.logical_and(ncTimes >= datefrom, ncTimes <= dateto)) != 0:
                if ii == 0:
                    # makes sure spinup dates are discared while filtering only period between desired dates
                    self.iTimes = np.logical_and.reduce((ncTimes >= spinupDate, ncTimes >= datefrom, ncTimes <= dateto))

                    # nummber of grid centered points
                    zdim = extract_dim(f2, 'bottom_top')
                    xdim = extract_dim(f2, 'west_east')
                    ydim = extract_dim(f2, 'south_north')

                    inear_WE, inear_SN = ll_to_xy(f2, lat_s, lon_s)

                    iWE, iSN = getneighbours(Nav, inear_WE, inear_SN)

                    self.min_i_lat = min(iWE)
                    self.min_i_lon = min(iSN)
                    iBT = np.arange(0, zdim)
                    self.lat_lon = build_lat_lon_table(f2, iWE, iSN)
                    timesSim = ncTimes[self.iTimes]

                    # Build z above ground
                    zagl = get_zagl(f2, self.iTimes, iBT, iSN, iWE)

                    self.vDict = readAllvars(f2, varsWRF, self.iTimes, iBT, iSN, iWE, lat_s)
                    self.input_file = f2

                else:

                    self.iTimes = np.logical_and.reduce((ncTimes > spinupDate, ncTimes >= datefrom, ncTimes <= dateto))

                    timesSim = np.concatenate((timesSim, ncTimes[iTimes]))

                    # concatenate z above ground
                    zaglTmp = get_zagl(f2, iTimes, iBT, iSN, iWE)
                    zagl = np.concatenate((zagl, zaglTmp), axis=0)

                    # concatenate every selected variable
                    vDictTmp = readAllvars(f2, varsWRF, iTimes, iBT, iSN, iWE, lat_s)

                    for iVar in vDictTmp.keys():
                        self.vDict[iVar] = np.concatenate((vDict[iVar], vDictTmp[iVar]), axis=0)
                    f2.close()
            else:
                print('File : ' + file + 'does not contain any valid dates')

        # average-out the height as the variability of the actual values is low
        self.heights = np.nanmean(zagl, axis=(0, 2, 3))

        # Reference Time for the netcdf files
        dateRef = np.datetime64(datetime.datetime(1970, 1, 1, 0, 0, 0))

        # convert time stamp to "seconds since <dateRef>"
        self.secsDiff = pd.Series(timesSim - dateRef).dt.total_seconds().values

        # TODO a hack here
        # variables_to_write = [['eastward_wind', 'U'],
        #                       ['northward_wind', 'V'],
        #                       ['upward_air_velocity', 'W'],
        #                       ['air_potential_temperature', 'Th'],
        #                       ['TKE', 'TKE']]

        # rename one invented variable name into another
        self.vDict['Th'] = self.vDict.pop('POT')
        self.vDict['TKE'] = self.vDict.pop('TKE_PBL')

    def __del__(self):
        if hasattr(self, 'input_file'):
            self.input_file.close()

    def get_indexes(self, lat, lon, height):
        i_lat, i_lon = ll_to_xy(self.input_file, lat, lon)
        i_h = np.searchsorted(self.heights, height)
        i_lat = i_lat - self.min_i_lat
        i_lon = i_lon - self.min_i_lon
        return i_lat, i_lon, i_h

    def get_nearest_points(self, lat, lon, height):
        i_lat, i_lon, i_h = self.get_indexes(lat, lon, height)
        x, y, _, _ = utm.from_latlon(lat, lon)
        z = height

        x_near, y_near = self.utm_from_indexes((i_lat, i_lon))

        if x_near > x:  # East from sampling point
            i_lat = i_lat - 1

        if y_near > y:  # North from sampling point
            i_lon = i_lon - 1

        sw_point = (i_lat, i_lon, i_h)
        se_point = (i_lat + 1, i_lon, i_h)
        nw_point = (i_lat, i_lon + 1, i_h)
        ne_point = (i_lat + 1, i_lon + 1, i_h)

        top_sw_point = (i_lat, i_lon, i_h + 1)
        top_se_point = (i_lat + 1, i_lon, i_h + 1)
        top_nw_point = (i_lat, i_lon + 1, i_h + 1)
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
        [lat, lon] = self.lat_lon[p[0], p[1]]
        return lat, lon

    def utm_from_indexes(self, p):
        [lat, lon] = self.lat_lon[p[0], p[1]]
        x, y, _, _ = utm.from_latlon(lat, lon)
        return x, y

    def utm_h_from_indexes(self, p):
        x, y = self.utm_from_indexes(p)
        z = self.heights[p[2]]
        return x, y, z

    def get_all_vals(self, indexes):
        out_dict = {}
        for _, var_name in self.variables_to_write:
            if len(self.vDict[var_name].shape) == 4:
                out_dict[var_name] = self.vDict[var_name][:, indexes[2], indexes[0], indexes[1]]
            else:
                out_dict[var_name] = self.vDict[var_name][:, indexes[0], indexes[1]]
        return pd.DataFrame(out_dict, index=self.secsDiff)

    def interp_3d(self, lat, lon, height):
        x, y, _, _ = utm.from_latlon(lat, lon)

        nearest_points_indexes, weights = self.get_nearest_points(lat, lon, height)

        out = self.get_all_vals(nearest_points_indexes[0]) * weights[0]
        for i in range(1, len(nearest_points_indexes)):
            out = out + self.get_all_vals(nearest_points_indexes[i]) * weights[i]

        return out

    def get_time(self, ):
        return self.secsDiff

    def get_height(self, ):
        return self.heights

    def get_point(self, lat, lon, height):

        out = self.interp_3d(lat, lon, height)

        coords = {'time': self.secsDiff,
                  'variables': [x[1] for x in self.variables_to_write]}
        dims = ('time', 'variables')
        out = xr.DataArray(out, coords=coords, dims=dims)
        return out

    def get_column(self, lat, lon):
        shape = (len(self.secsDiff),
                 len(self.heights),
                 len(self.variables_to_write))
        coords = {'time': self.secsDiff,
                  'height': self.heights,
                  'variables': [x[1] for x in self.variables_to_write]}
        dims = ('time', 'height', 'variables')
        out = xr.DataArray(np.empty(shape), coords=coords, dims=dims)

        for i_h in range(len(self.heights)-1):
            out[:, i_h, :] = self.get_point(lat, lon, self.heights[i_h])
        return out

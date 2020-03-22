#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Set of auxliary functions for the WRF_read script
author:    Roberto Chavez <roberto.chavez@ul.com>  march/2020

'''

import numpy as np
import pandas as pd
import datetime
from scipy.spatial import cKDTree
import utm
from wrf import extract_dim, ll_to_xy, xy_to_ll


def getVarEff(ncfid, varName, iTimes, iBT, iSN, iWE):
    '''
	Memmory efficient method of getting the nc data and destagered (if
	that is the case ) into a centered grid.
	NOTE: it assumes that the variables stored in the netcdf have as
	first dimension of the Time dimension

	Parameters
	----------
	ncfid : file identifier class
		netcdf file identifier.
	varName : str
		name of the variable to extract from the netcdf. Should match exactly the name
		of the variable as stored in the netcdf.
	iTimes : int, logic
		index of the times to extract form the netcdf data.
	iBT : int
		CENTERED (i.e. unstaggered) indexes of desired bottom-top levels.
	iSN : int
		CENTERED (i.e. unstaggered) indexes of desired south-north coordinates.
	iWE : int
		CENTERED (i.e. unstaggered) indexes of desired weast-east coordinates.

	Returns
	-------
	out : ndarray
		numpy array with the extracted 3D or 4D data for the given time and coords dimensions

	'''
    if not varName in ncfid.variables.keys():
        raise Exception("ERROR: " + varName + " variable does not exist in the netcdf file")

    varObj = ncfid.variables.get(varName)

    ncDims = varObj.dimensions

    # it is assuming that first dimension is Time!
    logicSlices = [iTimes]
    stageredDim = -1

    for ii, iStr in enumerate(ncDims[1:]):

        if iStr.find('stag') > 0:

            stageredDim = ii + 1

            if iStr.startswith('bottom'):
                logicSlices.append(np.append(iBT, iBT[-1] + 1))
            elif iStr.startswith('south'):
                logicSlices.append(np.append(iSN, iSN[-1] + 1))
            elif iStr.startswith('west'):
                logicSlices.append(np.append(iWE, iWE[-1] + 1))

        else:

            if iStr.startswith('bottom'):
                logicSlices.append(iBT)
            elif iStr.startswith('south'):
                logicSlices.append(iSN)
            elif iStr.startswith('west'):
                logicSlices.append(iWE)

    # Extract and unstager the data
    varData = varObj[logicSlices]

    if (len(ncDims) == 3) and (stageredDim > 0):
        if stageredDim == 1:
            varData = (varData[:, 0:-1, :] + varData[:, 1:, :]) * 0.5
        elif stageredDim == 2:
            varData = (varData[:, :, 0:-1] + varData[:, :, 1:]) * 0.5

    elif (len(ncDims) == 4) and (stageredDim > 0):
        if stageredDim == 1:
            varData = (varData[:, 0:-1, :, :] + varData[:, 1:, :, :]) * 0.5
        elif stageredDim == 2:
            varData = (varData[:, :, 0:-1, :] + varData[:, :, 1:, :]) * 0.5
        elif stageredDim == 3:
            varData = (varData[:, :, :, 0:-1] + varData[:, :, :, 1:]) * 0.5

    return varData.data


def readAllvars(ncfid, varsWRF, iTimes, iBT, iSN, iWE):
    '''
	Function to loop over the same netcdf to extract several variables for
	a given set of times and subdomain.
	It is more efficient as it doesn't need to open and close the nc
	file several times

	Parameters
	----------
	ncfid : file identifier class
		netcdf file identifier..
	varsWRF : dict
		dictionary with the list of string of the variables aimed to extract from the netcdf.
	iTimes : int, logic
		index of the times to extract form the netcdf data.
	iBT : int
		CENTERED (i.e. unstaggered) indexes of desired bottom-top levels.
	iSN : int
		CENTERED (i.e. unstaggered) indexes of desired south-north coordinates.
	iWE : int
		CENTERED (i.e. unstaggered) indexes of desired weast-east coordinates.

	Returns
	-------
	dOut : dict
		dictionary with the variableName:nparray of the subset of data from the netcdf
		for each variable.

	'''

    # create output dictionary
    dOut = {}

    # loop over all variables
    for v2extract in varsWRF:
        # print('  extracting from netcdf: '+v2extract)

        if v2extract == 'L':
            # RMOL stands for the 1/ Obukhov length, thus we directly save it as L
            dOut['L'] = 1. / getVarEff(ncfid, 'RMOL', iTimes, iBT, iSN, iWE)
            dOut['L'][np.isinf(dOut['L'])] = np.nan
        elif v2extract == 'TENDENCIES':
            # coriolis as it is specified in nc file
            fc = getVarEff(ncfid, 'F', iTimes, iBT, iSN, iWE)
            fc = fc.mean()

            dOut['Vg'] = -(1 / fc) * getVarEff(ncfid, 'RU_TEND_PGF', iTimes, iBT, iSN, iWE)
            dOut['Ug'] = (1 / fc) * getVarEff(ncfid, 'RV_TEND_PGF', iTimes, iBT, iSN, iWE)
            dOut['UADV'] = (1 / fc) * getVarEff(ncfid, 'RU_TEND_ADV', iTimes, iBT, iSN, iWE)
            dOut['VADV'] = (1 / fc) * getVarEff(ncfid, 'RV_TEND_ADV', iTimes, iBT, iSN, iWE)
            dOut['POT_ADV'] = getVarEff(ncfid, 'T_TEND_ADV', iTimes, iBT, iSN, iWE)
        elif v2extract == 'Th':
            # The Potential temperature is extracted as the perturbation potential temperature + base state
            # temperature (t00
            T00 = ncfid.variables.get('T00')[iTimes]
            dOut['Th'] = getVarEff(ncfid, 'T', iTimes, iBT, iSN, iWE)
            for iT in range(len(T00)):
                dOut['Th'][iT, :, :, :] = dOut['Th'][iT, :, :, :] + T00[iT]
        elif v2extract == 'TKE':
            dOut['TKE'] = getVarEff(ncfid, 'TKE_PBL', iTimes, iBT, iSN, iWE)
        else:
            dOut[v2extract] = getVarEff(ncfid, v2extract, iTimes, iBT, iSN, iWE)

    # dummy matrix with zeros, to be added as surface values to 4D data
    zSfc = np.zeros((sum(iTimes), 1, len(iSN), len(iWE)))

    # add a layer at the surface for some variables
    for vN in dOut:
        if vN == 'Th':
            # for better consistency the skin temperature is added as sfc level
            tsk = getVarEff(ncfid, 'TSK', iTimes, iBT, iSN, iWE)
            dOut['Th'] = np.concatenate((np.reshape(tsk, zSfc.shape), dOut['Th']), axis=1)

        elif len(dOut[vN].shape) == 4:
            dOut[vN] = np.concatenate((zSfc, dOut[vN]), axis=1)

    return dOut


def get_nc_coordinates(ncfid, iTimes, iBT, iSN, iWE):
    '''
	Memory efficient function to extract the height (averaged in time) of WRF output
	The height is provided in 3D (i.e. bottom-top, north-south, west-east)

	Parameters
	----------
	ncfid : file id
		file of the netcdf.
	iTimes : int, logic
		index of the times to extract form the netcdf data.
	iBT : int
		CENTERED (i.e. unstaggered) indexes of desired bottom-top levels.
	iSN : int
		CENTERED (i.e. unstaggered) indexes of desired south-north coordinates.
	iWE : int
		CENTERED (i.e. unstaggered) indexes of desired weast-east coordinates.

	Returns
	-------
	out : ndarray
		numpy array with the height above sea level

	'''
    nz = len(iBT) + 1

    ltmp = getVarEff(ncfid, 'XLAT', iTimes, iBT, iSN, iWE).mean(axis=0)
    LAT = np.tile(ltmp, (nz, 1, 1))

    ltmp = getVarEff(ncfid, 'XLONG', iTimes, iBT, iSN, iWE).mean(axis=0)
    LON = np.tile(ltmp, (nz, 1, 1))

    PHT = (getVarEff(ncfid, 'PH', iTimes, iBT, iSN, iWE) + getVarEff(ncfid, 'PHB', iTimes, iBT, iSN, iWE)) / 9.81
    HGT = getVarEff(ncfid, 'HGT', iTimes, iBT, iSN, iWE)

    # time-average height
    zT = PHT.mean(axis=0)

    # time-average surface height
    hSfc = np.zeros((1, len(iSN), len(iWE)))
    hSfc[0, :, :] = HGT.mean(axis=0)

    Z = np.concatenate((hSfc, zT), axis=0)

    return LON, LAT, Z


def get_zagl(ncfid, iTimes, iBT, iSN, iWE):
    '''
	Memory efficient function to extract the height above ground of WRF output
	The height is provided in 4D (i.e. with time) dimensions and cell-centered

	Parameters
	----------
	ncfid : file id
		file of the netcdf.
	iTimes : int, logic
		index of the times to extract form the netcdf data.
	iBT : int
		CENTERED (i.e. unstaggered) indexes of desired bottom-top levels.
	iSN : int
		CENTERED (i.e. unstaggered) indexes of desired south-north coordinates.
	iWE : int
		CENTERED (i.e. unstaggered) indexes of desired weast-east coordinates.

	Returns
	-------
	out : ndarray
		numpy array with the height above ground
	'''

    # dummy matrix with zeros, to be added as value at the surface
    zSfc = np.zeros((sum(iTimes), 1, len(iSN), len(iWE)))

    PHT = (getVarEff(ncfid, 'PH', iTimes, iBT, iSN, iWE) + getVarEff(ncfid, 'PHB', iTimes, iBT, iSN, iWE)) / 9.81
    HGT = np.reshape(getVarEff(ncfid, 'HGT', iTimes, iBT, iSN, iWE), zSfc.shape)

    zaglT = PHT - np.repeat(HGT, PHT.shape[1], axis=1)

    return np.concatenate((zSfc, zaglT), axis=1)


def coriolis(lat):
    '''
	Compute the Coriolis frequency based on the latitude parameter

	Parameters
	----------
	lat : float
		latitude .

	Returns
	-------
	fc : float
		Coriolis frequency.
	'''
    # angular speed of the Earth [rad/s]
    omega = 7.2921159e-5

    return 2 * omega * np.sin(lat * np.pi / 180)


def getneighbours(Nxy, inear, jnear):
    '''
	Function to obtain the indexes of the neighours to a given central index
	from WRF in a given spatial box
	INPUS:
		Nxy = number of points to include in the spatial averaging box
		inear = given index in the west-east direction
		jnear = given index in the south-north direction
	'''
    if Nxy == 1:  # nearest grid point
        ixav = np.array([inear]).astype(int)
        iyav = np.array([jnear]).astype(int)

    elif Nxy == 2:  # four nearest grid points
        ixav = np.array([inear, inear + 1]).astype(int)
        iyav = np.array([jnear, jnear + 1]).astype(int)
    else:
        if Nxy % 2 == 1:  # Nxy (odd) nearest points
            ixav = np.arange(inear - 0.5 * (Nxy - 1), inear + 0.5 * (Nxy - 1) + 1).astype(int)
            iyav = np.arange(jnear - 0.5 * (Nxy - 1), jnear + 0.5 * (Nxy - 1) + 1).astype(int)
        else:
            ixav = np.arange(inear - 0.5 * Nxy + 1, inear + 0.5 * Nxy + 1).astype(int)
            iyav = np.arange(jnear - 0.5 * Nxy + 1, jnear + 0.5 * Nxy + 1).astype(int)

    return ixav, iyav


def get_index_of_subset_domain(ncfid, lat_s, lon_s, L=None):
    lat = ncfid.variables.get('XLAT')[0, :, 0]
    lon = ncfid.variables.get('XLONG')[0, 0, :]

    # makes sure central box coordinates lie inside wrf domain box
    if (min(abs(lat - lat_s)) > max(np.diff(lat))) | (min(abs(lon - lon_s)) > max(np.diff(lon))):
        raise Exception("ERROR: lat | lon chosen is outside wrf domain box")

    # get the grid-spacing assuming dx=dy
    dxy = ncfid.getncattr('DX')

    # Extract only a box of LxL of all wrf domain to make the process more memory efficient
    # the selection is made by the grid-spacing of wrf output
    if L is None:
        if dxy >= 9e3:
            L = 55e3
        elif dxy >= 1e3:
            L = 10e3
        else:
            L = 1e3

    print(' Only data in a box of ' + str(L) + 'x' + str(L) + ' is extracted as subset')
    # number of points to include in the spatial sampling box.
    Nxy = int(L / dxy) + 1

    # number of grid centered points in the bottom-top (vertical) dimension
    nz = extract_dim(ncfid, 'bottom_top')
    # get indexes of the bottom-top levels to extract (so far is all levels)
    iBT = np.arange(0, nz)

    # get indexes of nearest coordinate to the given point
    inear_WE, inear_SN = ll_to_xy(ncfid, lat_s, lon_s)

    # get indexes of nearest + Nxy points to the given point
    iWE, iSN = getneighbours(Nxy, inear_WE, inear_SN)

    return iBT, iSN, iWE


class nc_results:
    '''
	Class (lean version) to query any given point or set of points in from a previously
	created WRF dictionary

'''

    def __init__(self, timesSim, zagl):
        '''
		Initialize the class
		'''

        self.times = timesSim

        # Reference Time for the netcdf files
        self.referenceDate = np.datetime64(datetime.datetime(1970, 1, 1, 0, 0, 0))

        # convert time stamp to "seconds since <dateRef>"
        self.seconds = pd.Series(self.time - self.referenceDate).dt.total_seconds().values

        self.heights = np.nanmean(zagl, axis=(0, 2, 3))

        self.nt, self.nz, self.ny, self.nx = zagl.shape

    def to_timeseries(self, vDict, X, Y, Z, qCoords):
        '''
		Function that writes the variables extracted and procesed from WRF into a new nc file of time series stlye
		Parameters
		----------
		vDict : dictionary
			Python dictionary with the variables extracted from WRF

		Returns
		-------
		out : ndarray
				numpy array with the height above ground

		'''

        index = self._find_nearest_index(X, Y, Z, qCoords)

        self.timeseries = {}
        for iVar in vDict:
            print('  variable: ' + iVar)

            vtmp = vDict[iVar]
            if len(vDict[iVar].shape) == 3:
                vtmp = np.reshape(vtmp, (self.nt, 1, self.ny, self.nx))
                vtmp = np.tile(vtmp, (1, self.nz, 1, 1))

            vtmp = np.reshape(vtmp, (self.nt, self.nz * self.ny * self.nx))

            self.timeseries[iVar] = pd.DataFrame(vtmp[:, index], index=self.time)
            self.timeseries[iVar].name = iVar

    def _find_nearest_index(self, X, Y, Z, qCoords):

        # avoid creating the tree every time the function is called
        if not hasattr(self, 'tree'):
            coords = np.column_stack((X.flatten(), Y.flatten(), Z.flatten()))
            self.tree = cKDTree(coords)

        dist, index = self.tree.query(qCoords)

        return index

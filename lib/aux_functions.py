#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Set of auxliary functions for the WRF_2_SCD script
'''
import numpy as np

def getVarEff(ncfid, varName, iTimes, iBT, iSN, iWE):
	'''
	Memmory efficient method of getting the nc data and destagered (if
	that is the case ) into a centered grid.

	INPUTS:
		nctif   = netcdf file identifier
	    varName = name of the variale iTimes, iBT, iSN, iSW
		          CENTERED (i.e. unstaggered) indexes of desired
		          bottom-up, south-north and west-east dimensions respe
	'''
	varObj = ncfid.variables.get(varName)

	ncDims = varObj.dimensions

	# it is assuming that first dimension is Time!
	logicSlices = [iTimes]
	stageredDim = -1

	for ii,iStr in enumerate(ncDims[1:]):

		if iStr.find('stag') > 0:

			stageredDim = ii+1

			if iStr.startswith('bottom'):
				logicSlices.append( np.append(iBT , iBT[-1]+1 ) )
			elif iStr.startswith('south'):
				logicSlices.append( np.append(iSN , iSN[-1]+1) )
			elif iStr.startswith('west'):
				logicSlices.append( np.append(iWE , iWE[-1]+1) )

		else:

			if iStr.startswith('bottom'):
				logicSlices.append( iBT )
			elif iStr.startswith('south'):
				logicSlices.append( iSN )
			elif iStr.startswith('west'):
				logicSlices.append( iWE )

	# Extgract and unstager data
	varData = varObj[logicSlices]

	if (len(ncDims) == 3) and (stageredDim > 0):
		if stageredDim == 1:
			varData = (varData[:,0:-1,:] + varData[:,1:,:]) * 0.5
		elif stageredDim == 2:
			varData = (varData[:,:,0:-1] + varData[:,:,1:]) * 0.5

	elif (len(ncDims) == 4) and (stageredDim > 0):
		if stageredDim == 1:
			varData = (varData[:,0:-1,:,:] + varData[:,1:,:,:]) * 0.5
		elif stageredDim == 2:
			varData = (varData[:,:,0:-1,:] + varData[:,:,1:,:]) * 0.5
		elif stageredDim == 3:
			varData = (varData[:,:,:,0:-1] + varData[:,:,:,1:]) * 0.5

	return varData.data


def readAllvars(ncfid, varsWRF,iTimes, iBT, iSN, iWE, lat):
	'''
	Auxiliary function to extract the WRF results of the variables set in varsWRF list
	'''
	# output dictionary
	dOut={}
	for v2extract in varsWRF:
		print('  extracting from netcdf: '+v2extract)

		if (v2extract == 'RMOL'):
			# RMOL stands for the inverso of the Obukhov length, thus we directly save it as L
			dOut['L'] = 1./getVarEff(ncfid, 'RMOL', iTimes, iBT, iSN, iWE)
			dOut['L'][np.isinf(dOut['L'])] = np.nan
		elif (v2extract == 'TENDENCIES'):
			fc = coriolis(lat)
			dOut['Vg']      = -(1/fc) * getVarEff(ncfid, 'RU_TEND_PGF', iTimes, iBT, iSN, iWE)
			dOut['Ug']      =  (1/fc) * getVarEff(ncfid, 'RV_TEND_PGF', iTimes, iBT, iSN, iWE)
			dOut['UADV']    =  (1/fc) * getVarEff(ncfid, 'RU_TEND_ADV', iTimes, iBT, iSN, iWE)
			dOut['VADV']    =  (1/fc) * getVarEff(ncfid, 'RV_TEND_ADV', iTimes, iBT, iSN, iWE)
			dOut['POT_ADV'] =           getVarEff(ncfid, 'T_TEND_ADV',  iTimes, iBT, iSN, iWE)
		elif (v2extract == 'POT'):
			# The Potential temperature is extracted as the perturbation potential temperature + base state temperature (t00
			T00 = ncfid.variables.get('T00')[iTimes]
			dOut['POT'] = getVarEff(ncfid, 'T', iTimes, iBT, iSN, iWE)
			for iT in range(len(T00)):
				dOut['POT'][iT,:,:,:] = dOut['POT'][iT,:,:,:] + T00[iT]
		else:
			dOut[v2extract] = getVarEff(ncfid, v2extract, iTimes, iBT, iSN, iWE)

	#dummy matrix with zeros, to be added as surface values to 4D data
	zSfc = np.zeros((sum(iTimes),1,len(iWE),len(iSN)))

	# add a layer at the surface for some variables
	for vN in dOut:
		if (vN == 'POT'):
			# for better consistency the skin temperature is added as sfc level
			tsk = getVarEff(ncfid,'TSK', iTimes, iBT, iSN, iWE)
			dOut['POT'] = np.concatenate((np.reshape(tsk, zSfc.shape), dOut['POT']), axis=1)

		elif (len(dOut[vN].shape) == 4):
			dOut[vN] = np.concatenate((zSfc, dOut[vN]), axis=1)

	return dOut

def get_zagl(ncfid, iTimes, iBT, iSN, iWE):
	'''
	memory efficient function to extract the height above ground of WRF output
	'''
	#dummy matrix with zeros, to be added as value at the surface
	zSfc = np.zeros((sum(iTimes),1,len(iWE),len(iSN)))

	PHT = ( getVarEff(ncfid,'PH',iTimes,iBT,iSN,iWE) + getVarEff(ncfid,'PHB',iTimes,iBT,iSN,iWE) )/9.81
	HGT = getVarEff(ncfid,'HGT',iTimes,iBT,iSN,iWE)

	zaglT = np.empty(PHT.shape, np.float32)
	for jj in range(PHT.shape[1]):
		zaglT[:,jj,:,:] = PHT[:,jj,:,:] - HGT

	return np.concatenate((zSfc,zaglT),axis=1)

def coriolis(lat):
	'''
	Compute the Coriolis frequency based on the latitude parameter
	'''
	omega  = 7.2921159e-5    # angular speed of the Earth [rad/s]
	fc = 2 * omega * np.sin(lat * np.pi / 180)
	return fc

def getneighbours(Nav,inear,jnear):
	'''
	Function to obtain the indexes of the neighours to a given central index
	from WRF in a given spatial box
	INPUS:
		Nav = number of points to include in the spatial averaging box
		inear = given index in the west-east direction
		jnear = given index in the south-north direction
	'''
	if Nav == 1:  # nearest grid point
		ixav = np.array([inear]).astype(int)
		iyav = np.array([jnear]).astype(int)

	elif Nav == 2:  # four nearest grid points
		ixav = np.array([inear, inear + 1]).astype(int)
		iyav = np.array([jnear, jnear + 1]).astype(int)
	else:
		if Nav % 2 == 1:  # Nav (odd) nearest points
			ixav = np.arange(inear - 0.5 * (Nav - 1), inear + 0.5 * (Nav - 1) + 1).astype(int)
			iyav = np.arange(jnear - 0.5 * (Nav - 1), jnear + 0.5 * (Nav - 1) + 1).astype(int)
		else:
			ixav = np.arange(inear - 0.5 * Nav + 1, inear + 0.5 * Nav + 1).astype(int)
			iyav = np.arange(jnear - 0.5 * Nav + 1, jnear + 0.5 * Nav + 1).astype(int)

	return ixav,iyav
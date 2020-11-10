import pandas as pd
import numpy as np
import datetime
import netCDF4
import utm
import xarray as xr
import math
import rasterio as rio
import matplotlib.pyplot as plt 
from matplotlib.colors import LightSource

from lib.variables_dictionary.variables import Variables
from lib.variables_dictionary.variables import nc_global_attributes_from_yaml

def read_sim(filename):
    print(filename)
    M = xr.open_dataset(filename)
    try: 
        U = M.eastward_wind
        V = M.northward_wind
        S = (U**2 + V**2)**0.5
        WD = (270-np.rad2deg(np.arctan2(V,U)))%360
        M = M.assign(wind_speed = S, wind_direction = WD)
    except AttributeError:
        print('No U and V fields found')
        
    return M

def read_obs(filename):
    M = xr.open_dataset(filename)
    #M = M.where(M <= 360) # nan when is greater than 360
    M = M.rename({'wind_from_direction': 'wind_direction'})
    M = M.rename({'specific_turbulent_kinetic_energy': 'turbulent_kinetic_energy'})
    M = M.rename({'specific_upward_sensible_heat_flux_in_air': 'heat_flux'})
    M = M.interpolate_na(dim='height') # Fill nans with interpolated values between observational levels
    return M

def vector_mean_std_wind_direction(WD):
    U = -np.sin(2*np.pi*WD/360)
    V = -np.cos(2*np.pi*WD/360)
    Umean = np.mean(U)
    Vmean = np.mean(V)
    WDmean = 180 + math.atan2(Umean,Vmean)*180/np.pi
    eps = np.sqrt(1 - (Umean**2 + Vmean**2))
    WDstd = math.asin(eps)*(1 + (2./np.sqrt(3.) - 1)*eps**3)*180/np.pi # Yamartino (1984)
    return WDmean, WDstd

def events_qois_vs_masts_table(events,masts_obs,height):
    qois = ['wind_speed', 'wind_direction', 'turbulence_intensity',
                            'wind_shear','stability'] 
    events_name = list(events.keys())

    index = pd.MultiIndex.from_product([list(events.keys()), qois],
                                       names=['event', 'QoI'])
    columns = pd.MultiIndex.from_product([masts_obs.id.values.tolist(), ['mean', 'std']],
                                         names=['mast', ''])

    # create the DataFrame
    df = pd.DataFrame(np.random.randn(len(events_name)*len(qois), len(masts_obs.id)*2), index=index, columns=columns)
    idx = pd.IndexSlice

    for e in events_name:
        for q in qois:
            for m in masts_obs.id.values.tolist():
                if q == 'wind_speed' or q == 'turbulence_intensity':
                    df.loc[idx[e, q], idx[m, 'mean']] = masts_obs[q].sel(id = m, height = height, 
                                                                         time = slice(events[e][0],events[e][1])).mean()
                    df.loc[idx[e, q], idx[m, 'std']] = masts_obs[q].sel(id = m, height = height, 
                                                                        time = slice(events[e][0],events[e][1])).std()
                elif q == 'wind_direction':
                    WD = masts_obs[q].sel(id = m, height = height, time = slice(events[e][0],events[e][1]))
                    df.loc[idx[e, q], idx[m, 'mean']], df.loc[idx[e, q], idx[m, 'std']] = vector_mean_std_wind_direction(WD)
                elif q == 'wind_shear':
                    df.loc[idx[e, q], idx[m, 'mean']] = masts_obs[q].sel(id = m, 
                                                                        time = slice(events[e][0],events[e][1])).mean()
                    df.loc[idx[e, q], idx[m, 'std']] = masts_obs[q].sel(id = m, 
                                                                        time = slice(events[e][0],events[e][1])).std()
                elif q == 'stability':
                    df.loc[idx[e, q], idx[m, 'mean']] = masts_obs[q].sel(id = m, height = 10, 
                                                                         time = slice(events[e][0],events[e][1])).mean()
                    df.loc[idx[e, q], idx[m, 'std']] = masts_obs[q].sel(id = m, height = 10, 
                                                                        time = slice(events[e][0],events[e][1])).std()
    return df 

                
def mast_sims_vs_obs_timeseries_plot(mast, h, masts_obs, masts_sim, sims, datefrom, dateto, events):
    fig, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(6,1,figsize = (14,14), sharex = True)
    masts_obs.wind_speed.sel(id = mast, height = h).plot(x = 'time', label = 'obs', color = 'k', ax = ax1)
    masts_obs.wind_direction.sel(id = mast, height = h).plot(x = 'time', label = 'obs', color = 'k', ax = ax2)
    masts_obs.turbulence_intensity.sel(id = mast, height = h).plot(x = 'time', label = 'obs', color = 'k', ax = ax3)
    masts_obs.wind_shear.sel(id = mast).plot(x = 'time', label = 'obs', color = 'k', ax = ax4)
    masts_obs.stability.sel(id = mast, height = 10).plot(x = 'time', label = 'obs', color = 'k', ax = ax5)
    masts_obs.heat_flux.sel(id = mast, height = 10).plot(x = 'time', label = 'obs', color = 'k', ax = ax6)
    for i_sim in range (0,len(masts_sim)):
        masts_sim[i_sim].wind_speed.sel(id = mast).interp(height= h).plot(x = 'time', label = sims['ID'][i_sim], ax = ax1)
        masts_sim[i_sim].wind_direction.sel(id = mast).interp(height= h).plot(x = 'time', label = sims['ID'][i_sim], ax = ax2)
        masts_sim[i_sim].turbulence_intensity.sel(id = mast).interp(height= h).plot(x = 'time', label = sims['ID'][i_sim], ax = ax3)
        masts_sim[i_sim].wind_shear.sel(id = mast).plot(x = 'time', label = sims['ID'][i_sim], ax = ax4)
    color_events = {'neutral': 'silver', 'unstable': 'salmon','stable': 'lightblue', 'very stable': 'deepskyblue'}
    for e in events:
        ax1.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
        ax2.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
        ax3.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
        ax4.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
        ax5.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
        ax6.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
    ax1.set_xlim([datefrom, dateto])
    ax1.legend(bbox_to_anchor=(1.13, 1))
    ax1.grid(); ax1.set_xlabel(''); ax1.set_ylabel(r'wind speed [$m s^{-1}$]')
    ax2.grid(); ax2.set_xlabel(''); ax2.set_ylabel(r'wind direction [º]')
    ax3.grid(); ax3.set_xlabel(''); ax3.set_ylabel(r'turbulence intensity'); ax3.set_ylim([0,0.6])
    ax4.grid(); ax4.set_xlabel(''); ax4.set_ylabel(r'wind shear $\alpha(80/40)$'); ax4.set_ylim([-0.5,0.5])
    ax5.grid(); ax5.set_xlabel(''); ax5.set_ylabel(r'stability $z/L$ ($z$=10m)'); ax5.set_ylim([-1,1])
    ax6.grid(); ax6.set_xlabel(''); ax6.set_ylabel(r'heat flux [$m K s^{-1}$]');# ax6.set_ylim([-1,1])
    return [ax1, ax2, ax3, ax4, ax5, ax6]

def compare_masts_timeseries_plot(mast, h, masts_obs, datefrom, dateto, events):
    fig, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(6,1,figsize = (14,14), sharex = True)
    masts_obs.wind_speed.sel(id = mast, height = h).plot(x = 'time', hue = 'id', ax = ax1)
    masts_obs.wind_direction.sel(id = mast, height = h).plot(x = 'time', hue = 'id', ax = ax2)
    masts_obs.turbulence_intensity.sel(id = mast, height = h).plot(x = 'time', hue = 'id', ax = ax3)
    masts_obs.wind_shear.sel(id = mast).plot(x = 'time', hue = 'id', ax = ax4)
    masts_obs.stability.sel(id = mast, height = 10.).plot(x = 'time', hue = 'id', ax = ax5)
    masts_obs.heat_flux.sel(id = mast, height = 10.).plot(x = 'time', hue = 'id', ax = ax6)
    color_events = {'neutral': 'silver', 'unstable': 'salmon','stable': 'lightblue', 'very stable': 'deepskyblue'}
    for e in events:
        ax1.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
        ax2.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
        ax3.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
        ax4.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
        ax5.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
        ax6.axvspan(events[e][0], events[e][1], alpha=0.5, color=color_events[e])
    ax1.set_xlim([datefrom, dateto])
    #ax1.legend(bbox_to_anchor=(1.13, 1))
    ax2.get_legend().remove(); ax3.get_legend().remove(); ax4.get_legend().remove(); ax5.get_legend().remove()
    ax1.grid(); ax1.set_xlabel(''); ax1.set_ylabel(r'wind speed [$m s^{-1}$]')
    ax2.grid(); ax2.set_xlabel(''); ax2.set_ylabel(r'wind direction [º]')
    ax3.grid(); ax3.set_xlabel(''); ax3.set_ylabel(r'turbulence intensity [-]'); ax3.set_ylim([0,0.6])
    ax4.grid(); ax4.set_xlabel(''); ax4.set_ylabel(r'wind shear $\alpha(80/40)$'); ax4.set_ylim([-0.5,0.5])
    ax5.grid(); ax5.set_xlabel(''); ax5.set_ylabel(r'stability $z/L$ ($z$=10m)'); ax5.set_ylim([-1,1])
    ax6.grid(); ax6.set_xlabel(''); ax6.set_ylabel(r'heat flux [$m K s^{-1}$]');# ax6.set_ylim([-1,1])
    return [ax1, ax2, ax3, ax4]

def masts_sims_vs_obs_profiles_plot(event, masts_obs, masts_sim, sims):
    fig, axes = plt.subplots(2,4,figsize = (14,8), sharey = True, sharex = True)
    masts = masts_obs.coords['id'].values.tolist()
    for i, mast in enumerate(masts):
        index = np.unravel_index(i,(2,4))
        ax = axes[index]
        for i_sim in range (0,len(masts_sim)):
            h_sim = masts_sim[i_sim].wind_speed.sel(time = slice(event[0],event[1]), id = mast).mean(dim = 'time', skipna = True).plot(y = 'height', ax = axes[index], label = sims['ID'][i_sim])
        h_obs = masts_obs.wind_speed.sel(time = slice(event[0],event[1]), id = mast).mean(dim = 'time', skipna = True).plot(y = 'height', ax = axes[index], label = 'obs',
                                                        marker = 'o', linestyle='none', color = 'silver')
        ax.set_ylabel(''); ax.set_xlabel(''); 
        ax.set_ylim(1,1000)
        ax.grid()
    plt.yscale('symlog')
    axes[(0,0)].set_ylabel('z [m]'); axes[(1,0)].set_ylabel('z [m]')
    axes[(1,0)].set_xlabel('wind speed [$m s^{-1}$]'); axes[(1,1)].set_xlabel('wind speed [$m s^{-1}$]'); axes[(1,2)].set_xlabel('wind speed [$m s^{-1}$]')
    axes[(1,3)].axis('off')
    axes[(1,2)].legend(bbox_to_anchor=(1.13, 1))
    return axes

def Ztransect_sims_vs_obs_plot(t, Ztransect_obs, Ztransect_sim, masts, Ztransect):
    fig, (ax1,ax2) = plt.subplots(2,1, figsize = (8,6), sharex = True)
    Zprofile = Ztransect_plot(masts, Ztransect, ax2)
    ax2.set_title('')
    #for i_sim in range (0,n_sim):
        #Ztransect_sim[i_sim].wind_speed.sel(height = h).plot(x = 'id', label = sims['ID'][i_sim], ax = ax1)
    #ax1.legend(bbox_to_anchor=(1, 1))
    masts_inZ = [] # index of Z_transect position nearest to each mast
    for i, row in masts.iterrows():
        d = np.sqrt((Ztransect['x'] - masts['x'][i])**2 + (Ztransect['y'] - masts['y'][i])**2)
        masts_inZ.append(d[d == d.min()].index[0])
    for x in masts_inZ:
        ax1.axvline(x, color = 'silver', linestyle = '--', zorder = 0)
    return [ax1, ax2]

def basemap_plot(src, masts, Ztransect, ref, ax, coord = 'utm'):
    # Add overviews to raster to plot faster at lower resolution (https://rasterio.readthedocs.io/en/latest/topics/overviews.html)
    #from rasterio.enums import Resampling
    #factors = [2, 4, 8, 16]
    #dst = rio.open('./inputs/DTM_Alaiz_2m.tif', 'r+')
    #dst.build_overviews(factors, Resampling.average)
    #dst.update_tags(ns='rio_overview', resampling='average')
    #dst.close()
    A_ind = Ztransect['Name'].str.contains('A') # Tajonar ridge scan
    B_ind = Ztransect['Name'].str.contains('B') # Elortz valley scan
    C_ind = Ztransect['Name'].str.contains('C') # Alaiz ridge scan
    oview = src.overviews(1)[2] # choose overview (0 is largest, -1 is the smallest)
    topo = src.read(1, out_shape=(1, int(src.height // oview), int(src.width // oview)))
    if coord == 'xy':
        spatial_extent = [src.bounds.left - ref[0], src.bounds.right - ref[0], src.bounds.bottom - ref[1], src.bounds.top - ref[1]]
    else:
        spatial_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
        
    topo_ma = np.ma.masked_where(topo == 0 , topo, copy=True) 
#    ls = LightSource(azdeg=315, altdeg=60)
#    rgb = ls.shade(topo_ma, cmap=plt.cm.terrain, blend_mode='overlay') 
#    h_topo = ax.imshow(rgb, extent=spatial_extent, vmin=400, vmax=1200)
    h_topo = ax.imshow(topo_ma, cmap = plt.cm.terrain, extent=spatial_extent, vmin=300, vmax=1200)
    if coord == 'xy':
        h_masts = ax.scatter(masts['x'], masts['y'], s = 10, marker='s', c='k', label = 'Masts')
        h_A = ax.scatter(Ztransect[A_ind]['x'], Ztransect[A_ind]['y'], s = 2, marker='.', c='blue', label = 'A-transect')
        h_B = ax.scatter(Ztransect[B_ind]['x'], Ztransect[B_ind]['y'], s = 2, marker='.', c='black', label = 'B-transect')
        h_C = ax.scatter(Ztransect[C_ind]['x'], Ztransect[C_ind]['y'], s = 2, marker='.', c='red', label = 'C-transect')
        for i, txt in enumerate(masts['Name']):
            ax.annotate(txt, (masts['x'][i]+50, masts['y'][i]+50))   
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
    else:
        h_masts = ax.scatter(masts['easting[m]'], masts['northing[m]'], s = 10, marker='s', c='k', label = 'Masts')
        h_A = ax.scatter(Ztransect[A_ind]['easting[m]'], Ztransect[A_ind]['northing[m]'], s = 2, marker='.', c='blue', label = 'A-transect')
        h_B = ax.scatter(Ztransect[B_ind]['easting[m]'], Ztransect[B_ind]['northing[m]'], s = 2, marker='.', c='black', label = 'B-transect')
        h_C = ax.scatter(Ztransect[C_ind]['easting[m]'], Ztransect[C_ind]['northing[m]'], s = 2, marker='.', c='red', label = 'C-transect')
        for i, txt in enumerate(masts['Name']):
            ax.annotate(txt, (masts['easting[m]'][i]+50, masts['northing[m]'][i]+50))   
        ax.set_xlabel('Easting [m]')
        ax.set_ylabel('Northing [m]')      
            
    ax.set_title('ALEX17 sites')

    ax.legend(handles = [h_masts,h_A,h_B,h_C])
    plt.colorbar(h_topo, ax = ax)
    return [h_masts,h_A,h_B,h_C]

    
def Ztransect_plot(masts, Ztransect, ax):
    A_ind = Ztransect['Name'].str.contains('A') # Tajonar ridge scan
    B_ind = Ztransect['Name'].str.contains('B') # Elortz valley scan
    C_ind = Ztransect['Name'].str.contains('C') # Alaiz ridge scan
    h_topoZ = (Ztransect['z']-125.).plot.area(stacked=False, color = 'lightgrey', alpha = 1, ax = ax)
    h_topoA = Ztransect[A_ind]['z'].plot.line(style = '.', color = 'blue', ms = 3, ax = ax)
    h_topoB = Ztransect[B_ind]['z'].plot.line(style = '.', color = 'black', ms = 3, ax = ax)
    h_topoC = Ztransect[C_ind]['z'].plot.line(style = '.', color = 'red', ms = 3, ax = ax)
    masts_inZ = [] # index of Ztransect position nearest to each mast
    for i, row in masts.iterrows():
        d = np.sqrt((Ztransect['x'] - masts['x'][i])**2 + (Ztransect['y'] - masts['y'][i])**2)
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
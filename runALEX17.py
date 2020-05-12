import lib.alex17_functions
import lib.WrfReader
import numpy as np
import utm
import time

# ['standard_name', 'wrf_short_name']
variables_to_write = np.array([['eastward_wind', 'U'],
                               ['northward_wind', 'V'],
                               ['upward_air_velocity', 'W'],
                               ['air_potential_temperature', 'Th'],
                               ['specific_turbulent_kinetic_energy', 'TKE']])

# SETUP THE WRF FILE READER
# a class for reading the simulation results
wrf_files_path = './inputs/'
wrf_domain_number = 3
dateLimits = ["2018-09-30 00:00:00", "2018-10-04 00:00:00"]
variables_to_extract = variables_to_write[:, 1]
(lat, lon) = utm.to_latlon(612000, 4726000, 30, 'T')
domain_of_interest = [lat, lon, 30e3]

wrf_inp = lib.WrfReader.WrfReader(wrf_files_path, wrf_domain_number, dateLimits, variables_to_extract,
                                  subset_of_wrfDomain=domain_of_interest)

# DEFINE THE GET FUNCTIONS TO BE USED BY THE SCRIPT
f_get_column = wrf_inp.get_column  # (lat, lon) -> (time, height, variables)
f_get_point = wrf_inp.get_point  # (lat, lon, height) -> (time, variables)

# TEST IF IT WORKS
# f_get_point(42.695, -1.558, 100)
# f_get_column(42.695, -1.558)

utm_zone_number = 30
utm_zone_letter = 'T'
(lat, lon) = lib.alex17_functions.utm.to_latlon(615100, 4730900, utm_zone_number, utm_zone_letter)
f_get_point(lat, lon, 40)

sim_id = 'alex17_00a'
output_folder = './outputs/'
# ------------------------------- 1. simID_zt.nc:
# Vertical profiles extending 1km from the surface at the met mast locations (see Validation Data for UTM coordinates);
file_path = output_folder + sim_id + '_masts.nc'

start_time = time.time()
lib.alex17_functions.create_alex17_file_1(file_path, f_get_column, variables_to_write)
print("--- %s seconds -1--" % (time.time() - start_time))

# ------------------------------- 2. simID_xyt.nc: Horizontal planes at 125m and 40m a.g.l. with horizontal
# resolution of 100 m, covering the area [612–622; 4726–4736] km;
start_time = time.time()
file_path = output_folder + sim_id + '_box.nc'
#lib.alex17_functions.create_alex17_file_2(file_path, f_get_point, variables_to_write)
print("--- %s seconds -2--" % (time.time() - start_time))

# ------------------------------- 3. simID_yzt.nc: Vertical planes along the Transect Line (see Validation Data for
# UTM coordinates) extending 1 km in the vertical and until [4726–4736] km in the N-S direction
start_time = time.time()
file_path = output_folder + sim_id + '_Ztransect.nc'
lib.alex17_functions.create_alex17_file_3(file_path, f_get_column, variables_to_write)
print("--- %s seconds -3--" % (time.time() - start_time))

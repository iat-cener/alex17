import lib.alex17_functions


# ['standard_name', 'short_name']
variables_to_write = [['eastward_wind', 'U'],
                      ['northward_wind', 'V'],
                      ['upward_air_velocity', 'W'],
                      ['air_potential_temperature', 'Th'],
                      ['TKE', 'TKE']]

# a class for reading the simulation results
wrf_inp = lib.alex17_functions.WrfReader(variables_to_write)
f_get_column = wrf_inp.get_column  # (lat, lon) -> (time, height, variables)
f_get_point = wrf_inp.get_point  # (lat, lon, height) -> (time, variables)
# f_get_point(42.695, -1.558, 100)
# f_get_column(42.695, -1.558)

utm_zone_number = 30
utm_zone_letter = 'T'
(lat, lon) = lib.alex17_functions.utm.to_latlon(615100, 4730900, utm_zone_number, utm_zone_letter)
f_get_point(lat, lon, 40)


sim_id = 'wrf_demo'
output_folder = './out/'
# ------------------------------- 1. simID_zt.nc:
# Vertical profiles extending 1km from the surface at the met mast locations (see Validation Data for UTM coordinates);
file_path = output_folder + sim_id + '_zt.nc'
lib.alex17_functions.create_alex17_file_1(file_path, f_get_column, variables_to_write)

# ------------------------------- 2. simID_xyt.nc: Horizontal planes at 125m and 40m a.g.l. with horizontal
# resolution of 100 m, covering the area [612–622; 4726–4736] km;
file_path = output_folder + sim_id + '_hxyt.nc'
lib.alex17_functions.create_alex17_file_2(file_path, f_get_point, variables_to_write)

# ------------------------------- 3. simID_yzt.nc: Vertical planes along the Transect Line (see Validation Data for
# UTM coordinates) extending 1 km in the vertical and until [4726–4736] km in the N-S direction
file_path = output_folder + sim_id + '_yzt.nc'
lib.alex17_functions.create_alex17_file_3(file_path, f_get_column, variables_to_write)





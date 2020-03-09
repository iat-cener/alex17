# Variables dictionary for wind energy
The dictionary allows definiton of the best practice of variable naming in the wind energy community. It is intended to be a dictionarry allowing translation betwean different naming conventions used by the community and we welcome colaborators willing to extend it.

# Not yet-another naming convention!
While the dictionary allows definition of a "default" name for variables, it is not intended to define a new naming convention. In general, the default variable name used come from the CF conventions. If unaviable, other existing standards are used. Non standard names are only used when no other option is available.

# Use
```python
#import the class
import variables
#intialise the dictionary - it will test the .json files accordin to the scheema and build the database
var_dict = variables.Variables()
#fetch the data of a variable
metadata = var_dict.lookup('time')
```

# More than a dictionary
The dictionarry includes additional metadata necesarry for use with NetCDF files
```python
# def nc_create(self, output_dataset, name, dimensions, standard_name = ""):
var = var_dict.nc_create(output_dataset, 'time', ('time',))
```

# Colaboration
The dictionary is ment to be multidisciplinary, community driven and build while being used. To colaborate, write an email to pgancarski@cener.com

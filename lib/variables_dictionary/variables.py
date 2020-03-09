import jsonschema
import simplejson as json
# import datetime
import csv
import os
# from StringIO import StringIO
from io import StringIO
from sys import exit
from yaml import safe_load


def nc_global_attributes_from_yaml(nc_file, file_path):
    try:
        with open(file_path, 'r') as stream:
            config = safe_load(stream)
            for key, value in config.items():
                setattr(nc_file, key, value)
    except FileNotFoundError as e:
        print('bad_config_path', file_path)
    except Exception as e:
        print('bad_config_formatting', str(e))


class Variables:
    """
    Database of variables and related metadata
    """

    def __init__(self):
        path = os.path.abspath(__file__)  # i.e. /path/to/dir/foobar.py
        self.dictionaries_path = os.path.split(path)[0]  # i.e. /path/to/dir/

        ## Read the schema file 
        with open(os.path.join(self.dictionaries_path, 'schema.json'), 'r') as f:
            schema_file = f.read()
        schema = json.loads(schema_file)

        dictionaries = ["IEC_Met.json", "CENER.json"]
        self.dictionary = []

        for file_name in dictionaries:
            ## Read variables definitions
            with open(os.path.join(self.dictionaries_path, file_name), 'r') as f:
                json_file = f.read()
            variables_json_obj = json.loads(json_file)

            # Validate the variables definitions
            jsonschema.validate(variables_json_obj, schema)

            # Add to the global dictionary
            self.dictionary = self.dictionary + variables_json_obj

    def lookup(self, variable_name):
        """
        Finds a variable metadata based on its standardised name
        TODO Need an index for variables names to get results in O(1)
        """
        res = {}
        for obj in self.dictionary:
            if obj['name']['default'] == variable_name:
                res = obj
                break

        if len(res) == 0:
            print("WARNING: Unrecognised variable name " + variable_name)
            print(
                "If using for the first time, then you should add the metadata of the variable to an appropriate "
                ".json file in the lidaco/variables folder.")
            print("It is YOUR responsibility to do that.")
            res = self.lookup("undefined")
            res['name']['default'] == variable_name
        return res

    def nc_create_dimension(self, output_dataset, var_name, data, standard_name=""):
        output_dataset.createDimension(var_name, len(data))
        dim = self.nc_create(output_dataset, var_name, (var_name,), standard_name)
        dim[:] = data
        return dim

    def nc_create(self, output_dataset, name, dimensions, standard_name=""):
        """
        Creates a variable based on a standardised variable name in an nc object 
        and returns a reference to that variable.
        TODO Choice of naming conventions for variables names
        TODO copy other properties
        """
        if standard_name == "":
            standard_name = name

        metadata = self.lookup(name)
        var = output_dataset.createVariable(standard_name, metadata["netcdf"]["var_type"], dimensions, zlib=True)
        var.units = metadata["units"]
        var.long_name = metadata["name"]["default"]
        var.comment = metadata["description"]
        # var.accuracy = metadata["accuracy"]
        # var.accuracy_info = metadata["accuracy_info"]

        ### add other, variable specific attributes
        if metadata["netcdf"]["other"] != "":
            ## use csv reader as it handles well quotations
            other_args = csv.reader(StringIO(metadata["netcdf"]["other"]), delimiter=',')
            # print(other_args)
            for row in other_args:
                for arg in row:
                    # print(arg)
                    key, val = arg.split('=')
                    key = key.strip()
                    val = val.strip().strip('"')
                    setattr(var, key, val)
        return var

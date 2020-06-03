# alex17
Diurnal Cycles Benchmark: a large domain in complex terrain


## Data conversion instructions
A script has been provided to make the data conversion process as easy as possible. The standard way of using the scripts consists of the following steps:
1. Clone or download the repository https://github.com/iat-cener/alex17/archive/master.zip
1. Install the requirements
```bash
pip install -r requirements.txt
```
3. Update the [config/Marinet.yaml](https://github.com/iat-cener/alex17/blob/master/config/Marinet.yaml) file.
3. Edit the [runALEX17.py](https://github.com/iat-cener/alex17/blob/master/runALEX17.py) file, making sure to provide it with functions that can extract your data for a given point and column.
```python
wrf_inp = lib.WrfReader.WrfReader(variables_to_write)
f_get_column = wrf_inp.get_column  # (lat, lon) -> (time, height, variables)
f_get_point = wrf_inp.get_point  # (lat, lon, height) -> (time, variables)
```
The expected output from the functions are labeled, xarray tables. An example of how to define those functions can be found here:
[get_point](https://github.com/iat-cener/alex17/blob/5f1fc540065f1e4b23114e42930fa5f5c7ca4965/lib/WrfReader.py#L322)
and [get_column](https://github.com/iat-cener/alex17/blob/5f1fc540065f1e4b23114e42930fa5f5c7ca4965/lib/WrfReader.py#L332).

(If you don't feel confortable with xarrays, you can try hacking the script and copy the numbers directly to the generated output files [file1](https://github.com/iat-cener/alex17/blob/5f1fc540065f1e4b23114e42930fa5f5c7ca4965/lib/alex17_functions.py#L82), [file2](https://github.com/iat-cener/alex17/blob/5f1fc540065f1e4b23114e42930fa5f5c7ca4965/lib/alex17_functions.py#L130), [file3](https://github.com/iat-cener/alex17/blob/5f1fc540065f1e4b23114e42930fa5f5c7ca4965/lib/alex17_functions.py#L174). This approach is not advised as it will be prone to errors and most likelly it will be more time consuming then understanding the suggested approach).

5. Finally, edit your [simID](https://github.com/iat-cener/alex17/blob/e6b7ec4a7205d97402201430537001fafb154044/runALEX17.py#L24) representing your simulation identifier (should be provided to you).


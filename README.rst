ALEX17 Diurnal Cycles Benchmark: A Large Domain in Complex Terrain
===========================
`Javier Sanz Rodrigo <mailto:jsrodrigo@cener.com>`_, `Pawel Gancarski <mailto:pgancarski@cener.com>`_, `Pedro Alvim De Azevedo Santos <mailto:paas@dtu.dk>`_


Background 
----------
This repository holds model evaluation scripts for the  `ALEX17 Diurnal Cycles benchmark <https://thewindvaneblog.com/alex17-diurnal-cycles-benchmark-a-large-domain-in-complex-terrain-b5029e94485>`_ to test mesoscale-to-microscale flow models in the assessment of wind conditions during several diurnal cycles happening during the `NEWA ALEX17 experiment <https://thewindvaneblog.com/the-alaiz-experiment-alex17-revealing-mountain-valley-large-scale-flow-patterns-6176416dbf2>`_.

The project is integrated in the Phase 3 of the `IEA Task 31 Wakebench <https://community.ieawind.org/task31/home>`_ international framework for wind farm modeling and evaluation.

Scope and Objectives
--------------------
This benchmark is intended for flow models that can reproduce wind conditions at microscale, relevant for wind turbine siting and energy yield assessment, with realistic large-scale forcing characterized by mesoscale modelling. Therefore, the benchmark is designed with transient meso-to-micro modelling in mind. Nevertheless, steady-state modellers are also welcomed to participate and identify situations where the flow is quasi-steady at different levels of thermal stratification.

The ALEX17 intensive campaign, where all the instruments were operational, lasted 4 months. During this period numerous interesting phenomena were captured like gravity waves, hydraulic jumps, etc., but we would like to focus first on testing models in a simpler environment to make sure models can deal with more fundamental challenges, such as:

* How to set up an efficient meshing strategy for a large microscale domain while capturing all relevant topographic features;
* How to deal with surface boundary conditions in the presence of thermal heating/cooling, roughness changes and heterogeneous forest canopies;
* How to parameterize turbulence models;
* How to couple mesoscale and microscale.

Benchmark Guides
----------------
The following blog posts were used to guide benchmark participants:
* `Benchmark guide <https://thewindvaneblog.com/alex17-diurnal-cycles-benchmark-a-large-domain-in-complex-terrain-b5029e94485>`_  
* `ALEX17 background <https://thewindvaneblog.com/the-alaiz-experiment-alex17-revealing-mountain-valley-large-scale-flow-patterns-6176416dbf2>`_  

Data
----
Benchmark input data and simulation data is published open-access in the following data repository: [zenodo dataset]

Citation
--------
You can cite the github repo in the following way:

[zenodo github release]

Installation
------------
We use Jupyter notebooks based on Python 3. We recomend the `Anaconda distribution <https://www.anaconda.com/distribution/>`_ to install python. The libraries used by the notebooks can be installed with 

```bash
pip install -r requirements.txt
```

Benchmark participants can have access to a cloud benchmark environment at CSC's virtual research environment. To this end, please ask for an account to `Pawel Gancarski <mailto:pgancarski@cener.com>`_. Then: 

1. Log in at ``https://notebooks.csc.fi/#/account`` through *Click here for alternate login*. 
2. Click *Join Group* and use the following code: **alex17-xmcp2**
3. Go to the dashboard ´´https://notebooks.csc.fi/#/admin-dashboard´´ and lauch a test notebook. It will take a few minutes to fetch all the data from the b2drop repository to make a local copy. Then it will give you an URL to the Jupyter instance. 

Please note that the instance will expire eventually and will not host your data. Any data or changes you have made during the session have to be downloaded is you want to keep them. 


Data conversion instructions
----------------------------
A script has been provided to make the data conversion process as easy as possible. The standard way of using the scripts consists of the following steps:

1. Update the ``config/Marinet.yaml`` file.
2. Edit the ``runALEX17.py`` file, making sure to provide it with functions that can extract your data for a given point and column.

.. code:: python

	wrf_inp = lib.WrfReader.WrfReader(variables_to_write)
	f_get_column = wrf_inp.get_column  # (lat, lon) -> (time, height, variables)
	f_get_point = wrf_inp.get_point  # (lat, lon, height) -> (time, variables)

The expected output from the functions are labeled, xarray tables. An example of how to define those functions can be found here:
[get_point](https://github.com/iat-cener/alex17/blob/5f1fc540065f1e4b23114e42930fa5f5c7ca4965/lib/WrfReader.py#L322)
and [get_column](https://github.com/iat-cener/alex17/blob/5f1fc540065f1e4b23114e42930fa5f5c7ca4965/lib/WrfReader.py#L332).

(If you don't feel confortable with xarrays, you can try hacking the script and copy the numbers directly to the generated output files [file1](https://github.com/iat-cener/alex17/blob/5f1fc540065f1e4b23114e42930fa5f5c7ca4965/lib/alex17_functions.py#L82), [file2](https://github.com/iat-cener/alex17/blob/5f1fc540065f1e4b23114e42930fa5f5c7ca4965/lib/alex17_functions.py#L130), [file3](https://github.com/iat-cener/alex17/blob/5f1fc540065f1e4b23114e42930fa5f5c7ca4965/lib/alex17_functions.py#L174). This approach is not advised as it will be prone to errors and most likelly it will be more time consuming then understanding the suggested approach).

3. Finally, edit your [simID] representing your simulation identifier (should be provided to you).


License
-------
Copyright 2020 CENER
Licensed under the GNU General Public License v3.0

Acknowledgements
----------------
The authors would like to thank the benchmark participants for their simulations and in-kind support in fine-tuning the benchmark set-up and evaluation methodology. The benchmark is run under the umbrella of IEA-Wind Task 31 with support from the H2020-MARINET2 project, where it is used as a pilot case-study for its virtual research environment. 

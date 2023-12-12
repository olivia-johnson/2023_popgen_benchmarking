# Population genetic simulation: Benchmarking frameworks for non-standard models of natural selection

This repository contains the code to replicate the benchmarking conducted in the aforementioned article. 

The [yaml](benchmarking_env.yml) file can be used to create a conda environment with the python programs required to conduct this benchmarking.

Once cloning the repository, the environment can be created using the following code

```ruby
conda env create -f benchmarking_env.yml

conda activate benchmarking_env
```

To run the simulations, execute the [run_benchmarking.py](run_benchmarking.py) script. This will run all the simulations and generate the data for the burn-in, the forward simulation, and the analysis. 

```ruby
python run_benchmarking.py
```

Figures can be recreated using the [benchmarking_plots.R](benchmarking_plots.R) script. Significance tests and collation of the files created during the simulations are also included in this script.

Any packages not already installed in your R library will be installed by the first line of the script. 



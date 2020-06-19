- Fully manual execution (if e.g. your cluster requires it) is given below

# Semi Manual Installation and Execution
- Using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
package manager all dependencies can be downloaded into a virtual environment, here named `bond_lattice` :
```
conda create --name bond_lattice -f binder/environment.yml
pip install .
```

## Execution
- Modify parameter arrays in `input.json` as described [here](Readme.md)
- Copy files to `run` subfolder
```
mkdir -p run; cp mpi_sweep.py run; cp input.json run;
```
- Run in parallel
```
cd run
mpirun -np ${NPROCS} python mpi_sweep.py
cd ..
````
where `NPROCS` is an integer multiple of the `workers_per_value` parameter

## Analysis
- Run the conversion notebook
```
jupyter nbconvert --to notebook --ExecutePreprocessor.timeout=10000 --execute analyze/convert_to_dataframe.ipynb;
```
- Load and run the free energy calculations cell-by-cell
```
jupyter-notebook analyze/bond_lattice_free_energy.ipynb
```
- If code ran with `JointHist=1` then correlation analysis can also be performed with
```
jupyter-notebook analyze/bond_lattice_correlation.ipynb
```

# Fully Manual Installation and Execution

## Sufficient requirements for execution (i.e. built with `==`):
- `numpy>=1.17`
- `scipy>=1.4.1`
- `mpi4py>=3.0.0`
- `gcc>=7.5`

## Additional requirements for analysis
- `scipy>=1.4.1`
- `jupyter-notebook>=5.7`

This can all  be done with `pip install`

## Compilation
```
cd bond_lattice
make
```
## Execution
- Copy files and source code to `run` subfolder
```
mkdir -p run; cp mpi_sweep.py run; cp input.json run; cp -R bond_lattice run/
```
- Run in parallel
```
cd run
mpirun -np ${NPROCS} python mpi_sweep.py
cd ..
````
where `NPROCS` is an integer multiple of the `workers_per_value` parameter


## Analysis
- Run the conversion notebook
```
jupyter nbconvert --to notebook --ExecutePreprocessor.timeout=10000 --execute analyze/convert_to_dataframe.ipynb;
```
- Load and run the free energy calculations cell-by-cell
```
jupyter-notebook analyze/bond_lattice_free_energy.ipynb
```
- If code ran with `JointHist=1` then correlation analysis can also be performed with
```
jupyter-notebook analyze/bond_lattice_correlation.ipynb
```

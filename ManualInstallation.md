# Manual Installation and Execution

## Sufficient requirements (i.e. built with `==`):
Either install the conda environment:
```
conda create --name bond_lattice -f binder/environment.yml
```
Or install the dependences manually:
- `numpy>=1.17`
- `scipy>=1.4.1`
- `mpi4py>=3.0.0`
- `jupyter-notebook>=5.7`
- `gcc>=7.5`

## Compilation
Inside the new environment you can either use pip to install the `bond_lattice` MD code using:
```
pip install .
```
or you can compile it on your own:
```
cd bond_lattice
make
```

## Execution
- Modify parameter arrays in `input.json` as described [here](Readme.md)
- Copy files to `run` subfolder
```
mkdir -p run; cp mpi_sweep.py run; cp input.json run;
```
- If you did *not* install the `bond_lattice` MD code using pip copy the source code also
```
cp -R bond_lattice run/
```
- Run in parallel
```
cd run
mpirun -np ${NPROCS} python mpi_sweep.py
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

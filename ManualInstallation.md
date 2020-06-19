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
- Modify parameter arrays in `input.json`
- Run in parallel with `mpirun -np ${NPROCS} python mpi_sweep.py`. Reminder: If you installed the `bond_lattice` MD code using pip you have to either remove/rename the bond_lattice module in the source directory or copy the `mpi_sweep.py` file and the `input.json` file to a separate directory. Otherwise python tries to load the local module first.

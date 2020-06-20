rule md:
    input:
        "input.json"
    output:
        "analyze/bond_lattice_free_energy.nbconvert.ipynb",
        "analyze/bond_lattice_correlation.nbconvert.ipynb",
    conda:
        "binder/environment.yml"
    shell:
        # Execution of the MD code 
        # "pip install .;"
        # "mkdir -p run; cp mpi_sweep.py run; cp input.json run; cd run;"
        # "mpirun --oversubscribe -np 10 python mpi_sweep.py;"
        # "cd ..;"
        # 
        # Convert format
        "jupyter nbconvert --to notebook"
        "    --ExecutePreprocessor.timeout=10000"
        "    --execute analyze/convert_to_dataframe.ipynb;"
        # Bond based approximation of the free energy 
        "jupyter nbconvert --to notebook"
        "    --ExecutePreprocessor.timeout=10000"
        "    --execute analyze/bond_lattice_free_energy.ipynb;"
        # Bond correlation 
        "jupyter nbconvert --to notebook"
        "    --ExecutePreprocessor.timeout=10000"
        "    --execute analyze/bond_lattice_correlation.ipynb"

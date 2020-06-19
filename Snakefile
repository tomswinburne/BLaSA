rule md:
    input:
        "input.json"
    output:
        "analyze/bond_lattice_free_energy.nbconvert.ipynb",
        "analyze/dEnergy_dFreeEnergy.pdf"
    conda:
        "binder/environment.yml"
    shell:
        "pip install .;"
        "mkdir -p run; cp mpi_sweep.py run; cp input.json run; cd run;"
        "mpirun --oversubscribe -np 10 python mpi_sweep.py;"
        "cd ..;"
        "jupyter nbconvert --to notebook"
        "    --ExecutePreprocessor.timeout=10000"
        "    --execute analyze/convert_to_dataframe.ipynb;"
        "jupyter nbconvert --to notebook"
        "    --ExecutePreprocessor.timeout=10000"
        "    --execute analyze/bond_lattice_free_energy.ipynb"

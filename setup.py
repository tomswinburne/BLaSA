from distutils.core import setup, Extension

bond_lattice = Extension('bond_lattice.libmcsim', sources = ['bond_lattice/chain_lattice.cpp'])

setup (name = 'bond_lattice',
       version = '0.0.1',
       packages=["bond_lattice"],
       ext_modules = [bond_lattice])
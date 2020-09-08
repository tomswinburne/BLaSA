from ctypes import *
import os
import glob
import numpy as np
from numpy.ctypeslib import ndpointer
from scipy.optimize import minimize
from scipy.special import iv,erf

class bond_lattice:
    def __init__(self,\
        a0=1.0, strain=0.0, N=6,seed=40,rank=0,fcc=True,\
        AL=1.0, D=1.0, RT=0.0, CentralForce=False,\
        min_r=-0.5,max_r=4.5,bins=501,CorrelationType=0,CorrelationTarget=0.5,\
        steps=1000,therm_steps=500,libname="libmcsim"):

        """relative path to shared library"""
        self.lib_path = glob.glob(os.path.join(os.path.dirname(__file__), libname + "*.*so"))[-1]

        self.mclib = cdll.LoadLibrary(self.lib_path)
        self.D = D
        self.AL = AL
        self.RT = RT
        self.a0 = a0
        self.N = N
        self.seed = seed
        self.steps = steps
        self.therm_steps = therm_steps
        self.fcc = fcc
        self.CorrelationType = CorrelationType
        self.CorrelationTarget = CorrelationTarget #
        self.rank = rank
        self.max_r = max_r
        self.min_r = min_r
        self.bins = bins
        self.CentralForce = CentralForce

    def run(self,strain=0.0,T=0.1):
        H = np.zeros(4*self.bins,np.int32)
        if self.CorrelationType==1:
            CH = np.zeros(self.bins*self.bins,np.int32)
        elif self.CorrelationType==2:
            CH = np.zeros(3*self.bins*self.bins,np.int32)
        else:
            CH = np.zeros(1,np.int32)
        V = np.zeros(self.bins+6)
        r = np.zeros(self.bins)
        if self.rank==0:
            print(CH.shape,r.shape)

        """
        bond_lattice_sim(
        double a0, double epsilon, double T, int N,
        double AL, double D, double RT, bool CentralForce,
        double min_r, double max_r, int bins, int CorrelationType, double CorrelationTarget,
        int *H, double *V, int *CH, double *r,
        unsigned seed, int tsample, int ttherm, int rank, bool fcc)
        """

        self.mclib.bond_lattice_sim.argtypes = \
                    [c_double, c_double, c_double, c_int,\
                    c_double, c_double, c_double, c_bool,
                    c_double, c_double, c_int, c_int, c_double,\
                    ndpointer(c_int32, flags="C_CONTIGUOUS"),\
                    ndpointer(c_double, flags="C_CONTIGUOUS"),\
                    ndpointer(c_int32, flags="C_CONTIGUOUS"),\
                    ndpointer(c_double, flags="C_CONTIGUOUS"),\
                    c_uint,c_int,c_int,c_int,c_bool]

        self.mclib.bond_lattice_sim.restypes = None
        self.mclib.bond_lattice_sim(\
                self.a0,strain*0.01,T,self.N,\
                self.AL,self.D,self.RT,self.CentralForce,\
                  self.min_r,self.max_r,self.bins,self.CorrelationType,self.CorrelationTarget,\
                  H,V,CH,r,\
                  self.seed,self.steps,self.therm_steps,self.rank,self.fcc,)

        if self.CorrelationType>0:
            return r,H,V,CH
        else:
            return r,H,V

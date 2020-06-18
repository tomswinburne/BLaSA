from ctypes import *
import os
import glob
import numpy as np
from numpy.ctypeslib import ndpointer
from scipy.optimize import minimize
from scipy.special import iv,erf

class bond_lattice:
    def __init__(self,D=1.0,AL=1.0,RT=0.0,a0=1.0,N=6,\
        min_r=-0.5,max_r=4.5,bins=501,am=1.0,steps=1000,therm_steps=500,\
        seed=40,rank=0,fcc=True,CorrelationType=0,CentralForce=False,libname="libmcsim"):
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
        self.CorrelationType = CorrelationType #
        self.rank = rank
        self.max_r = max_r
        self.min_r = min_r
        self.bins = bins
        self.am = am
        self.CentralForce = CentralForce

    def morse(self,r):
        return self.D*(1.0+np.exp(-2.0*self.AL*(r-self.a0))-2.0*np.exp(-self.AL*(r-self.a0)))
    def dmorse(self,r):
        return -2.0*self.D*self.AL*(np.exp(-2.0*self.AL*(r-self.a0))-np.exp(-self.AL*(r-self.a0)))

    def find_linear_correction(self,am=1.0,T=0.1,quad=1.0,nstep=2000,maxr=2.0):
        # base on quad value
        r = np.linspace(-0.5,maxr,nstep,endpoint=True)
        brho = np.exp(-self.morse(r)/T - 0.5*quad*(r-am*self.a0)**2/T)
        brho /= brho.sum()
        def f(m):
            rm = (r-am*self.a0).mean()
            rho = brho * np.exp(-m*(r-am*self.a0)/T + m*rm/T)
            lhs = am * self.a0* rho.sum()
            rhs = (rho*r).sum()
            res = (np.log(lhs/rhs)**2).sum()
            return res
        solver = minimize(f, [0.1], method='Nelder-Mead', tol=1e-20)
        return solver.x

    def run(self,am=1.0,T=0.1):
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
            print(CH.shape)
        """
        void simpos3D(double a0, double am,
        int N,
        double AL, double D, double RT, double T, double min_r, double max_r,
        int bins, int *H, double *V, int *CH,
        double *r,
        unsigned seed, int SAMPLE, int THERM, int rank,
        bool fcc,
        int CorrelationType, bool CentralForce)
        """
        self.mclib.simpos3D.argtypes = [c_double, c_double, c_int,\
                    c_double,c_double,c_double,c_double,c_double,c_double,\
                    c_int,\
        ndpointer(c_int32, flags="C_CONTIGUOUS"),\
        ndpointer(c_double, flags="C_CONTIGUOUS"),\
        ndpointer(c_int32, flags="C_CONTIGUOUS"),\
        ndpointer(c_double, flags="C_CONTIGUOUS"),\
        c_uint,c_int,c_int,c_int,c_bool,c_int32,c_bool]
        print(self.fcc,self.CorrelationType,self.CentralForce)
        self.mclib.simpos3D.restypes = None
        self.mclib.simpos3D(self.a0,am,self.N,\
                  self.AL,self.D,self.RT,T,\
                  self.min_r,self.max_r,\
                  self.bins,\
                  H,V,CH,r,\
                  self.seed,self.steps,\
                  self.therm_steps,self.rank,\
                  self.fcc,self.CorrelationType,self.CentralForce)

        if self.CorrelationType>0:
            return r,H,V,CH
        else:
            return r,H,V

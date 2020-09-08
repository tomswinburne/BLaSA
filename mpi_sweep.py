import os,json
import numpy as np
from mpi4py import MPI
from bond_lattice import bond_lattice

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# Read input file
if rank == 0:
    with open("input.json", "r") as f:
        input_data = json.load(f)
else:
    input_data = None
input_data = comm.bcast(input_data, root=0)

procs_per_worker = input_data['workers_per_value']
if size % procs_per_worker!=0:
    comm.Barrier()
    print("NONFACTORIZING PROCS OR WORKERS! EXITING!")
    MPI.Finalize()
    exit()

bins = input_data['md']['bins']
N = input_data['md']['n']
THERM = input_data['md']['therm']
STEPS = input_data['md']['steps']
D0 = input_data['potential']['D0']
a0 = input_data['potential']['a0']
AL = input_data['potential']['AL']
rlim = [input_data['potential']['r_min'], input_data['potential']['r_max']]
strain_array = input_data['strain_array']
T_array = input_data['t_array']
RT_array = input_data['RT_array']

seed = int(1000.0*np.random.uniform())


# Histgram properties
dump_folder = input_data['dump_folder']
potname = "M"
if rank==0:
    os.system("mkdir -p %s" % dump_folder)



kb = 8.617e-5

# job splitting. Will iterate over all T,strain,AL combinations

E_array=[]
strain_RT_T_array = []
for strain in strain_array:
    for RT in RT_array:
        for T in T_array:
            strain_RT_T_array.append([strain,RT,T])

color = rank // procs_per_worker

lrank = rank % procs_per_worker

nworkers = size // procs_per_worker

njobs = len(strain_RT_T_array)
jobs_per_worker = njobs // nworkers
c=0
jl=np.arange(nworkers+1) * jobs_per_worker
jl[-1] = njobs

if lrank==0:
    for i in np.arange(jl[color],jl[color+1]):
        strain,RT,T = strain_RT_T_array[i]

lcomm = comm.Split(color,lrank)

comm.Barrier()


"""
CorrelationType = 0 : Just |b_perp|,b_para,|b| marginals
CorrelationType = 1 : += |b_perp|,b_para joint distribution
CorrelationType = 2 : += b_para joint dist. for (l,l+7) and (l,l+1)
CorrelationType > 0 can give quite large files, GB of data in a batch
"""

CorrelationType = 2 * input_data["JointHist"]

for ji in np.arange(jl[color],jl[color+1]):
    strain,RT,T = strain_RT_T_array[ji]
    sim_pbc = bond_lattice(a0=a0,N=N,seed=seed+rank,rank=rank,fcc=True,\
                        AL=1.0, D=1.0, RT=0.0, CentralForce=False,\
                        min_r=-0.5,max_r=4.5,bins=bins,\
                        CorrelationType=0,CorrelationTarget=0.5,\
                        steps=1000,therm_steps=500,libname="libmcsim")

    dump_string = "N%d_S%dk_RT%f_T%f_a%f_%s" % (N,STEPS//1000,RT,T,strain,potname)

    # collect data from each worker
    if CorrelationType>0:
        r,H,UmU,CH = sim_pbc.run(T=T*kb,strain=strain)
    else:
        r,H,UmU = sim_pbc.run(T=T*kb,strain=strain)

    # Take value and square value for ensemble averages
    Ums = np.zeros(12)
    Ums[:6] = UmU[-6:]
    Ums[-6:] = UmU[-6:]**2

    BondPotential  =  UmU[:-6]

    lcomm.Barrier()

    data = lcomm.gather(H, root=0)
    Udata = lcomm.gather(Ums, root=0)
    if CorrelationType>0:
        cdata = lcomm.gather(CH, root=0)
    if lrank == 0:
        H = data[0]
    
    meanstdU = Udata[0] / procs_per_worker

    for i in range(1,procs_per_worker):
        meanstdU += Udata[i] / procs_per_worker
        H += data[i]
    meanstdU[-6:] = np.sqrt(meanstdU[-6:]-meanstdU[:6]**2)

    equ = 1.5*kb

    print("pre T:\n",T,"Teq:",meanstdU[:2].sum()/equ,"Tvir",meanstdU[2]/kb,end="| ")
    print("Teq:",meanstdU[6:8].sum()/equ,"Tvir",meanstdU[8]/kb)

    print("post T:\n",T,"Teq:",meanstdU[3:5].sum()/equ,"Tvir",meanstdU[5]/kb)
    print("Teq:",meanstdU[9:11].sum()/equ,"Tvir",meanstdU[11]/kb)

    for ss in [0,3,6,9]:
        meanstdU[ss:ss+2] /= equ
        meanstdU[ss+2] /= kb

    ead = np.zeros(15)

    ead[0],ead[1],ead[2] = T,strain,RT
    ead[3:] = meanstdU

    print(r.shape,H.shape,bins)
    E_array.append(ead)
    np.savetxt(os.path.join(dump_folder,'r_H_Hp_Hd_U_%s' % dump_string),
    np.vstack((r,H[:bins],H[bins:2*bins],H[2*bins:3*bins],H[3*bins:],
    BondPotential)).T,fmt="%f %d %d %d %d %f")

    if CorrelationType>0:
        CH = cdata[0]
        for i in range(1,procs_per_worker):
            CH += cdata[i]
        np.savetxt(os.path.join(dump_folder,'H_TL_%s' % dump_string),CH,fmt='%d')
g_E_array = comm.gather(E_array,root=0)
comm.Barrier()
#if lrank==0:
#	local_f.close()
if rank==0:
    f_E_array=[]
    for E_array in g_E_array:
        if len(E_array)>0:
            for s_E_array in E_array:
                if len(s_E_array)>0:
                    f_E_array.append(s_E_array)
                    print(s_E_array)
    np.savetxt(os.path.join(dump_folder,"E_data_2"),np.r_[f_E_array])
    f = open(os.path.join(dump_folder,"E_data"),'ba')
    np.savetxt(f,np.r_[f_E_array])
    f.close()

MPI.Finalize()

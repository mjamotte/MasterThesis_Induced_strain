import numpy as np
import scipy as sc
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from numpy import linalg
import cmath
import time
from mpl_toolkits import mplot3d
import sys
sys.path.append('/home/maxime/Desktop/Codes/')
import lib_GrossPitaevskii as GP
import lib_Manage_files as Mf
import lib_Honeycomb_Lattice as FHoney

############## Trigonal strain ##################

tau_list = [0.02]
Nx_s = [101]
Ny_s = [100]

gamma2 = 1
gamma3 = 1

BoundCond = ['Trigonal']
BC = BoundCond[0]

t = 1
a1 = np.array([1,0])
a2 = np.array([0,1])
kxa = 0
kya = 0

start = time.time()
for Nx in Nx_s:
	for Ny in Ny_s:

		sites_dic = FHoney.square_lattice(Nx,Ny,a1,a2)

		for tau in tau_list:

			system = [BC,Nx,Ny,t,tau,gamma2,gamma3]

			if tau<=4/(Nx-1) and tau<=2/(Ny-1) or kxa==0: # if no Lifshitz transition
			
				H_trig = FHoney.H_honey(system,kxa,kya,sites_dic)
				eigVal= linalg.eigh(H_trig)[0]
				#H_trig = sc.sparse.coo_matrix(FHoney.H_honey(system,kxa,kya,sites_dic))
				#print(H_trig)
				#eigVal = sc.sparse.linalg.eigsh(H_trig,k=600,which='SM')[0]
				#print(eigVal)

				# Save in a file
				data = [system,eigVal] #eigVectors]
				dataID = '%s_Nx_%.1iNy_%.1i_tau_%.3f' % (system[0],Nx,Ny,system[4])
				np.save('Datas/'+dataID+'.npy',data)

print("Execution time:",time.time()-start,"secondes")

#pyplot.pcolormesh(abs(H_trig)) # doesn't work for sparse matrix
#pyplot.colorbar()
#pyplot.show()


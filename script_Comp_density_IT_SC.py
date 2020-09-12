import numpy as np
import matplotlib.colors as color
import matplotlib.axes as axes
import matplotlib.pyplot as pyplot
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib import rc
from numpy import linalg
import cmath
import time
from mpl_toolkits import mplot3d
from random import randint

import os
import sys
path_above_Codes = (sys.path[0]).split('Codes')[0]
sys.path.append(path_above_Codes+'Codes/')
## Import the libraries in the parent folder "Codes"
import lib_GrossPitaevskii as GP
import lib_Manage_files as Mf
import lib_Optical_Honeycomb_Lattice as OHL

rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
rc('text', usetex=True)

#### COMPARISON OF DENSITIES FOR IT AND SC ######

## Initialize some variables

mu_IT_s = np.array([])
mu_SC_s = np.array([])
n0_IT_s = []
n0_SC_s = []
U_IT_s = np.array([])
U_SC_s = np.array([])
n_TF_IT_s = []
n_TF_SC_s = []

Efun_SC = np.array([]) # to compute mu = dE_SC/dN
Efun_IT = np.array([]) # to compute mu = dE_IT/dN

## Parameters to pick the datas from Datas_SC nd Datas_IT

list_args = ['N','V','U']
split = ['N_','V_','U_']
params_interv = {'U': [1e-5,1.3e-4],'V': [3.5e-5],'N': [1e5]}
extension = '.npy'

## Import datas from Datas_IT

path = '/home/maxime/Desktop/Codes/Optical_1D_Lattice_BEC/Imaginary_time/Datas_IT/'
files_names = Mf.select_files(path,extension,list_args,params_interv,split)
split = ['U_']
list_args = ['U']
files_names = Mf.reorder_filesnames(files_names,list_args,split)

for filename in files_names[::5]:

	data = np.load(filename+'.npy',allow_pickle=True)
	args_syst = data[0]
	J = args_syst['J']
	Nx = args_syst['Nx']
	V = args_syst['V']
	U_IT = args_syst['U']
	N = args_syst['N']

	args_init = data[1]
	#mu_all = data['mu_all']	
	mu_IT = data[2]
	psi0_IT = data[3]
	E_funct_IT = data[4]
	n0_IT = np.abs(psi0_IT)**2*N

	## Analytical gaussian 

	Trap = GP.trap_1D(args_syst)
	m = 1/(2*J)
	w0 = np.sqrt(V/m)
	x0 = (Nx-1)/2
	positions = np.arange(Nx)

	## Stock datas

	Efun_IT = np.append(Efun_IT,E_funct_IT)
	U_IT_s = np.append(U_IT_s,U_IT)
	mu_IT_s = np.append(mu_IT_s,mu_IT)
	n0_IT_s.append(n0_IT)

	## Thomas-Fermi IT

	n_TF_IT = (mu_IT-Trap)/U_IT/N# divide by because here, normalised to 1
	n_TF_IT_s.append(n_TF_IT-np.max(n_TF_IT)+np.max(n0_IT))

n0_IT_s = np.array(n0_IT_s)
n_TF_IT_s = np.array(n_TF_IT_s)

## Import datas for SC

path = '/home/maxime/Desktop/Codes/Optical_1D_Lattice_BEC/Self_consistent/Datas_SC/'
list_args = ['N','V','U']
split = ['N_','V_','U_']
extension = '.npy'
files_names = Mf.select_files(path,extension,list_args,params_interv,split)
split = ['U_']
list_args = ['U']
files_names = Mf.reorder_filesnames(files_names,list_args,split)

for filename in files_names[::5]:

	data = np.load(filename+'.npy',allow_pickle=True)
	args_syst = data[0]
	J = args_syst['J']
	Nx = args_syst['Nx']
	V = args_syst['V']
	U_SC = args_syst['U']
	N = args_syst['N']

	args_init = data[1]
	#mu_all = data['mu_all']	
	mu_SC = data[2]
	psi0_SC = data[3]
	E_funct_SC = data[4]
	n0_SC = np.abs(psi0_SC)**2*N

	Trap = GP.trap_1D(args_syst)

	Efun_SC = np.append(Efun_SC,E_funct_SC)
	U_SC_s = np.append(U_SC_s,U_SC)
	mu_SC_s = np.append(mu_SC_s,mu_SC)
	n0_SC_s.append(n0_SC)

	## Thomas-Fermi SC

	n_TF_SC = (mu_SC-Trap)/U_SC/N# divide by because here, normalised to 1
	n_TF_SC_s.append(n_TF_SC-np.max(n_TF_SC)+np.max(n0_SC))

n0_SC_s = np.array(n0_SC_s)
n_TF_SC_s = np.array(n_TF_SC_s)

## Gaussian (via Hermite functions)

V0 = params_interv['V'][0]
m = 1/(2*J)
w0 = np.sqrt(V0/m)
x0 = (Nx-1)/2
positions = np.arange(Nx)
Hermite_0 = (m*w0/np.pi)**0.25*np.exp(-m*w0*(positions-x0)**2/2)

## Plots

fig_OL1D = pyplot.figure(figsize=(8,8))
ax1 = pyplot.subplot(111)
ax1.set_xlabel('$x/a$',fontsize=20)


for i in range(len(n0_IT_s)):
	ax1.plot(positions,n0_IT_s[i],'-',label="$n_{IT}$, $U = %.2e$"%(U_IT_s[i]))
	ax1.plot(positions,n0_SC_s[i],'o',fillstyle='none',label="$n_{SC}$, $U = %.2e$"%(U_SC_s[i]))
	#ax1.plot(positions,n_TF_IT_s[i],'--',fillstyle='none',label="$IT(TF), U = {}$".format(U_IT_s[i]))
	#ax1.plot(positions,n_TF_SC_s[i],'--',fillstyle='none',label="$SC(TF), U = {}$".format(U_SC_s[i]))
	#ax1.plot(positions,Hermite_0**2,'-')

# pyplot.suptitle('Comparison of densities IT and SC for %s, %s,\
# 	 Nx = %.1i, J = %.2f, V = %.3e, U = %.3e' % \
# 	 (args_syst['Trap'],args_syst['Symm'],\
# 	args_syst['Nx'],args_syst['J'],args_syst['V'],args_syst['U']))

ax1.legend(loc=1,fontsize=18);
ax1.grid(axis='both')
ax1.set_ylim(0)
ax1.set_xlim(x0)
pyplot.xticks(fontsize=16)
pyplot.yticks(fontsize=16)
pyplot.show()

## Grid, legend and save the figure

temp = '1D_comp_densities_ITSC_%s_%s_Nx_%.1i_J_%.2f_V_%.3e_U_%.3e' %\
		(args_syst['Trap'],args_syst['Symm'],\
		 args_syst['Nx'],args_syst['J'],args_syst['V'],args_syst['U'])
#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
#fig_OL1D.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
#path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
#fig_OL1D.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/fig'+temp+".pdf") # Linux
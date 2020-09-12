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

#### COMPARISON OF ENERGIES FOR IT AND SC ######

## Figure

fig_OL1D = pyplot.figure(figsize=(8,8))
ax1 = pyplot.subplot(111)
ax1.set_xlabel('U')

## Initialize some variables

U_IT_s = np.array([])
U_SC_s = np.array([])
mu_IT_s = np.array([])
mu_SC_s = np.array([])
E_kin_s_IT = np.array([])
E_tr_s_IT = np.array([])
E_U_s_IT = np.array([])
E_kin_s_SC = np.array([])
E_tr_s_SC = np.array([])
E_U_s_SC = np.array([])

Efun_SC = np.array([]) # to compute mu = dE_SC/dN
Efun_IT = np.array([]) # to compute mu = dE_IT/dN

## Parameters to pick the datas from Datas_SC nd Datas_IT

list_args = ['N','V','U']
split = ['N_','V_','U_']
params_interv = {'U': [2e-6,2e-4],'V': [3.5e-5],'N': [1e5-1,1e5,1e5+1]}
extension = '.npy'

## Import datas from Datas_IT

path = path_above_Codes+'Codes/Optical_1D_Lattice_BEC/Imaginary_time/Datas_IT/'
files_names = Mf.select_files(path,extension,list_args,params_interv,split)
split = ['U_']
list_args = ['U']
files_names = Mf.reorder_filesnames(files_names,list_args,split)

for m in np.arange(0,len(files_names)+1,3):
	split = ['N_']
	list_args = ['N']
	files_names[m:m+3] = Mf.reorder_filesnames(files_names[m:m+3],list_args,split)

for filename in files_names:

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
	n0_IT = np.abs(psi0_IT)**2 

	Trap = GP.trap_1D(args_syst)

	Efun_IT = np.append(Efun_IT,E_funct_IT)
	U_IT_s = np.append(U_IT_s,U_IT)
	mu_IT_s = np.append(mu_IT_s,mu_IT)
	E_kin_s_IT = np.append(E_kin_s_IT,GP.E_kin(psi0_IT,args_syst))
	E_tr_s_IT = np.append(E_tr_s_IT,GP.E_trap(psi0_IT,args_syst))
	E_U_s_IT = np.append(E_U_s_IT,GP.E_int(psi0_IT,args_syst))


## Import datas for SC

path = path_above_Codes+'Codes/Optical_1D_Lattice_BEC/Self_consistent/Datas_SC/'
list_args = ['N','V','U']
split = ['N_','V_','U_']
extension = '.npy'
files_names = Mf.select_files(path,extension,list_args,params_interv,split)
split = ['U_']
list_args = ['U']
files_names = Mf.reorder_filesnames(files_names,list_args,split)

for m in np.arange(0,len(files_names)+1,3):
	split = ['N_']
	list_args = ['N']
	files_names[m:m+3] = Mf.reorder_filesnames(files_names[m:m+3],list_args,split)


for filename in files_names:
	
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

	## Analytical gaussian 

	Trap = GP.trap_1D(args_syst)
	m = 1/(2*J)
	w0 = np.sqrt(V/m)
	x0 = (Nx-1)/2
	positions = np.arange(Nx)

	Efun_SC = np.append(Efun_SC,E_funct_SC)
	U_SC_s = np.append(U_SC_s,U_SC)
	mu_SC_s = np.append(mu_SC_s,mu_SC)
	E_kin_s_SC = np.append(E_kin_s_SC,GP.E_kin(psi0_SC,args_syst))
	E_tr_s_SC = np.append(E_tr_s_SC,GP.E_trap(psi0_SC,args_syst))
	E_U_s_SC = np.append(E_U_s_SC,GP.E_int(psi0_SC,args_syst))

## Compute mu = dE/dN

dN = 1
dEfun_SCdN = GP.dEdN_O2(Efun_SC,dN)
dEfun_ITdN = GP.dEdN_O2(Efun_IT,dN)	

## Plots

ax1.semilogx(U_IT_s,E_kin_s_IT/N+2*J,'bs',fillstyle='none',label="$E^{IT}_{kin}[\psi]/N$")
ax1.semilogx(U_IT_s,E_tr_s_IT/N,'gs',fillstyle='none',label="$E^{IT}_{trap}[\psi]/N$")
ax1.semilogx(U_IT_s,E_U_s_IT/N,'rs',fillstyle='none',label="$E^{IT}_U[\psi]/N$")
ax1.semilogx(U_SC_s,E_kin_s_SC/N+2*J,'b^',fillstyle='none',label="$E^{SC}_{kin}[\psi]/N$")
ax1.semilogx(U_SC_s,E_tr_s_SC/N,'g^',fillstyle='none',label="$E^{SC}_{trap}[\psi]/N$")
ax1.semilogx(U_SC_s,E_U_s_SC/N,'r^',fillstyle='none',label="$E^{SC}_U[\psi]/N$")
ax1.semilogx(U_SC_s,E_tr_s_SC*2/N-2*(E_kin_s_SC/N+2*J),\
				'k-',fillstyle='none',label="Vir. theor")

pyplot.suptitle('Comparison of energies IT and SC for %s, %s,\
	 Nx = %.1i, J = %.2f, V = %.3e' % \
	 (args_syst['Trap'],args_syst['Symm'],\
	args_syst['Nx'],args_syst['J'],args_syst['V']))

ax1.legend(loc=2);
ax1.grid(axis='both')
#pyplot.show()

temp = '1D_comp_energies_ITSC_%s_%s_Nx_%.1i_J_%.2f_V_%.3e' %\
		(args_syst['Trap'],args_syst['Symm'],\
		 args_syst['Nx'],args_syst['J'],args_syst['V'])

fig_OL1D.savefig("Figures_OL_1D_BEC_Comparison/fig_"+temp+".pdf")

#################################

fig_OL1D2 = pyplot.figure(figsize=(8,8))
ax1 = pyplot.subplot(111)
ax1.set_xlabel('$U$',fontsize='20')

#ax1.semilogx(U_IT_s,Efun_IT/N+2*J,'ks',fillstyle='none',label="$E^{IT}[\psi]/N$")
ax1.loglog(U_IT_s,mu_IT_s+2*J,'cs',fillstyle='none',label="$\mu_{IT}+2J$",markersize=12)
#ax1.semilogx(U_SC_s,Efun_SC/N+2*J,'k^',fillstyle='none',label="$E^{SC}[\psi]/N$")
ax1.loglog(U_SC_s,mu_SC_s+2*J,'g^',fillstyle='none',label="$\mu_{SC}+2J$",markersize=12)

ax1.loglog(U_SC_s[::3],dEfun_SCdN+2*J,'c+',fillstyle='none',label="$\partial E_{SC}/\partial N$",markersize=12)
ax1.loglog(U_IT_s[::3],dEfun_ITdN+2*J,'g*',fillstyle='none',label="$\partial E_{IT}/\partial N$")
ax1.loglog(U_SC_s,(U_SC_s*np.sqrt(V/2)*3/4*N)**(2/3),'-m',label="$\\propto (UN)^{2/3}$")

# pyplot.suptitle('Compar. $E[\psi]$ and $\mu$ with IT and SC for %s, %s,\
# 	 Nx = %.1i, J = %.2f, V = %.3e' % \
# 	 (args_syst['Trap'],args_syst['Symm'],\
# 	args_syst['Nx'],args_syst['J'],args_syst['V']))

ax1.legend(loc=2,fontsize=16);
pyplot.xticks(fontsize=16)
pyplot.yticks(fontsize=16)
ax1.grid(axis='both')
pyplot.show()

temp = '1D_comp_funct_mu_ITSC_%s_%s_Nx_%.1i_J_%.2f_V_%.3e' %\
		(args_syst['Trap'],args_syst['Symm'],\
		 args_syst['Nx'],args_syst['J'],args_syst['V'])

#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
#fig_OL1D2.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
fig_OL1D2.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/fig'+temp+".pdf") # Linux





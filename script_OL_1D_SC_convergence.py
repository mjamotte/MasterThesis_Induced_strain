import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from numpy import linalg
import cmath
import time
from matplotlib import rc
from scipy.integrate import ode

from numpy import linalg
import time

import os
import sys
path_above_Codes = (sys.path[0]).split('Codes')[0]
sys.path.append(path_above_Codes+'Codes/')
## Import the libraries in the parent folder "Codes"
import lib_GrossPitaevskii as GrPi
import lib_Manage_files as Mf
import lib_Optical_Honeycomb_Lattice as OHL

rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
rc('text', usetex=True)


# mu = -2
# t = 0
# U = 0.01

# Y, X = numpy.mgrid[-0.1:0.1:200j, -20:0:200j]
# U = Y
# V = -mu*X + t**2*X + U*X**3

# pyplot.streamplot(X,Y,U,V,density=1)
# pyplot.show()

def conv_mu_SC(args_syst,args_init):

	N = args_syst['N']
	U = args_syst['U']
	H_KV = args_init['H_KV']
	err = args_init['err_SC'] # condition to stop convergence
	mu_all = []
	psi0_all = [np.array([])]
	err_psi = []

	
	#H_KV = sc.sparse.coo_matrix(H_KV) # KV = Kinetic + Trap	

	if 'psi_init' in args_init:
		psi_old = args_init['psi_init']

	else:
		#E0_old,eigVecs = sc.sparse.linalg.eigsh(H_KV,which='SA',k=1)
		E0_old,eigVecs = linalg.eigh(H_KV)
		psi_old = np.matrix.transpose(eigVecs[:,0])

	counterSC = 0
	lam = 1e-3 #1/(100*N+1) # to avoid oscillations in energy that slow down the code
	#flag = 0

	while True:

		#H_U = sc.sparse.coo_matrix(H_int(psi_old,args_syst))
		#E0_new,eigVecs = sc.sparse.linalg.eigsh(H_KV+H_U,which='SA',k=1)
		H_U = GrPi.H_int(psi_old,args_syst)
		eigVal,eigVecs = linalg.eigh(H_KV+H_U)
		psi_new = np.matrix.transpose(eigVecs[:,0])

		psi_lam = np.sqrt((1-lam)*psi_old**2 + lam*psi_new**2)
		mu_all.append(eigVal[0])

		# err_psi = np.abs(np.max(psi_lam)\
		# 			-np.max(psi_old))\
		# 			/np.abs(np.max(psi_lam))
		err_psi.append(np.abs(np.max(np.abs(psi_lam)**2)\
					-np.max(np.abs(psi_old)**2))\
					/np.max(np.abs(psi_lam)**2))

		if err_psi[-1]<err:
			break

		psi_old = psi_lam
		counterSC += 1		

	print('Number of iterations of self-consistent method =',counterSC)

	return np.array(mu_all),counterSC,psi_lam,err_psi


###############################################
def conv_IT(args_syst,args_init):

	## Initialisation 
	H_KV = args_init['H_KV']
	#H_KV = sc.sparse.coo_matrix(H_KV) # KV = Kinetic + Trap

	if 'psi_init' in args_init:
		psi_old = args_init['psi_init'] 
		psi_old = np.concatenate((np.real(psi_old), np.imag(psi_old)))

	else:
		#gauss = sc.sparse.linalg.eigsh(H_KV)[1][:,0]
		gauss = linalg.eigh(H_KV)[1][:,0]
		psi_old = np.concatenate((np.real(gauss), np.imag(gauss)))

	## parameters for set_integrator and GP
	tol = 1e-9 # tolerance
	nsteps = np.iinfo(np.int32).max
	solver = ode(GrPi.GP)
	solver.set_f_params(args_syst,args_init) # parameters needed in GP_t_real
	solver.set_integrator('dop853', atol=tol, rtol=tol, nsteps=nsteps)

	## Evolution
	t = 0
	dt = args_init['dt']
	err = args_init['err_IT']
	counterIT = 0
	dim = len(psi_old)
	err_psi = []

	while True:
		## time evolution
		solver.set_initial_value(psi_old, t)
		solver.integrate(t+dt)
		t = t+dt
		psi_new = solver.y

		psi_old_co = psi_old[:int(dim/2)] + 1j*psi_old[int(dim/2):]
		psi_new_co = psi_new[:int(dim/2)] + 1j*psi_new[int(dim/2):] # not renorm. yet

		## renormalize
		psi_new = psi_new/np.sqrt(np.sum(abs(psi_new)**2))
		psi_new_co = psi_new[:int(dim/2)] + 1j*psi_new[int(dim/2):]

		err_psi.append(np.sqrt(np.abs(np.max(np.abs(psi_new_co)**2\
					-np.max(np.abs(psi_old_co)**2))))\
					/np.sqrt(np.max(np.abs(psi_new_co)**2)))

		if err_psi[-1]<err:
			break

		psi_old = psi_new
		counterIT += 1

	if solver.successful():
		sol = solver.y
		sol = sol/np.sqrt(sum(abs(sol)**2))
		sol_re = sol[:int(dim/2)]
		sol_im = sol[int(dim/2):]
		psi0 = sol_re + 1j*sol_im

	print('The number of iterations for IT =', counterIT)

	return psi0,err_psi,counterIT
###############################################

fig_OL1D = pyplot.figure(figsize=(8,8))
k = 0
for U in [1e-5]:
	args_syst = {
				'syst' : '1D',
				'J' : 1,
				'N' : 1e4,
				'V' : 0.5e-4,
				'Nx' : 201,
				'U' : U,
				'Method' : 'SC',
				'Trap' : 'Harmonic',
				'Symm' : 'Isotropic',
				}
	args_syst.update({'sites_dic' : GrPi.lattice_1D(args_syst)})
	H_KV = GrPi.H_1D(args_syst)

	if k==0:
	 	args_init = {
	 	'H_KV' : H_KV,
	 	'err_SC' : 1e-10,
	 	'dt' : 1e-1,
		'err_IT' : 1e-9
	 	}

	else:
		args_init.update({'psi_init' : psi0})

	mu_all,counterSC,psi0,err_psi_SC = conv_mu_SC(args_syst,args_init)
	psi0,err_psi_IT,counterIT = conv_IT(args_syst,args_init)

	k += 1

pyplot.semilogx(np.arange(counterSC+1),err_psi_SC/np.max(err_psi_SC),'-',\
		label="SC,$U = %.1e$"%(U),fillstyle='none')
pyplot.semilogx(np.arange(counterIT+1),err_psi_IT/np.max(err_psi_IT),'-',\
		label="IT,$U = %.1e$"%(U),fillstyle='none')

pyplot.grid(axis='both')
pyplot.legend(fontsize=18)
pyplot.xlabel('\\# iterations',fontsize=20)
pyplot.ylabel('$\\delta_r/\\max(\\delta_r)$',fontsize=20)
pyplot.xticks(fontsize=16)
pyplot.yticks(fontsize=16)

pyplot.show()

temp = 'conv_mu_Us_%s_%s_%s_Nx_%.1i_J_%.3f_V_%.3e' %\
		(args_syst['Method'],args_syst['Trap'],\
		args_syst['Symm'],args_syst['Nx'],args_syst['J'],\
		args_syst['V'])

#path_above_Maxime = (sys.path[0]).split('Maxime')[0]
#fig_OL1D.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\Figures_OL_1D_BEC_SC'+temp+".pdf") # Windows
path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
fig_OL1D.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/fig_conv_SC'+temp+".pdf") # Linux


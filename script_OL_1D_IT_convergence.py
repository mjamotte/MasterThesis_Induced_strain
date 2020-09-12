import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from numpy import linalg
import cmath
import time
from matplotlib import rc

from numpy import linalg
import time

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

############### TEST Imaginary Time ###################

GenerateData = True
Plot = False
conv_mu = True

if GenerateData==True:

	params_interv = {'U': [1e-5],'V': [0.5e-4],'N': [1e4],\
	'dt': [0.1]}
	extension = '.npy'

	U = params_interv['U'][0]
	V = params_interv['V'][0]
	N = params_interv['N'][0]

	args_syst = {
	'J' : 1,
	'N' : N,
	'V' : V,
	'Nx' : 201,
	'U' : U,
	'Method' : 'IT',
	'Trap' : 'Harmonic',
	'Symm' : 'Isotropic',
	}

	args_syst.update({'sites_dic' : GP.lattice_1D(args_syst)})

	## Kinetic + Trap part of Hamiltonian

	H_KV = GP.H_1D(args_syst)
	eigVal,eigVecs = linalg.eigh(H_KV)
	E_OH_max = np.max(eigVal)
	psi_Gauss = eigVecs[:,0]

	for dt in params_interv['dt']:

		start = time.time()

		args_init = {
	 	'H_KV' : H_KV,
		'dt' : dt,
		'err_IT' : 1e-9,
		'psi_init' : psi_Gauss
		}
		mu_all,psi0,E_funct_IT = GP.calc_psi0(args_syst,args_init)

		print(linalg.eigh(H_KV+GP.H_int(psi0,args_syst))[0][:2])

		if U*N*np.max(np.abs(psi0)**2)<1 or True: # J = 1
			mu = mu_all[-1,int(len(mu_all[0])/2)]
			Trap = GP.trap_1D(args_syst)
			U = args_syst['U']
			N = args_syst['N']
			Nx = args_syst['Nx']
			n_TF = (mu-Trap)/U/N
			diff_TF_n0 = np.sum(np.abs(n_TF-abs(psi0)**2))

			data = [args_syst,args_init,mu,psi0,E_funct_IT,H_KV,diff_TF_n0]
			dataID = '1D_%s_%s_%s_Nx_%.1i_J_%.2f_N_%.3e_V_%.3e_U_%.3e_dt_%.3e' %\
					(args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
					args_syst['Nx'],args_syst['J'],N,args_syst['V'],args_syst['U'],\
					args_init['dt'])
			np.save('Datas_IT_conv/'+dataID,data)

			print("For Nx = ", Nx, ", N = ", N, ", U = ", U, ", V = ", V,\
				'dt = ', dt, "it took",time.time()-start,"secondes")

############# PLOTS ###################

if Plot==True:

	fig_OL1D = pyplot.figure(figsize=(8,8))
	ax1 = pyplot.subplot(111)
	ax1.set_xlabel('$\log(dt)$')
	ax1.set_ylabel('$\log(d\phi)$')

	## Get the data
	path = \
	'/home/maxime/Desktop/Codes/Optical_1D_Lattice_BEC/Imaginary_time/Datas_IT_conv/'
	list_args = ['N','V','U','dt']
	split = ['N_','V_','U_','dt_']
	params_interv = {'U': [1e-4],'V': [0.5e-4],'N': [1e4],\
						'dt': [5e-2,5e-1]}
	extension = '.npy'

	files_names = Mf.select_files(path,extension,list_args,params_interv,split)

	split = ['dt_']
	list_args = ['dt']
	files_names = Mf.reorder_filesnames(files_names,list_args,split)

	psi0_s = []
	dt_s = []
	for filename in files_names:

		data = np.load(filename+'.npy',allow_pickle=True)
		args_syst = data[0]
		J = args_syst['J']
		Nx = args_syst['Nx']
		V = args_syst['V']
		U = args_syst['U']
		N = args_syst['N']

		args_init = data[1]
		dt = args_init['dt']	
		mu = data[2]
		psi0 = data[3]
		E_funct_IT = data[4]

		psi0_s.append(psi0)
		dt_s.append(dt)

	psi0_s = np.array(psi0_s)
	dt_s = np.array(dt_s)

	## Compare to the maximum of the most precise psi
	err_max = []
	half = int(len(psi0_s[0])/2)
	
	for i in range(1,len(psi0_s)): ## avoid to compare two identic psi's
		err_max.append(abs(psi0_s[0,half]-psi0_s[i,half]))

	err_max = np.array(err_max)

	## Plot

	ax1.plot(np.log(dt_s[1:]),np.log(err_max),'.-')
	
	ax1.grid(axis='both')
	pyplot.suptitle('Test conv. IT: %s, %s,\
	 Nx = %.1i, N = %.3e, J = %.2f, V = %.3e, U = %.3e' % \
	 (args_syst['Trap'],args_syst['Symm'],\
	args_syst['Nx'],args_syst['N'],args_syst['J'],args_syst['V']\
	,args_syst['U']))

	pyplot.show()

	temp = 'conv_%s_%s_%s_Nx_%.1i_J_%.3f_V_%.3e_U_%.3e' %\
		(args_syst['Method'],args_syst['Trap'],\
		args_syst['Symm'],args_syst['Nx'],args_syst['J'],\
		args_syst['V'],args_syst['U'])
	fig_OL1D.savefig("Figures_OL_1D_BEC_IT_conv/fig_"+temp+".pdf")

#######################################
#######################################

def conv_mu_IT(args_syst,args_init):

	## Initialisation 
	mu_all = []
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

	while True:
		## time evolution
		solver.set_initial_value(psi_old, t)
		solver.integrate(t+dt)
		t = t+dt
		psi_new = solver.y

		## Compute mu
		psi_old_co = psi_old[:int(dim/2)] + 1j*psi_old[int(dim/2):]
		psi_new_co = psi_new[:int(dim/2)] + 1j*psi_new[int(dim/2):] # not renorm. yet

		## renormalize
		psi_new = psi_new/np.sqrt(np.sum(abs(psi_new)**2))
		psi_new_co = psi_new[:int(dim/2)] + 1j*psi_new[int(dim/2):]

		err_psi = np.sqrt(np.abs(np.max(np.abs(psi_new_co)**2\
					-np.max(np.abs(psi_old_co)**2))))\
					/np.sqrt(np.max(np.abs(psi_new_co)**2))

		mu_all.append(GrPi.compute_mu(psi_old_co,dt,args_syst,args_init))

		if err_psi<err:
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

	mu_all = np.array(mu_all)

	return mu_all,counterIT,psi0
'''
if conv_mu==True:
	k = 0
	for U in [1e-6]:
		args_syst = {
					'syst' : '1D',
					'J' : 1,
					'N' : 1e4,
					'V' : 0.5e-4,
					'Nx' : 201,
					'U' : U,
					'Method' : 'IT',
					'Trap' : 'Harmonic',
					'Symm' : 'Isotropic',
					}
		args_syst.update({'sites_dic' : GrPi.lattice_1D(args_syst)})
		H_KV = GrPi.H_1D(args_syst)

		if k==0:
		 	args_init = {
		 	'H_KV' : H_KV,
		 	'err_IT' : 1e-6,
		 	'dt' : 1e-1
		 	}

		else:
			args_init.update({'psi_init' : psi0})

		dt = args_init['dt']
		mu_all,counterIT,psi0 = conv_mu_IT(args_syst,args_init)
	
np.save('mu_all',mu_all)
'''
args_syst = {
			'syst' : '1D',
			'J' : 1,
			'N' : 1e4,
			'V' : 0.5e-4,
			'Nx' : 201,
			'U' : 1e-5,
			'Method' : 'IT',
			'Trap' : 'Harmonic',
			'Symm' : 'Isotropic',
			}
fig_OL1D = pyplot.figure(figsize=(8,8))
mu_all = np.load('mu_all.npy',allow_pickle=True)
# for x in range(len(mu_all[0])):
# 	pyplot.plot(np.arange(len(mu_all)),mu_all[:,x].real) # mu.imag is zero
pyplot.semilogx(np.arange(len(mu_all)),mu_all[:,100],'-',label="$U = %.1e$"%(args_syst['U']),fillstyle='none')

# pyplot.suptitle('Conv. $\mu$: %s, %s,\
# 	 Nx = %.1i, N = %.3e, J = %.2f, V = %.3e, U = %.3e' % \
# 	 (args_syst['Trap'],args_syst['Symm'],\
# 	args_syst['Nx'],args_syst['N'],args_syst['J'],args_syst['V']\
# 	,args_syst['U']))

pyplot.grid(axis='both')
pyplot.legend(fontsize=20)
pyplot.xlabel('\\# iterations',fontsize=18)
pyplot.ylabel('$\\mu$',fontsize=20)
pyplot.xticks(fontsize=16)
pyplot.yticks(fontsize=16)
pyplot.ylabel('$\mu$')
pyplot.show()

temp = 'conv_mu_%s_%s_%s_Nx_%.1i_J_%.3f_V_%.3e_U_%.3e' %\
		(args_syst['Method'],args_syst['Trap'],\
		args_syst['Symm'],args_syst['Nx'],args_syst['J'],\
		args_syst['V'],args_syst['U'])

#path_above_Maxime = (sys.path[0]).split('Maxime')[0]
#fig_OL1D.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\Figures_OL_1D_BEC_SC'+temp+".pdf") # Windows
path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
fig_OL1D.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/fig_conv_SC'+temp+".pdf") # Linux
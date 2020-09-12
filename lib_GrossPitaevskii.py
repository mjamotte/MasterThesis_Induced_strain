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
from scipy.integrate import ode
from numba import jit
import sys
sys.path.append('/home/maxime/Desktop/Codes/')
import lib_Manage_files as Mf
import lib_Optical_Honeycomb_Lattice as OHL
############# FUNCTIONS OPTICAL LATTICES #############

'''

^ y	  |	 |	|  |  |
|	  A	 B--A  B--A
|	  |  |  |  |  |
|	  B--A  B--A  B
|	  |  |  |  |  |
|	  A  B--A  B--A
|	(0,0)
|-------------------> x

'''

############# 1D OPTICAL LATTICE #############

''' 
	The goal of coding this is the understanding of the physics 
	of the Hamiltonian presented in the Lühmann's paper to 
	aquire feelings with correlated hoppings. 
'''

def lattice_1D(args_syst):
	Nx = args_syst['Nx']

	sites_dic = {}
	n = 0
	for i in range(Nx): # fill line by line, from left to right
		sites_dic[n] = np.array([n])
		n += 1
			
	return sites_dic

def trap_1D(args_syst):

	V0 = args_syst['V']
	sites_dic = args_syst['sites_dic']
	Nx = args_syst['Nx']
	x_max = np.max(OHL.extract_coord_from_dict(sites_dic)[0])

	x_center = x_max/2
	V = np.zeros(len(sites_dic)) 

	for n in range(len(sites_dic)):
		x = sites_dic[n][0]
		V[n] = 1/2*V0*(x-x_center)**2

	return V


def H_1D(args_syst):

	J = args_syst['J']
	Nx = args_syst['Nx']

	H = 1j*np.zeros((Nx,Nx))

	if args_syst['Symm']=='Isotropic':

		if args_syst['Trap']=='Harmonic':
			H = np.diag(trap_1D(args_syst))

		for i in range(Nx-1):
			H[i,i+1] = -J
			H[i+1,i] = np.conj(H[i,i+1])

	return H


def calc_psi0(args_syst,args_init):

	if args_syst['Method']=='SC': #SC = self-consistent method
		psi0 = solveGP_SC(args_syst,args_init)
		return psi0

	elif args_syst['Method']=='IT': #Imaginary time
		psi0 = solve_GP_IT(args_syst,args_init)
		#mu_all,psi0,E_funct_IT = solve_GP_IT_RK4(args_syst,args_init)
		return psi0


def solveGP_SC(args_syst,args_init):

	N = args_syst['N']
	U = args_syst['U']
	H_KV = args_init['H_KV']
	err = args_init['err_SC'] # condition to stop convergence
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
		H_U = H_int(psi_old,args_syst)
		eigVal,eigVecs = linalg.eigh(H_KV+H_U)
		psi_new = np.matrix.transpose(eigVecs[:,0])

		psi_lam = np.sqrt((1-lam)*psi_old**2 + lam*psi_new**2)

		err_psi.append(np.abs(np.max(np.abs(psi_lam)**2)\
					-np.max(np.abs(psi_old)**2))\
					/np.max(np.abs(psi_lam)**2))

		if err_psi[-1]<err:
			break

		psi_old = psi_lam
		counterSC += 1			

	#print('Number of iterations of self-consistent method =',counterSC)

	return psi_lam


def GP(t,psi_old,args_syst,args_init):

	psi_old = psi_old/np.sqrt(np.sum(np.abs(psi_old)**2))
	H_KV = args_init['H_KV']
	#H_KV = sc.sparse.coo_matrix(H_KV) # KV = Kinetic + Trap
	N = args_syst['N']
	dim = len(psi_old)
	psi_old_co = psi_old[:int(dim/2)] + 1j*psi_old[int(dim/2):]

	# Hopping part + trap part of the GP
	y1 = H_KV.dot(psi_old_co)
	# Interacting part of the GP
	y2 = H_int(psi_old_co,args_syst).dot(psi_old_co)

	y = y1 + y2
	# -d_tau psi(tau) = H psi(tau)
	y = -y
	return np.concatenate((np.real(y),np.imag(y)))



def compute_mu(psi_old,dt,args_syst,args_init):

	H_KV = args_init['H_KV']
	psi_new = sc.linalg.expm(-(H_KV+H_int(psi_old,args_syst))*dt).dot(psi_old)
	## Here, psi_new isn't renormalized

	mu = -np.log(psi_new/psi_old)/dt
	mu = -np.log(np.sqrt(np.sum(np.abs(psi_new)**2)))/dt

	return mu#,psi_new


def solve_GP_IT(args_syst,args_init):

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
	solver = ode(GP)
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

	return psi0


def H_int(psi,args_syst):

	U = args_syst['U']
	N = args_syst['N']
	return U*(N-1)*np.diag(np.abs(psi)**2)


def E_kin(psi,args_syst):

	J = args_syst['J']
	N = args_syst['N']
	sites_dic = args_syst['sites_dic']

	x_s = OHL.extract_coord_from_dict(sites_dic)[0]
	x_s,psi = OHL.reorder(x_s,psi)

	

	if args_syst['syst']=='Honey':

		if args_syst['BC']=='Yper':
			E_kin = 0
			for k in range(1,len(sites_dic)-1):
				if k%2==0:
					E_kin += -2*J*np.conj(psi[k])*psi[k+1]-J*psi[k]*np.conj(psi[k-1])
				if k%2!=0:
					E_kin += -J*np.conj(psi[k+1])*psi[k]-2*J*psi[k-1]*np.conj(psi[k])

			E_kin = N*E_kin

	if args_syst['syst']=='1D':
		E_kin = 0
		for k in range(len(sites_dic)-1):
			E_kin += -J*N*np.conj(psi[k])*psi[k+1]-J*N*np.conj(psi[k+1])*psi[k]

	return E_kin


def E_int(psi,args_syst):

	N = args_syst['N']
	U = args_syst['U']
	Nx = args_syst['Nx']

	E_U = U*N*(N-1)/2*np.sum(np.abs(psi)**4)  
	return E_U


def E_trap(psi,args_syst):

	N = args_syst['N']
	
	if args_syst['syst']=='1D':
		E_tr = N*np.sum(trap_1D(args_syst)*abs(psi)**2)
	if args_syst['syst']=='Honey':
		E_tr = N*np.sum(OHL.trap(args_syst)*abs(psi)**2)

	return E_tr


def energy_functional(psi,args_syst):

	E_U = E_int(psi,args_syst)
	E_tr = E_trap(psi,args_syst)
	E_k = E_kin(psi,args_syst)

	return E_U+E_tr+E_k


def impose_sym(vector): 
	## imposes symmetry from the middle of the vector
	Len = len(vector)

	for i in range(int(Len/2)):
		vector[i] = (vector[i]+vector[Len-i-1])*0.5
		vector[Len-i-1] = vector[i]

	norm = np.sqrt(np.sum(np.abs(vector)**2))
	vector = vector/norm
	return vector


def dEdN_O2(E,dN):#(path,N,dN,filesnames):

	dEdN = np.array([])
	for i in range(int(len(E)/3)): # 3 = "order of approx"+1
		dEdN = np.append(dEdN,(E[i*3+2]-E[i*3])/(2*dN))
	return dEdN


def gauss(xs,x0,sigma): # analytical

	Gauss = np.zeros(len(xs))
	i = 0
	for x in xs:
		Gauss[i] = 1/np.sqrt(sigma**2*2*np.pi)*np.exp(-(x-x0)**2/(2*sigma**2))
		i += 1

	return Gauss


def vector2matrix(vector,Nx,Ny):

	# Shape a vector into a Nx x Ny matrix
	L = len(vector)
	if L==Nx*Ny:
		mat = np.zeros((Ny,Nx))

		for i in range(Ny):
			for j in range(Nx):
				mat[i,j] = vector[i*Nx+j]

	else:
		print('Sizes do not match')
		mat = 0

	return mat


def extract_coord_from_dict(dictio):

	L = len(dictio)
	x_s = np.zeros(L)
	y_s = np.zeros(L) # revlevant for 2D systems

	for n in range(L):

		if len(dictio[n])==1:
			x_s[n] = dict[n]
		else:
			x_s[n],y_s[n] = dictio[n]

	return x_s,y_s


def size_cloud(args_syst,mu,n0):

	n0_copy = n0.copy()
	Lx = np.max(extract_coord_from_dict(args_syst['sites_dic'])[0])
	Nx = args_syst['Nx']
	Ny = args_syst['Ny']
	V = args_syst['V']
	J = args_syst['J']
	k = 0
	for i in range(0,len(n0),2):
		if n0_copy[i]>1e-5:
			k += 1

	R = k*3/2

	#if args_syst['syst']=='Honey' and args_syst['BC']=='Yper':
	#	R = np.sqrt(2*mu/V)*2 # distance between the two roots of the parabola

	return R,Lx
		
def mean_x(psi,positions):
	x_center = np.max(positions)/2
	return np.sum(np.abs(psi)**2*(positions-x_center))

def nearest_value(liste,x0):

	out = np.min(np.abs(liste))
	index = np.argmin(np.abs(liste))

	for i in range(len(liste)):
		x = liste[i]
		
		if abs(x0-x)<=abs(x0-out):
			out = liste[i]

			index = i

	return out,index

#A = np.array([-2.23470900e-01, -2.00863960e-01, -1.70701945e-01, -4.64587897e-16, -1.39059170e-16, 1.70701945e-01,  2.00863960e-01])
#print(nearest_value(A,0.001))

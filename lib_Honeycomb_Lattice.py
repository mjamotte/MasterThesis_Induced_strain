import numpy as np
import sympy as sy
from sympy import *
#from sympy import Matrix
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from numpy import linalg
#import cmath
import numba
from numba import jit

##############################
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

def square_lattice(Nx,Ny,a1,a2):

####################################################################################
#																				   #	
#	Constructs a square lattice with Bravais vectors a1 and a2, and composed of Nx # 
# 	sites in the x direction and Ny sites in the y direction.					   #
#																				   #	
####################################################################################
	 
	sites_dic = {}
	n = 0
	for j in range(Ny):
		for i in range(Nx): # fill line by line, from left to right
			sites_dic[n] = np.array([i,j])
			n += 1
			
	return sites_dic

def honey_lattice(Nx,Ny): # boundaries = zigzag

	####################################################################################
	#																				   #	
	#	Constructs a square lattice with Bravais vectors a1 and a2, and composed of Nx # 
	# 	sites in the x direction and Ny sites in the y direction.					   #
	#																				   #	
	####################################################################################

	sites_dic = {}
	n = 0
	x = 0

	for j in range(Ny):

		## create the scheme to repeat along y axis
		x = 0
		y = j*np.sqrt(3)/2

		for i in range(Nx): # fill line by line, from left to right
			
			if i%2==0:
				sites_dic[n] = np.array([x,y])
				x = x+2

			else:
				sites_dic[n] = np.array([x,y])
				x = x+1
				
			n += 1	

		if j%2!=0:
			for k in range(Nx):
				L = len(sites_dic)
				if k%2==0:
					sites_dic[L-Nx+k][0] += 0.5
				else:
					sites_dic[L-Nx+k][0] -= 0.5
			
	return sites_dic


def honey_hopping(n1,n2,sites_dic,system,delta_link):

####################################################################################
# 																				   #
#	Computes the Hamiltonian element H[n1,n2] where n1 and n2 are the index of the #
#	sites in the sites dictionnary; 								   #
#																				   #
####################################################################################
	
	t = system[3]
	tau = system[4]

	if system[0]=='realB_Ax': #the vector potential is parallel to x-axis
		Ny = system[2]
		Ly = Ny*np.sqrt(3)/2
		yA = sites_dic[n1][1]*np.sqrt(3)/2
		phi = 0.0346/2/np.pi

		if delta_link=='delta_1':
			hop = t*np.exp(-1j*2*np.pi*phi*(yA-Ly/2))

		elif delta_link=='delta_2':
			hop = t

		elif delta_link=='delta_3':
			hop = t

	if system[0]=='Yper_Lifsh':

		Nx = system[1]
		Lx = 3/2*(Nx-1) # Nx must be odd
		xA = (sites_dic[n1][0])*3/2

		if delta_link=='delta_1': # horizontal link, n1 is always even
			#hop = 2*t # to observe Lifshitz transition
			hop = t*1.5

		elif delta_link=='delta_2':
			hop = t

		elif delta_link=='delta_3':
			hop = t

	if system[0]=='Yper':

		Nx = system[1]
		Lx = 3/2*(Nx-1) # Nx must be odd
		xA = (sites_dic[n1][0])*3/2

		if delta_link=='delta_1': # horizontal link, n1 is always even
			#hop = 2*t # to observe Lifshitz transition
			hop = t*(1+tau/3*(xA-Lx/2))

		elif delta_link=='delta_2':
			hop = t

		elif delta_link=='delta_3':
			hop = t

	if system[0]=='Yper_opp_strain':

		Nx = system[1]
		Lx = 3/2*(Nx-1) # Nx must be odd
		xA = (sites_dic[n1][0])*3/2

		if delta_link=='delta_1': # horizontal link, n1 is always even
			#hop = 2*t # to observe Lifshitz transition
			hop = t*(1-tau/3*(xA-Lx/2))

		elif delta_link=='delta_2':
			hop = t

		elif delta_link=='delta_3':
			hop = t

	if system[0]=='Yper_abs_strain':

		Nx = system[1]
		Lx = 3/2*(Nx-1) # Nx must be odd
		xA = (sites_dic[n1][0])*3/2

		if delta_link=='delta_1': # horizontal link, n1 is always even
			#hop = 2*t # to observe Lifshitz transition
			hop = t*(1+tau/3*abs(xA-Lx/2))

		else:
			hop = t

	if system[0]=='Ypert2t3':

		Nx = system[1]
		Lx = 3/2*(Nx-1) # Nx must be odd
		xA = (sites_dic[n1][0])*3/2
		gamma2 = system[5]
		gamma3 = system[6]

		if delta_link=='delta_1': # horizontal link, n1 is always even
			#hop = 2*t # to observe Lifshitz transition
			hop = t*(1-tau/3*(xA-Lx/2))

		if delta_link=='delta_2':
			hop = t*(1-gamma2*tau/3*(xA-Lx/2))

		if delta_link=='delta_3':
			hop = t*(1-gamma3*tau/3*(xA-Lx/2))
	
	if system[0]=='Xper':

		yA = sites_dic[n1][1]*np.sqrt(3)/2
		Ny = system[2]
		Ly = Ny*np.sqrt(3)/2

		if delta_link=='delta_1':
			hop = t

		elif delta_link=='delta_2':
			hop = t*(1+tau*(yA-Ly/2)*np.sqrt(3)/3)

		elif delta_link=='delta_3':
			hop = t*(1-tau*(yA-Ly/2)*np.sqrt(3)/3)

	if system[0]=='Trigonal':

		xA = sites_dic[n1][0]*3/2
		yA = sites_dic[n1][1]*np.sqrt(3)/2
		Nx = system[1]
		Ny = system[2]
		Lx = 3/2*(Nx-1)
		Ly = Ny*np.sqrt(3)/2
		gamma2 = system[5]
		gamma3 = system[6]

		if delta_link=='delta_1':
			hop = t*(1-tau/3*(xA-Lx/2))

		#elif delta_link=='delta_2':
		#	hop = t*(1+tau*(yA-Ly/2)*np.sqrt(3)/3)

		#elif delta_link=='delta_3':
		#	hop = t*(1-tau*(yA-Ly/2)*np.sqrt(3)/3)

		elif delta_link=='delta_2':
			hop = t*(1+gamma2*tau*(1/2*(xA-Lx/2)+(yA-Ly/2)*np.sqrt(3)/2)/3)

		elif delta_link=='delta_3':
			hop = t*(1+gamma3*tau*(1/2*(xA-Lx/2)-(yA-Ly/2)*np.sqrt(3)/2)/3)

	return hop


def H_honey(system,kxa,kya,sites_dic):

	d1 = np.array([-1,0])
	d2 = np.array([0.5,np.sqrt(3)/2])
	d3 = np.array([0.5,-np.sqrt(3)/2])
	deltas = [d1,d2,d3]
	Nx = system[1]
	Ny = system[2]
	t = system[3]
	tau = system[4]
	
	
	if system[0]=='open':
		H = 1j*np.zeros((Nx*Ny,Nx*Ny))

		for n in range(Nx*Ny):
			if n%Nx!=0 and n%2==0: # this A isn't on the left boundry
				H[n,n-1] = honey_hopping(n,n-1,sites_dic,system)
				H[n-1,n] = np.conj(H[n,n-1])

			if (n<Nx*(Ny-1)): # this A is not on the upper boundary
				H[n,n+Nx] = honey_hopping(n,n+Nx,sites_dic,system)
				H[n+Nx,n] = np.conj(H[n,n+Nx])

	if system[0]=='realB_Ax':
		H = 1j*np.zeros((Nx*Ny,Nx*Ny))
		
		for n in range(0,Nx*Ny,2):

			if n%(2*Nx)==0:
				# n = A site to n+1 = B site
				H[n,n+1] = honey_hopping(n,n+1,sites_dic,system,'delta_1')*np.exp(+1j*kxa*3) # link t1 
				H[n+1,n] = np.conj(H[n,n+1])

				if n<Nx*(Ny-1):
					# n = A site to n+Nx = B site
					H[n,n+Nx] = honey_hopping(n,n+Nx,sites_dic,system,'delta_2') # link t2
					H[n+Nx,n] = np.conj(H[n,n+Nx])
					
					# n+1 = B site to n+1+Nx = A site
					H[n+1+Nx,n+1] = honey_hopping(n+1+Nx,n+1,sites_dic,system,'delta_3') # link t3
					H[n+1,n+1+Nx] = np.conj(H[n+1+Nx,n+1])


			elif n%(2*Nx)!=0:
				# n = B site to n+1 = A site
				H[n+1,n] = honey_hopping(n+1,n,sites_dic,system,'delta_1')
				H[n,n+1] = np.conj(H[n+1,n])

				if n<Nx*(Ny-1):
					# n = B site to n+Nx = A site
					H[n+Nx,n] = honey_hopping(n+Nx,n,sites_dic,system,'delta_3') # link t3
					H[n,n+Nx] = np.conj(H[n+Nx,n])

					# n+1 = A site to n+1+Nx = B site
					H[n+1,n+1+Nx] = honey_hopping(n+1,n+1+Nx,sites_dic,system,'delta_2') # link t2
					H[n+1+Nx,n+1] = np.conj(H[n+1,n+1+Nx])

	if system[0]=='Yper' or system[0]=='Ypert2t3' or 'Yper_opp_strain' or 'Yper_abs_strain':
		H = 1j*np.zeros((Nx*Ny,Nx*Ny))
		
		for n in range(Nx*Ny):
			if n%Nx!=0 and n%2==0 and n>0: # this A isn't on the left boundry
				H[n,n-1] = honey_hopping(n,n-1,sites_dic,system,'delta_1')
				H[n-1,n] = np.conj(H[n,n-1])

			if n<Nx*(Ny-1) and n%2==0:  # this A is not on the top boundary
				H[n,n+Nx] = honey_hopping(n,n+Nx,sites_dic,system,'delta_2')*np.exp(1j*kya*deltas[1][1])+\
							+honey_hopping(n,n+Nx,sites_dic,system,'delta_3')*np.exp(1j*kya*deltas[2][1])	
				H[n+Nx,n] = np.conj(H[n,n+Nx])

			if n<Nx*(Ny-1) and n%2!=0:  # this B is not on the top boundary
				H[n,n+Nx] = honey_hopping(n,n+Nx,sites_dic,system,'delta_3')*np.exp(-1j*kya*deltas[2][1])+\
							+honey_hopping(n,n+Nx,sites_dic,system,'delta_2')*np.exp(-1j*kya*deltas[1][1])
				H[n+Nx,n] = np.conj(H[n,n+Nx])

	if system[0]=='Xper': # ARMCHAIR
		H = 1j*np.zeros((Nx*Ny,Nx*Ny))

		for n in range(0,Nx*Ny,2):

			if n%(2*Nx)==0:
				# n = A site to n+1 = B site
				H[n,n+1] = honey_hopping(n,n+1,sites_dic,system,'delta_1')*np.exp(-1j*kxa*2) # link t1 
				H[n+1,n] = np.conj(H[n,n+1])

				if n<Nx*(Ny-1):
					# n = A site to n+Nx = B site
					H[n,n+Nx] = honey_hopping(n,n+Nx,sites_dic,system,'delta_2') # link t2
					H[n+Nx,n] = np.conj(H[n,n+Nx])
					
					# n+1 = B site to n+1+Nx = A site
					H[n+1+Nx,n+1] = honey_hopping(n+1+Nx,n+1,sites_dic,system,'delta_3') # link t3
					H[n+1,n+1+Nx] = np.conj(H[n+1+Nx,n+1])


			elif n%(2*Nx)!=0:
				# n = B site to n+1 = A site
				H[n+1,n] = honey_hopping(n+1,n,sites_dic,system,'delta_1')*np.exp(-1j*kxa)
				H[n,n+1] = np.conj(H[n+1,n])

				if n<Nx*(Ny-1):
					# n = B site to n+Nx = A site
					H[n+Nx,n] = honey_hopping(n+Nx,n,sites_dic,system,'delta_3') # link t3
					H[n,n+Nx] = np.conj(H[n+Nx,n])

					# n+1 = A site to n+1+Nx = B site
					H[n+1,n+1+Nx] = honey_hopping(n+1,n+1+Nx,sites_dic,system,'delta_2') # link t2
					H[n+1+Nx,n+1] = np.conj(H[n+1,n+1+Nx])

	if system[0]=='Trigonal':

		H = np.zeros((Nx*Ny,Nx*Ny))

		for n in range(Nx*Ny-1):

			if n%(2*Nx)<Nx: # line beginning with A site

				if n%2!=0 and (n+1)%Nx!=0: # not on the right boarder
					# n = B site to n+1 = A site
					H[n+1,n] = honey_hopping(n+1,n,sites_dic,system,'delta_1') #link t1
					H[n,n+1] = np.conj(H[n+1,n])

				if n<Nx*(Ny-1) and n%2==0:
					# n = A site to n+Nx = B site
					H[n,n+Nx] = honey_hopping(n,n+Nx,sites_dic,system,'delta_2') #link t2
					H[n+Nx,n] = np.conj(H[n,n+Nx])

					if n<Nx*(Ny-1)-1:
						# n+1 = B site to n+1+Nx = A site
						H[n+1+Nx,n+1] = honey_hopping(n+1+Nx,n+1,sites_dic,system,'delta_3') # link t3
						H[n+1,n+1+Nx] = np.conj(H[n+1+Nx,n+1])

			elif n%(2*Nx)>=Nx: # line beginning with B site

				if Nx%2!=0: # line begins and ends with same site-type

					if n%2!=0 and (n+1)%Nx!=0:
						# n = B site to n+1 = A site
						H[n+1,n] = honey_hopping(n+1,n,sites_dic,system,'delta_1')
						H[n,n+1] = np.conj(H[n+1,n])

					if n<Nx*(Ny-1) and n%2!=0:

						# n = B site to n+Nx = A site
						H[n+Nx,n] = honey_hopping(n+Nx,n,sites_dic,system,'delta_3')
						H[n,n+Nx] = np.conj(H[n+Nx,n])

						if n<Nx*(Ny-1)-1:
							# n+1 = A site to n+1+Nx = B site
							H[n+1,n+1+Nx] = honey_hopping(n+1,n+1+Nx,sites_dic,system,'delta_2')
							H[n+1+Nx,n+1] = np.conj(H[n+1,n+1+Nx])

				if Nx%2==0: # line doesn't begin and end with same site-type
				
					if n%2==0 and (n+1)%Nx!=0:
						# n = B site to n+1 = A site
						H[n+1,n] = honey_hopping(n+1,n,sites_dic,system,'delta_1')
						H[n,n+1] = np.conj(H[n+1,n])

					if n<Nx*(Ny-1) and n%2==0:

						# n = B site to n+Nx = A site
						H[n+Nx,n] = honey_hopping(n+Nx,n,sites_dic,system,'delta_3')
						H[n,n+Nx] = np.conj(H[n+Nx,n])

						if n<Nx*(Ny-1)-1:
							# n+1 = A site to n+1+Nx = B site
							H[n+1,n+1+Nx] = honey_hopping(n+1,n+1+Nx,sites_dic,system,'delta_2')
							H[n+1+Nx,n+1] = np.conj(H[n+1,n+1+Nx])


	if system[0]=='FullPeriodic':
		ka = np.array([kxa,kya])
		H = 1j*np.zeros((2,2))

		H[0,1] = t*(np.exp(-1j*np.vdot(ka,deltas[0]))+np.exp(-1j*np.vdot(ka,deltas[1]))\
					+np.exp(-1j*np.vdot(ka,deltas[2])))
		H[1,0] = np.conj(H[0,1])

	return H	

def E_relat_LL(n,system,kxa,kya,DP): # n is the number of the level
	
	t = system[3]
	tau = system[4]

	if system[0]=='Yper' or system[0]=='Ypert2t3' or 'Yper_opp_strain' or 'Yper_abs_strain':
		E_p = +t*np.sqrt(abs(n)*tau)*np.sqrt(1-DP[2]*(kya-DP[1]))
		#E_p = t*np.sqrt(n*tau)*np.sqrt(2-2*(kya-DP[1])+\
		#	n**2*tau/18+tau/12+(kya-DP[1])**2*(n+1/2))/np.sqrt(2)
		#E_p = np.sqrt(abs(n)*tau)*(1+3/2*(kya-DP[1])**3/np.sqrt(tau))
		E_m = -E_p

	if system[0]=='Xper':
		E_p = +t*np.sqrt(3*abs(n)*tau)*np.sqrt(1)
		E_m = -E_p

	if system[0]=='Trigonal':
		E_p = t*np.sqrt(3*tau*n)
		E_m = -E_p

	return E_p,E_m


def nearest_value(liste,x0):

	out = np.min(np.abs(liste))
	index = np.argmin(np.abs(liste))

	for i in range(len(liste)):
		x = liste[i]
		
		if abs(x0-x)<=abs(x0-out):
			out = liste[i]

			index = i

	return out,index

def LL_eigVect(x_list,N,x0,alpha):
	out = np.zeros(len(x_list))

	N = abs(N) # N could be negative 
	x = sy.symbols('x')
	eigVect = sy.exp(-((x-x0)/alpha)**2/2)*(-1)**N*sy.exp(((x-x0)/alpha)**2)*\
				sy.diff(sy.exp(-((x-x0)/alpha)**2),x,N)
	eigVect = lambdify(x,eigVect,"numpy")

	i = 0
	for x in x_list:
		out[i] = eigVect(x)
		i += 1
	return out/np.sqrt(np.sum(np.abs(out)**2))

# psi = LL_eigVect(np.arange(0,100,1),1,50,15)
# pyplot.plot(np.arange(0,100,1),np.abs(psi))
# pyplot.show()
# exit()

def polynom(x_list,coeff):

	out = np.zeros(len(x_list))

	i = 0
	for x in x_list:
		for n in range(len(coeff)):
			out[i] += coeff[len(coeff)-1-n]*x**n
		
		i += 1
	return out


def extract_coord_from_dict(dictio):

	L = len(dictio)
	x_s = np.zeros(L)
	y_s = np.zeros(L) # revlevant for 2D systems

	for n in range(L):

		if len(dictio[n])==1:
			x_s[n] = dictio[n]
		else:
			x_s[n],y_s[n] = dictio[n]

	return x_s,y_s



def reorder(x_s,psi): 
	'''
	The "x" coord. are reordered from smallest to greatest and to avoid 
	to lose the correspondant psi, we conjointly keep the couples (x,psi(x))
	intact.
	'''
	
	Len = len(x_s)
	indices_values = []
	for i,value in enumerate(x_s):
		index_value = [i,value]
		indices_values.append(index_value)

	copy = np.array(indices_values.copy())
	temp = np.zeros((Len,2))
	out_x = x_s.copy()
	out_psi = psi.copy()
	for i in range(Len):
		index_min = int(np.argmin(copy[:,1]))
		temp[i] = indices_values[index_min]
		index_max = int(np.argmax(copy[:,1]))
		copy[index_min,1] = copy[index_max][1]+1

	temp = np.array(temp) ## to be able to use int()
	for i in range(Len):
		out_x[i] = x_s[int(temp[i,0])]
		out_psi[i] = psi[int(temp[i,0])]

	return out_x,out_psi



'''
class site:

	Type = site_type # 'A' or 'B
	def AorB(self,i,j,Type):

		if (i%2==0 and j%2==0) or (abs(i%2)==1 and abs(j%2)==1):
			self.Type = 'A'
		else:
			self.Type = 'B'

	def nneigh_hop(self,m,t,tau):
		if (m==0):



	def index_neigh(self,nA,Nx,Ny): 

		 Returns the indices of the nearest neighbors of the nA^th site 
		   A of the lattice.

			Inputs: Nx,Ny = (integers) numbers of sites in x,y directions
					nA = (integer) 
			Ouputs: List of indices (integers)
		

		neigh_nA = np.array([])
		if (nA>Nx): # nA is not on the bottom boundary
			self.neigh_nA.append(nA-Nx)
			self.nneigh_hop(2)

		if (nA<Nx*Ny): # nA is not on the top boundary
			self.neigh_nA.append(nA+Nx)
			self.nneigh_hop(1)
			
		if (nA%Nx!=0): # nA is not on the left boundary
			self.neigh_nA.append(nA-1)
			self.nneigh_hop(0)
		
'''
	
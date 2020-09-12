import numpy as np
import sympy as sy
from sympy import *
#from sympy import Matrix
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from numpy import linalg
#import cmath
#from numba import jit

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

def honey_hopping_realB(n1,n2,sites_dic,system,link):

	if system[0]=='realB_open':

		xA = sites_dic[n1][0]*3/2
		yA = sites_dic[n1][1]*np.sqrt(3)/2
		Nx = system[1]
		Ny = system[2]
		Ax = system[7]
		Ay = system[8]
		Lx = 3/2*(Nx-1)
		Ly = Ny*np.sqrt(3)/2

		if link=='delta_1':
			hop = t*np.exp(1j)

		elif link=='delta_2':
			hop = t

		elif link=='delta_3':
			hop = t





def H_honey_realB(system,kxa,kya,sites_dic):

	d1 = np.array([-1,0])
	d2 = np.array([0.5,np.sqrt(3)/2])
	d3 = np.array([0.5,-np.sqrt(3)/2])
	deltas = [d1,d2,d3]
	Nx = system[1]
	Ny = system[2]
	t = system[3]
	tau = system[4]

	if system[0]=='realB_open': # Open system

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
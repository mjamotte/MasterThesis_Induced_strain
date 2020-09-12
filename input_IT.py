import numpy as np
import scipy as sc
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from numpy import linalg
import time

import sys
sys.path.append('/home/maxime/Desktop/Codes/')
import lib_GrossPitaevskii as GP
import lib_Manage_files as Mf
###############################################

## What to do:
testConvIT = {
				'Run' : True,
				'GenerateData' : True,
				'Plot' : True
			}

energIT = {
			'Run' : False,
			'GenerateData' : True,
			'Plot' : True
			}

densIT = {
			'Run' : False,
			'GenerateData' : True,
			'Plot' : True
			}

list_args = ['N','V','U']
split = ['N_','V_','U_']
params_interv = {'U': np.logspace(-5,-4,25),'V': [1e-5],\
					'N': [1e4]}
extension = '.npy'

if testConvIT['Run']==True:

	if testConvIT['GenerateData']==True:


	if testConvIT['Plot']==True:
		path = '/home/maxime/Desktop/Codes/Optical_1D_Lattice_BEC/Imaginary_time/Datas_IT_conv/'
		files_names = Mf.select_files(path,extension,list_args,params_interv,split)

if energIT['Run']==True:

	if energIT['GenerateData']==True:

	if energIT['Plot']==True:

if densIT['Run']==True:

	if densIT['GenerateData']==True:

	if densIT['Plot']==True:
#!/usr/bin/env python
# encoding: utf-8
#
# One dimensional solver for linear 1D Maxwell's equation via semi-spectral methods
import numpy as np
from params import Params as parameters
from geometry import Dimension
from geometry import Grid

fft = np.fft.fft
z = 1j

def qinit(grid,parameters):
    """
    Initial conditions in simulation grid for electromagnetic components q
    """
    x = grid.x_spectral.grid_points
    if parameters.source_type=='off':
        dd = np.max(x) - np.min(x)
        q[0,:] = 0.0
        q[1,:] = fft(parameters.source_amp[1]*np.exp(-(x)**2/(dd/10)**2)*np.sin(parameters.source_k*x)) 
    else:
        q[0,:] = 0.0
        q[1,:] = 0.0
	
	return q

# import parameters and set the default values
params = parameters(aux_shape='homogeneous',num_dim=1,num_aux=2)

# create dimension object and assign grid properties
X = Dimension(params)
grid = Grid(X,method='spectral')
ratio = grid._spectral_ratio
params.source_lambda = 1e-6*ratio
params.set_source()

q = qinit(grid,params)
params.solver_dt = 0.01



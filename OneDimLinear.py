#!/usr/bin/env python
# encoding: utf-8
#
# One dimensional solver for linear 1D Maxwell's equation via semi-spectral methods
import numpy as np
from params import Params as parameters
from geometry import Dimension
from geometry import Grid
from solver import State as state
from solver import Solver


def qinit(state):
    """
    Initial conditions in simulation grid for electromagnetic components q
    """
    x = state.grid.x_spectral.grid_points
    if state.parameters.source_type=='off':
        dd = np.max(x) - np.min(x)
        q[0,:] = 0.0
        q[1,:] = fft(state.parameters.source_amp[1]*np.exp(-(x)**2/(dd/10)**2)*np.sin(state.parameters.source_k*x)) 
    else:
        q[0,:] = 0.0
        q[1,:] = 0.0
	
	return q

# import parameters and set the default values
params = parameters(aux_shape='homogeneous',num_dim=1,num_aux=2)
params.dimension_lower[0] = 0.0
params.dimension_upper[0] = 100e-6
params.dimension_resolution[0] = 5 # set resolution in points per wavelength
params.dtcfl=0
params.set_dims()

# create dimension object and assign grid properties
X = Dimension(params)
state.grid = Grid(X,method='spectral')

# Update source to work correctly with spectral grid
params.source_lambda = 1e-6*state.grid._spectral_ratio
params.t_final = 1
params.solver_dt = 0.01
params.set_source()
state.parameters = params

# set initial conditions
state.q = qinit(state)


# instantiate solver
solver = Solver(state)

# set the solver kernel
solver.kernel = 'python'



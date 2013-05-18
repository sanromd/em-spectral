#!/usr/bin/env python
# encoding: utf-8
#
# One dimensional solver for linear 1D Maxwell's equation via semi-spectral methods
import numpy as np
from params import Params as parameters
from geometry import Dimension
from geometry import Grid
from state import State
from solver import Solver


def qinit(state):
    """
    Initial conditions in simulation grid for electromagnetic components q
    """
    x = state.grid.x_spectral.grid_points
    state.q = np.zeros([2,len(x)],dtype=np.complex)
    dd = np.max(x) - np.min(x)
    print dd, state.parameters.source_k, state.parameters.source_amp[1]
    print state.q
    if state.parameters.source_type=='off':
        state.q[0,:] = 0.0
        state.q[1,:] = np.fft.fft(state.parameters.source_amp[1]*np.exp(-(x)**2/(dd/10)**2)*np.sin(state.parameters.source_k*x)) 
    else:
        state.q[0,:] = 0.0
        state.q[1,:] = 0.0
    
    return state

# import parameters and set the default values
params = parameters(aux_shape='homogeneous',num_dim=1,num_aux=2)
params.dimension_lower[0] = 0.0
params.dimension_upper[0] = 100e-6
params.dimension_resolution[0] = 5 # set resolution in points per wavelength
params.dtcfl=0
params.set_dims()
state = State()
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
qinit(state)

# instantiate solver
solver = Solver(state)
print solver._N, solver._k
# set the solver kernel
solver.kernel = 'python'


def prestep1(solver,solution):
    print solver,solution

solver.add_prestep(prestep1)


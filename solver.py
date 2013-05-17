#!/usr/bin/env python
# encoding: utf-8
#
import numpy as np

class State(object):

    _initialized = 1 # sets the initialized flag
    t = 0 #initializes at time t = 0 by default
    parameters  = None
    q           = None
    dimensions  = None
    grid        = None
    solver      = None
    aux         = None
    _q_old      = None
    def __init__(self, **kargs):
        cls_dict = self.__class__.__dict__
        for key in kargs:
            assert key in cls_dict
        self.__dict__.update(kargs)

    def etar(self,t=0):


class Solver(object):

    def before_each_step(self,solution):
        r"""
        dummy function user must replace it by the desired calculation before each step.
        """
        pass


    kernel = 'python'
    update_before_each_step = 0 
    _z = 1j

    def __init__(self,State,**kargs):
        cls_dict = self.__class__.__dict__
        for key in kargs:
            assert key in cls_dict
        self.__dict__.update(kargs)

        self.num_t_steps = np.floor(State.parameters.t_final/State.parameters.solver_dt)
        self.dt = State.parameters.solver_dt        
        self._q = State.q
        self._q_temp = State.q
        self._t_ini = State.t
        self.t = np.linspace(self._t_ini,State.parameters.t_final,self.num_t_steps)
        self.tolerance = State.parameters.solver_tolerance

        if State.grid.num_dim==1:
            self.Lx = 2.0
            self._N = State.grid.num_points[1]
            self._kr = self._z*np.linspace(0,grid.num_points/2,grid.num_points/2+1)
            self._kl = self._z*np.linspace(-grid.num_points/2,0,grid.num_points/2+1)
        
        if state.aux_keep_copy==1:
            setattr(State,'_aux_old',None)


    def run(self,state):
        for n in range(1,self.num_steps+1):
            self.step(state)


    def step(self,state):
        r"""
        advances the solution one step. If a method before_each_step is present the routine computes the desired function before running the next step. 
        """
        if self.update_before_each_step==1:
            state._aux_old = state.aux
            self.before_each_step(self,state)
        else:
            state._aux_old = state.aux

        if self.kernel=='python'
            self.python_spectral(state)

    def python_spectral(self,state):
        fft = np.fft.fft
        ifft = np.fft.ifft
        z = 1j
        difference = 1
        chi2 = chi3 = 0
        q_old = q_calc = q_temp = state.q
        while difference>self.tolerance:
            q_mean = 0.5*(state.q + state._q_old)
            q_calc[0] =(state._aux_old[1]*q_old[0]+ self.dt*ifft(self.k*fft(q_mean[1]))
                        +(2*chi2*q_mean[0] + 3*chi3*q_mean[0]**2)*q_old[0])/(state.aux[1] + 2*chi2*q_mean[0] + 3*chi3*q_mean[0]**2)
            q_calc[1] =(state._aux_old[1]*q_old[1] + self.dt*ifft(self.k*fft(q_mean[0])
                        +(2*chi2*q_mean[1] + 3*chi3*q_mean[1]**2)*q_old[1])/(state.aux[1] + 2*chi2*q_mean[1] + 3*chi3*q_mean[1]**2) 
            difference = np.norm(q_temp - q_calc)
            q_temp = q_calc

        state.q = q_calc
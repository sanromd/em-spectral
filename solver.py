#!/usr/bin/env python
# encoding: utf-8
#
import numpy as np

class State(object):

    _initialized = 1 # sets the initialized flag
    t = 0 #initializes at time t = 0 by default
    parameters = None
    q = None
    dimensions = None
    grid = None
    solver = None
    aux = None
    _q_new = None
    _aux_old = None
    def __init__(self, **kargs):
        cls_dict = self.__class__.__dict__
        for key in kargs:
            assert key in cls_dict
        self.__dict__.update(kargs)

    def etar(self,t=0):


class Solver(object):

    def before_each_step(self):
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
        self._t_ini = State.t
        self.t = np.linspace(self._t_ini,State.parameters.t_final,self.num_t_steps)
        self.tolerance = State.parameters.solver_tolerance
        if State.grid.num_dim==1:
            self.Lx = 2
            self._N = State.grid.num_points[1]
            self._kr = self._z*np.linspace(0,grid.num_points/2,grid.num_points/2+1)
            self._kl = self._z*np.linspace(-grid.num_points/2,0,grid.num_points/2+1)

    def run(self,state):
        for n in range(1,self.num_steps+1):
            self.adavance_t(state)


    def adavance_t(self,state):
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
        chi1 = chi2 = chi3 = 0
        state._q[0] = state.q[0]
        state._q[1] = state.q[1]
        electric = magnetic = 0
        while difference>self.tolerance:
            electrictemp = electric
            magnetictemp = magnetic
            electricmean = 0.5*(electric + state._q[0])
            magneticmean = 0.5*(magnetic + state._q[1])
            electric=(state._aux_old[1,:] * state._q[0] + self.dt*ifft(self.k*fft(magneticmean)) 
                    +(2*chi2*electricmean + 3*chi3*electricmean^2)*state._q[0])/(state.aux[1]+2*chi2*electricmean+3*chi3*electricmean.^2)
            magnetic=(state._aux_old[1,:] * state._q[1] + self.dt*ifft(self.k.*fft(electricmean))
                    +(2*chi2*magneticmean+3*chi3*magneticmean.^2).*state._q[1])./(state.aux[1]+2*chi2*magneticmean+3*chi3*magneticmean.^2)
            difference=max(abs(electrictemp-electric))+max(abs(magnetictemp-magnetic))

        state.q[0] = electric
        state.q[1] = magnetic

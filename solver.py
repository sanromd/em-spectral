#!/usr/bin/env python
# encoding: utf-8
#
import numpy as np
class Solver(object):

    kernel = 'python'
    update_before_each_step = 0
    aux_keep_copy = 0 
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
        print State.grid.num_dim
        if State.grid.num_dim==1:
            self.Lx = 2.0
            self._N = int(State.grid.num_points[0])
            print self._N
            if np.mod(self._N,2)==0:
                self._I = complex(0,1)
                self._k = np.array([self._I*n for n in range(0,self._N/2) + [0] + range(-self._N/2+1,0)])
            else:
                raise Exception('The number of points should be even')
                
        if self.aux_keep_copy==1:
            setattr(State,'_aux_old',None)

        self.prestep_callbacks = []

    def add_prestep(self, prestep):
        self.prestep_callbacks.append(prestep)

    def run(self,state):
        for n in range(1,self.num_steps+1):
            self.step(state)


    def step(self,state):
        r"""
        advances the solution one step. If a method before_each_step is present the routine computes the desired function before running the next step. 
        """
        if self.update_before_each_step==1:
            state._aux_old = state.aux
            #self.before_each_step(self,state)
            for prestep in self.prestep_callbacks:
                prestep(self, state)
        else:
            state._aux_old = state.aux

        if self.kernel=='python':
            self.SpectralSolver1D(state)

    def SpectralSolver1D(self,state):
        # define some variables and functions
        fft = np.fft.fft
        ifft = np.fft.ifft
        k = self._k
        difference = 1
        chi2 = chi3 = 0
        q_old = q_calc = q_temp = state.q
        while difference>=self.tolerance:
            q_mean = 0.5*(state.q + q_old)
            q_calc[0] = (state._aux_old[1]*q_old[0] + self.dt*ifft(k*fft(q_mean[1])) + (2*chi2*q_mean[0] + 3*chi3*q_mean[0]**2)*q_old[0])/(state.aux[1] + 2*chi2*q_mean[0] + 3*chi3*q_mean[0]**2)
            q_calc[1] = (state._aux_old[1]*q_old[1] + self.dt*ifft(k*fft(q_mean[0])) + (2*chi2*q_mean[1] + 3*chi3*q_mean[1]**2)*q_old[1])/(state.aux[1] + 2*chi2*q_mean[1] + 3*chi3*q_mean[1]**2)
            q_temp = np.abs(q_temp - q_calc)
            difference = np.max(q_temp)
            q_temp = q_calc

        state.q = q_calc

    def SpectralSolver2D(self,state):
        pass

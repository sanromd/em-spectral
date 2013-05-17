#!/usr/bin/env python
# encoding: utf-8
# -------- GLOBAL SCALAR DEFINITIONS -----------------------------
# ======== all definitions are in m,s,g unit system.
import numpy as np
class Params(object):
    """
    Class: Params,  
    Instances the required simulation parameters, it uses a list of values from to the defaults file.
    """
    @property
    def aux_axis(self):
        return np.diagonal(self.aux_base)

    @property
    def ref_ind(self):
        if self.aux_tensor_kind=='isotropic':
            self.ref_ind = np.sqrt(self.aux_base[0,0]*self.aux_base[2,2])
        elif self.aux_tensor_kind=='anisotropic':
            self.ref_ind[0] = np.sqrt(self.aux_base[0,0]*self.aux_base[2,2])
            self.ref_ind[1] = np.sqrt(self.aux_base[1,1]*self.aux_base[2,2])
            self.ref_ind[2] = np.sqrt(self.aux_base[0,1]*self.aux_base[2,2])
            self.ref_ind[3] = np.sqrt(self.aux_base[1,0]*self.aux_base[2,2])
        elif self.aux_tensor_kind=='bianisotropic':
            self.ref_ind[0] = np.sqrt(self.aux_base[0,0]*self.aux_base[2,2])
            self.ref_ind[1] = np.sqrt(self.aux_base[1,1]*self.aux_base[2,2])
            self.ref_ind[2] = np.sqrt(self.aux_base[0,1]*self.aux_base[2,2])
            self.ref_ind[3] = np.sqrt(self.aux_base[1,0]*self.aux_base[2,2])

        return self.ref_ind

    num_frames = 30
    max_steps = 250000
    solver_cfl = 0.9
    # ....... dimensions .............................................
    num_dim = 2
    solver_dt = .1
    # ....... aux config .............................................
    num_aux = 3
    mode = 'TE'
    vacuum_config = 'real'
    aux_shape = 'homogeneous'
    aux_tensor_kind = 'isotropic'
    aux_base = np.zeros([num_aux,num_aux])
    _is_rip = 1
    _nlayers = 2
    _Nlayers = 5
    # ........ Boundary settings .................
    bc_lower = ['scattering', 'none']
    bc_upper = ['none', 'none']
    aux_bc_lower = ['scattering','pml']
    aux_bc_upper = ['metallic','metallic']
    solver_tolerance = 1e-6

    # parameters needed for pml calculation
    pml_points = 8
    pml_norder = 3
    pml_Ro = 1.0e-6

    # ........ excitation - initial conditoons .......................
    source_lambda  = 1e-6 
    source_type  = 'off'



    def __init__(self, **kargs):
        cls_dict = self.__class__.__dict__
        for key in kargs:
            assert key in cls_dict
        self.__dict__.update(kargs)
        self.set_vacuum(vc_config=self.vacuum_config)
        self.set_material(shape=self.aux_shape)
        self.set_source()
        self.set_dims(shape=self.aux_shape)
        self._initialized = 1

    def set_source(self):
        if '_initialized' in self.__dict__:
            pass
        else:
            self.source_width = 10.0 # width in the y-direction 
            self.source_t_offset  = 0.0
            self.source_offset  = [0.0, 0.0]      # offset in the y-direction
            self.source_time_width = 1.0*self.source_lambda   
            if self.num_dim==2:
                self.source_amp = np.zeros([3])
                self.source_amp[0] = 0.0                           # Amplitude of q1, q2 and q3
                self.source_amp[1] = 1.0
                self.source_amp[2] = 1.0
            elif self.num_dim==1:
                self.source_amp = np.zeros([2])
                self.source_amp[0] = 0.0                           # Amplitude of q1, q2 and q3
                self.source_amp[1] = 1.0

        self.source_k        = 2.0*np.pi/self.source_lambda       # wave vector        
        self.source_omega = 2.0*np.pi*self.co/self.source_lambda
        return self    

    def set_dims(self,shape=aux_shape,dtcfl=1):
        """
        set_dims(shape=aux_shape)
        calculates y_upper and my based on the material shape
        """
        if self.num_dim==1:
            self.dimension_name = 'x'
        elif self.num_dim==2:
            self.dimension_name = ['x','y']
        elif self.num_dim==3:
            self.dimension_name = ['x','y','z']

        #dimension = np.zeros([self.num_dim])
        if '_initialized' in self.__dict__:
            pass
        else:
            self._initialized = 1
            self.dimension_lower = np.zeros([self.num_dim])
            self.dimension_upper = np.zeros([self.num_dim])
            self.dimension_points = np.zeros([self.num_dim])
            self.dimension_resolution = np.zeros([self.num_dim])
            self.dimension_steps = np.zeros([self.num_dim])
            self.dimension_lower[0] = 0.0e-6                     # fill array with some initial values
            self.dimension_upper[0] = 100e-6
            self.dimension_resolution[0] = 60
            self.__t_final = 1.0
            self.__dxdt = 1.0
            self.__dydt = 1.0
        
        self.dimension_points[0] = np.floor(self.dimension_resolution[0]*(self.dimension_upper[0] - self.dimension_lower[0])/self.source_lambda)
        self.dimension_steps[0] = (self.dimension_upper[0]-self.dimension_lower[0])/self.dimension_points[0]

        if self.num_dim>=2:
            if '_initialized' in self.__dict__:
                pass
            else:
                dimension_lower[1] = 0.0e-6
                dimension_upper[1] = 1.0e-6                     # notice that for multilayer this is value will be over-written
                dimension_resolution[1] = 20
            
            if shape=='multilayer':
                self.aux_N_layers = np.floor(np.sum(self.aux_num_layers[:])/self.aux_n_layers) + 1
                self.dim_upper[1] = self.aux_N_layers*np.sum(self.aux_layers_thickness[:])
                self.dimension_tlp = np.sum(self.aux_layers_thickness[:])
                self.dimension_mlp = np.floor(self.tlp/1e-9)
                self.dimension_points[1] = np.floor(self.dimension_resolution[1]*(self.dimension_upper[1] - self.dimension_lower[1])/1e-9)
            else:
                self.dimension_points[1] = np.floor(self.dimension_resolution[1]*(self.dimension_upper[1] - self.dimension_lower[1])/self.source_lambda)

            self.dimension_steps[1] = (self.dimension_upper[1] - self.dimension_lower[1])/self.dimension_points[1]
            if dtcfl==1:
                self.__ddt = self.solver_cfl/(self.co*np.sqrt(1.0/(self.dimension_steps[0]**2) + 1.0/(self.dimension_steps[1]**2)))
                self.solver_dt = self.__ddt
        elif self.num_dim==1:
            if dtcfl==1:
                self.__ddt = 0.90/(self.co*np.sqrt(1.0/(self.dimension_steps[0]**2)))
                self.solver_dt = self.__ddt
                
        return self

    def set_material(self,shape=aux_shape,n_layers=_nlayers,N_layers=_Nlayers,rip=_is_rip):
        """
        set_material(shape=aux_shape,vx=0,vy=0,x_offset=10e-6,y_offset=0,sigma=10e-6,n_layers=2,N_layers=5)

        Initializes the basic material parameters, user most then specify the values require for each simulation.

        :Implemented mappings (shape = ...)
        
        ..gaussian1dx:  stationary and moving gaussian shape for eps and mu
        ..homogeneous:  homogeneous refractive index in eps and mu
        ..interface:    simple interface (jump) acroos the 2d domain
        ..interfacex:   simple interface (jump) 1D in x-direction
        ..interfacey:   ibid in y-direction
        ..multilayer:   2D multilayers in x or y direction.
        
        : Perturbation velocity vx = ..., vy = ...
        : Offsets, x_offset = ..., y_offset = ...
        : Perturbation spred (sigma) sigma = ...
        : Multilayer struture, n_layers = ... N_layers = ... (total number of sets)

        """

        # create a material base matrix
        if self.aux_tensor_kind=='isotropic':
            # background configuration
            self.aux_base[0,0] = 1.
            self.aux_base[1,1] = 1.
            self.aux_base[2,2] = 1.
        elif self.aux_tensor_kind=='anisotropic':
           # background configuration
            self.aux_base[0,0:1] = 1.
            self.aux_base[1,0:1] = 1.
            self.aux_base[2,2] = 1.
        elif self.aux_tensor_kind=='bianisotropic':
            # background configuration
            self.aux_base[:,:] = 1.
        
        # set the modifiers to the refractive index
        # if interface declare position
        if shape=='homogeneous':
            pass
        elif shape=='xinterface':
            self.aux_material = self.aux_base+1
            self.aux_xinterface = np.zeros([self.num_aux,self.num_aux])
            if self.aux_tensor_kind=='isotropic':
                np.fill_diagonal(self.aux_xinterface,self.dimension_upper[0]/2)
            elif self.aux_tensor_kind=='anisotropic' and self.num_dim>=2:
                np.fill_diagonal(self.aux_xinterface,self.dimension_upper[0]/2)
                self.aux_xinterface[self.num_aux-2,self.num_aux-2] = 1.0
            elif self.aux_tensor_kind=='bianisotropic':                
                self.aux_xinterface.fill(1.0) 
        elif shape=='yinterface':
            self.aux_material = self.aux_base+1
            self.aux_yinterface = np.zeros([self.num_aux,self.num_aux])
            if self.aux_tensor_kind=='isotropic':
                np.fill_diagonal(self.aux_yinterface,self.dimension_upper[1]/2)
            elif self.aux_tensor_kind=='anisotropic' and self.num_dim>=2:
                np.fill_diagonal(self.aux_yinterface,self.dimension_upper[1]/2)
                self.aux_yinterface[self.num_aux-2,self.num_aux-2] = 1.0
            elif self.aux_tensor_kind=='bianisotropic':                
                self.aux_yinterface.fill(1.0) 
        elif shape=='interface':
            self.aux_material = np.append(self.aux_base+1,self.aux_base+1,0).reshape(2,3,3)
            self.aux_interface = np.zeros([self.num_dim,self.num_aux,self.num_aux])
            if self.aux_tensor_kind=='isotropic':
                for j in range(0,self.num_dim):
                    np.fill_diagonal(self.aux_interface[j],self.dimension_upper[j]/2)
            elif self.aux_tensor_kind=='anisotropic' and self.num_dim>=2:
                for j in range(0,self.num_dim):
                    np.fill_diagonal(self.aux_interface[j],self.dimension_upper[j]/2)
                self.aux_interface[:,self.num_aux-2,self.num_aux-2] = 1.0
            elif self.aux_tensor_kind=='bianisotropic':                
                self.aux_interface.fill(1.0) 
        elif shape=='gaussian1dx' or shape=='gaussian':
            self.material = self.aux_base+1
            self.aux_gaussian_sigma =  self.aux_gaussian_offset = np.zeros([self.num_dim,self.num_aux,self.num_aux])
            if self.aux_tensor_kind=='isotropic':
                for j in range(0,self.num_dim):
                    np.fill_diagonal(self.aux_gaussian_sigma[j],1)
                    np.fill_diagonal(self.aux_gaussian_offset[j],1)
            elif self.aux_tensor_kind=='anisotropic' and self.num_dim>=2:
                for j in range(0,self.num_dim):
                    np.fill_diagonal(self.aux_gaussian_sigma[j],1)
                    np.fill_diagonal(self.aux_gaussian_offset[j],1)
                self.aux_gaussian_sigma[:,0:self.num_aux-1,0:self.num_aux-2] = self.aux_gaussian_offset[:,0:self.num_aux-2,0:self.num_aux-2] = 1.0
            elif self.aux_tensor_kind=='bianisotropic':
                self.aux_gaussian_sigma.fill(1.0)
                self.aux_gaussian_offset.fill(1.0)
        elif shape=='jump':
            pass
        elif shape=='multilayer':
            self.aux_n_layers = n_layers
            self.aux_material = np.zeros([n_layers,self.num_aux,self.num_aux])
            self.aux_layers_thickness = np.ones([n_layers])
            self.aux_num_layers = np.ones([n_layers])
            if self.aux_tensor_kind=='isotropic':
                np.fill_diagonal(self.aux_material,1.0)
            elif self.aux_tensor_kind=='anisotropic' and self.num_dim>=2:
                np.fill_diagonal(self.aux_material,1.0)
                self.aux_material[self.num_aux-2,self.num_aux-2] = 1.0
            elif self.aux_tensor_kind=='bianisotropic':                
                self.aux_material.fill(1.0) 
            # if self.aux_tensor_kind=='old_method':
            #     self.aux_n_layers = n_layers
            #     self.aux_layers = np.zeros([n_layers,7]) # _layer:  eps mu N t chi2e chi2m chi3e chi3m
            #     self.aux_layers[0,0] = 1.5
            #     self.aux_layers[0,1] = 1.5
            #     self.aux_layers[0,2] = 10
            #     self.aux_layers[0,3] = 15e-9
            #     self.aux_layers[1,0] = 2.5
            #     self.aux_layers[1,1] = 2.5
            #     self.aux_layers[1,2] = self.aux_layers[0,2] - 1
            #     self.aux_layers[1,3] = 50e-9
            #     self.aux_N_layers = N_layers
            # elif self.aux_tensor_kind=='isotropic':
                
        if rip:
            self.rip = self._is_rip
            self.aux_rip_velocity = np.zeros([self.num_dim,self.num_aux,self.num_aux])
            if self.aux_tensor_kind=='isotropic':
                for j in range(0,self.num_dim):
                    np.fill_diagonal(self.aux_rip_velocity[j],1.0)
            elif self.aux_tensor_kind=='anisotropic':
                for j in range(0,self.num_dim):
                    np.fill_diagonal(self.aux_rip_velocity[j],1.0)
                self.aux_rip_velocity[:,0:self.num_aux-2,0:self.num_aux-2] = 1.0

        return self

    def set_vacuum(self,vc_config=vacuum_config,mode_config=mode):

        self.vacuum = np.zeros(self.num_aux)
        if vc_config=='real':
            self.eo = 8.854187817e-12            # vacuum permittivity   - [F/m]
            self.mo = 4e-7*np.pi                 # vacuum peremeability  - [V.s/A.m]
            self.co = 1/np.sqrt(self.eo*self.mo)           # vacuum speed of light - [m/s]
            self.zo = np.sqrt(self.eo/self.mo)
        else:
            self.eo = 1            # vacuum permittivity   - [F/m]
            self.mo = 1                 # vacuum peremeability  - [V.s/A.m]
            self.co = 1/np.sqrt(self.eo*self.mo)           # vacuum speed of light - [m/s]
            self.zo = np.sqrt(self.eo/self.mo)
        
        if mode_config=='TE' and self.num_dim>=2:
            self.vacuum[0] = self.mo
            self.vacuum[1] = self.mo
            self.vacuum[2]  = self.eo
        elif mode_config=='TM' and self.num_dim>=2:
            self.vacuum[0] = self.eo
            self.vacuum[1] = self.eo
            self.vacuum[2]  = self.mo

        if self.num_dim==1:
            self.vacuum[0] = self.eo
            self.vacuum[1] = self.mo

        return self

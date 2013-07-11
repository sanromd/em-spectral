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

    num_frames  = 30
    t_final     = 1
    solver_cfl  = 0.9
    # ....... dimensions .............................................
    num_dim     = 2
    dtcfl       = 1
    # ....... aux config .............................................
    num_aux     = 3
    mode        = 'TE'
    vacuum_config = 'real'
    aux_shape   = 'homogeneous'
    aux_tensor_kind = 'isotropic'
    aux_base    = np.zeros([num_aux,num_aux])
    _is_rip     = 1
    _nlayers    = 2
    _Nlayers    = 5
    # ........ Boundary settings .................
    
    solver_tolerance = 1e-6

    # parameters needed for pml calculation
    pml_points = 8
    pml_norder = 3
    pml_Ro = 1.0e-6

    # ........ excitation - initial conditoons .......................
    source_lambda  = 1e-6 
    source_type  = 'off'
    dtcfl = 1


    def __init__(self, **kargs):
        cls_dict = self.__class__.__dict__
        for key in kargs:
            assert key in cls_dict
        self.__dict__.update(kargs)
        self._initialized = False
        self.set_vacuum(vc_config=self.vacuum_config)
        self.set_material(shape=self.aux_shape)
        self.set_source()
        self.set_dims(shape=self.aux_shape)
        self._initialized = True

    def set_source(self):
        if self._initialized:
            pass
        else:
            self.source_width      = 10.0 # width in the y-direction 
            self.source_t_offset   = 0.0
            self.source_offset     = [0.0, 0.0]      # offset in the y-direction
            self.source_time_width = 1.0*self.source_lambda   
            if self.num_dim==2:
                self.source_amp    = np.zeros([3])
                self.source_amp[0] = 0.0                           # Amplitude of q1, q2 and q3
                self.source_amp[1] = 1.0
                self.source_amp[2] = 1.0
            elif self.num_dim==1:
                self.source_amp    = np.zeros([2])
                self.source_amp[0] = 0.0                           # Amplitude of q1, q2 and q3
                self.source_amp[1] = 1.0

        self.source_k     = 2.0*np.pi/self.source_lambda       # wave vector        
        self.source_omega = 2.0*np.pi*self.co/self.source_lambda
        return self    

    def set_dims(self,shape=aux_shape,dtcfl=1):
        """
        set_dims(shape=aux_shape)
        calculates y_upper and my based on the material shape
        """
        if self._initialized:
            pass
        else:
            if self.num_dim==1:
                self.dimension_name = ['x']
                self.bc_lower       = ['scattering']
                self.bc_upper       = ['none']
                self.aux_bc_lower   = ['scattering']
                self.aux_bc_upper   = ['pml']
            elif self.num_dim==2:
                self.dimension_name = ['x','y']
                self.bc_lower       = ['scattering', 'none']
                self.bc_upper       = ['none', 'none']
                self.aux_bc_lower   = ['scattering','metallic']
                self.aux_bc_upper   = ['pml','metallic']
            elif self.num_dim==3:
                self.dimension_name = ['x','y','z']
                self.bc_lower       = None
                self.bc_upper       = None
                self.aux_bc_lower   = None
                self.aux_bc_upper   = None

        if self._initialized:
            pass
        else:
            print self.num_dim
            self.dimension_lower      = np.zeros([self.num_dim])
            self.dimension_upper      = np.zeros([self.num_dim])
            self.dimension_num_points = np.zeros([self.num_dim])
            self.dimension_resolution = np.zeros([self.num_dim])
            self.dimension_step_size  = np.zeros([self.num_dim])
            self.dimension_lower[0]   = 0.0e-6                 
            self.dimension_upper[0]   = 100e-6
            self.dimension_resolution[0] = 10
            self.__t_final  = 1.0
            self.__dxdt     = 1.0
            self.__dydt     = 1.0
        
        self.dimension_num_points[0]    = np.floor(self.dimension_resolution[0]*(self.dimension_upper[0] - self.dimension_lower[0])/self.source_lambda)
        self.dimension_step_size[0]     = (self.dimension_upper[0]-self.dimension_lower[0])/self.dimension_num_points[0]

        if self.num_dim>=2:
            if self._initialized:
                pass
            else:
                self.dimension_lower[1]      = 0.0e-6
                self.dimension_upper[1]      = 10.0e-6                     # notice that for multilayer this is value will be over-written
                self.dimension_resolution[1] = 10
            
            if shape=='multilayer':
                self.aux_N_layers   = np.floor(np.sum(self.aux_num_layers[:])/self.aux_n_layers) + 1
                self.dim_upper[1]   = self.aux_N_layers*np.sum(self.aux_layers_thickness[:])
                self.dimension_tlp  = np.sum(self.aux_layers_thickness[:])
                self.dimension_mlp  = np.floor(self.tlp/1e-9)
                self.dimension_num_points[1] = np.floor(self.dimension_resolution[1]*(self.dimension_upper[1] - self.dimension_lower[1])/1e-9)
            else:
                self.dimension_num_points[1] = np.floor(self.dimension_resolution[1]*(self.dimension_upper[1] - self.dimension_lower[1])/self.source_lambda)

            self.dimension_step_size[1] = (self.dimension_upper[1] - self.dimension_lower[1])/self.dimension_num_points[1]
            if dtcfl==1:
                self._dt_cfl    = self.solver_cfl/(self.co*np.sqrt(1.0/(self.dimension_step_size[0]**2) + 1.0/(self.dimension_step_size[1]**2)))
                self.solver_dt  = self._dt_cfl
            else:
                self._dt_cfl    = 0.0
                self.solver_dt  = 0.1
        elif self.num_dim==1:
            if dtcfl==1:
                self.__dt_cfl  = 0.90/(self.co*np.sqrt(1.0/(self.dimension_step_size[0]**2)))
                self.solver_dt = self.__dt_cfl
            else:
                self._dt_cfl    = 0.0
                self.solver_dt  = 0.1

        return self

    def set_material(self,shape=aux_shape,n_layers=_nlayers,N_layers=_Nlayers,rip=_is_rip):
        r"""
        set_material(shape=aux_shape,rip=_is_rip)

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
            self.aux_base[2,2]   = 1.
        elif self.aux_tensor_kind=='bianisotropic':
            # background configuration
            self.aux_base[:,:] = 1.
        
            self.aux_chi2 = 0.
            self.aux_chi3 = 0.
        # set the modifiers to the refractive index
        # if interface declare position
        if shape=='homogeneous':
            pass
        elif shape=='xinterface':
            self.material_type = 1
            self.aux_material   = self.aux_base+1
            self.aux_interface = np.zeros([self.num_aux,self.num_aux])
            if self.aux_tensor_kind=='isotropic':
                np.fill_diagonal(self.aux_interface,self.dimension_upper[0]/2)
            elif self.aux_tensor_kind=='anisotropic' and self.num_dim>=2:
                np.fill_diagonal(self.aux_interface,self.dimension_upper[0]/2)
                self.aux_interface[self.num_aux-2,self.num_aux-2] = 1.0
            elif self.aux_tensor_kind=='bianisotropic':                
                self.aux_interface.fill(1.0) 
        elif shape=='yinterface':
            self.material_type = 2
            self.aux_material   = self.aux_base+1
            self.aux_interface = np.zeros([self.num_aux,self.num_aux])
            if self.aux_tensor_kind=='isotropic':
                np.fill_diagonal(self.aux_interface,self.dimension_upper[1]/2)
            elif self.aux_tensor_kind=='anisotropic' and self.num_dim>=2:
                np.fill_diagonal(self.aux_interface,self.dimension_upper[1]/2)
                self.aux_interface[self.num_aux-2,self.num_aux-2] = 1.0
            elif self.aux_tensor_kind=='bianisotropic':                
                self.aux_interface.fill(1.0) 
        elif shape=='interface':
            self.material_type = 3
            self.aux_material   = np.append(self.aux_base+1,self.aux_base+1,0).reshape(2,3,3)
            self.aux_interface  = np.zeros([self.num_dim,self.num_aux,self.num_aux])
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
            if shape=='gaussian1dx':
                self.material_type = 4
            else:
                self.material_type = 5

            self.aux_base  = 1.5*self.aux_base
            self.aux_delta = 0.1*self.aux_base
            self.aux_sigma = self.aux_offset = np.zeros([self.num_dim,self.num_aux,self.num_aux])
            if self.aux_tensor_kind=='isotropic':
                for j in range(0,self.num_dim):
                    np.fill_diagonal(self.aux_sigma[j],1)
                    np.fill_diagonal(self.aux_offset[j],1)
            elif self.aux_tensor_kind=='anisotropic' and self.num_dim>=2:
                for j in range(0,self.num_dim):
                    np.fill_diagonal(self.aux_sigma[j],1)
                    np.fill_diagonal(self.aux_offset[j],1)
                self.aux_sigma[:,0:self.num_aux-1,0:self.num_aux-2] = self.aux_offset[:,0:self.num_aux-2,0:self.num_aux-2] = 1.0
            elif self.aux_tensor_kind=='bianisotropic':
                self.aux_sigma.fill(1.0)
                self.aux_offset.fill(1.0)
        elif shape=='jump':
            self.material_type = 6
            pass
        elif shape=='multilayer':
            self.material_type = 7
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

        if rip:
            self.rip = self._is_rip
            self.aux_velocity = np.zeros([self.num_dim,self.num_aux,self.num_aux])
            if self.aux_tensor_kind=='isotropic':
                for j in range(0,self.num_dim):
                    np.fill_diagonal(self.aux_velocity[j],1.0)
            elif self.aux_tensor_kind=='anisotropic':
                for j in range(0,self.num_dim):
                    np.fill_diagonal(self.aux_velocity[j],1.0)
                self.aux_velocity[:,0:self.num_aux-2,0:self.num_aux-2] = 1.0

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
        
        if mode_config=='TE' and self.num_dim==2:
            self.vacuum[0] = self.mo
            self.vacuum[1] = self.mo
            self.vacuum[2] = self.eo
        elif mode_config=='TM' and self.num_dim==2:
            self.vacuum[0] = self.eo
            self.vacuum[1] = self.eo
            self.vacuum[2] = self.mo

        if self.num_dim==1:
            self.vacuum[0] = self.eo
            self.vacuum[1] = self.mo

        return self

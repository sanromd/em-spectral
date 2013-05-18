#!/usr/bin/python
import numpy as np
import simulation_parameters as parameters
class material:
#    ------- vacuum
    def vacuum(self):
        if self.vacuum_config=='real':
            self.eo = 8.854187817e-12            # vacuum permittivity   - [F/m]
            self.mo = 4e-7*np.pi                 # vacuum peremeability  - [V.s/A.m]
            self.co = 1/np.sqrt(self.eo*self.mo)           # vacuum speed of light - [m/s]
            self.zo = np.sqrt(self.eo/self.mo)
        elif self.vacuum_config=='one':
            self.eo = 1            # vacuum permittivity   - [F/m]
            self.mo = 1                 # vacuum peremeability  - [V.s/A.m]
            self.co = 1/np.sqrt(self.eo*self.mo)           # vacuum speed of light - [m/s]
            self.zo = np.sqrt(self.eo/self.mo)
        
        return self

    def __init__(self):
        self.vacuum_config = parameters.vacuum_config
        self.bkg_er = parameters.bkg_er
        self.bkg_mr = parameters.bkg_mr
        self.vacuum()

    # # ------- material
    # def build_material(self,shape=parameters.shape):
    #     """
    #         buiild_material(shape=parameters.shape)

    #         material definition: 
            
    #         ..: homogeneous
    #         ..: interface 
    #         ..: rip (moving perturbation) 
    #         ..: multilayered
    #     """

    #     # build vacuum properties
    #     vacuum(self,self.vacuum_config)

    #     # build background
    #     self.bkg_n  = self.np.sqrt(self.bkg_er*self.bkg_mr)
    #     self.bkg_e  = self.eo*self.bkg_er
    #     self.bkg_m  = self.mo*self.bkg_mr
        
    #     # get material shape
    #     self.mat_shape = shape   

    #     # if interface declare position
    #     material.x_change = parameters.x_upper/2

    #     # set moving refractive index parameters
    #     material.rip_vx_e    = 0.0*co    # replace here the value of x
    #     rip_vx_m    = rip_vx_e
    #     rip_vy_e    = 0.0*co
    #     rip_vy_m    = rip_vy_e

    #     rip_xoff_e  = 10e-6
    #     rip_xoff_m  = rip_xoff_e
    #     rip_yoff_e  = rip_xoff_e
    #     rip_yoff_m  = rip_xoff_e

    #     rip_xsig_e  = 10.0e-6
    #     rip_xsig_m  = rip_xsig_e
    #     rip_ysig_e  = .9*y_upper/2
    #     rip_ysig_m  = rip_ysig_e
    #     s_x_e       = rip_xsig_e**2
    #     s_x_m       = rip_xsig_m**2
    #     s_y_e       = rip_ysig_e**2
    #     s_y_m       = rip_ysig_m**2

    #     prip        = 0.1
    #     deltan      = prip*(bkg_n) # assumes epsilon = mu
    #     d_e         = deltan #*(2.0*1.5+deltan)
    #     d_m         = deltan #*(2.0*1.5+deltan)

    #     # set multilayer parameters

    #     # multilayered definition
    #     n_layers = 2
    #     layers = np.zeros([n_layers,7]) # _layer:  eps mu N t chi2e chi2m chi3e chi3m
    #     layers[0,0] = 1.5
    #     layers[0,1] = 1.5
    #     layers[0,2] = 10
    #     layers[0,3] = 15e-9
    #     layers[1,0] = 2.5
    #     layers[1,1] = 2.5
    #     layers[1,2] = layers[0,2] - 1
    #     layers[1,3] = 50e-9
    #     N_layers = 5
    #     if mat_shape=='multilayer':
    #         y_upper = N_layers*np.sum(layers[:,3])+layers[0,3]
    #         tlp = np.sum(layers[:,3])
    #         mlp = np.floor(tlp/1e-9)

    # def etar(da,ddx,ddy,t=0):
    #     """
    #     eta = etar(num_aux,xi,xf,yi,yf,ddx,ddy)

    #     Sets the auxiliary arrays for permittivity and permeability.

    #     Implemented mappings

    #     ..gaussian1dx:  stationary and moving gaussian shape for eps and mu
    #     ..homogeneous:  homogeneous refractive index in eps and mu
    #     ..interface:    simple interface (jump) acroos the 2d domain
    #     ..interfacex:   simple interface (jump) 1D in x-direction
    #     ..interfacey:   ibid in y-direction
    #     ..multilayer:   2D multilayers in x or y direction.


    #     y,x are the point coordinates of the grid.
    #     t is the time coordinate

    #     on output aux holds:

    #                                EM equivalent

    #          idim   curvilinear  |   TE      TM
    #          0:     eta_1        |   mu1     eps1
    #          1:     eta_2        |   mu2     eps2
    #          2:     eta_3        |   eps3    mu3

    #     """
    #     nx, ny = da.getSizes()
    #     (xi, xf), (yi, yf) = da.getRanges()
    #     X = np.linspace(xi,xf,xf-xi)*ddx
    #     Y = np.linspace(yi,yf,yf-yi)*ddy
    #     y,x = np.meshgrid(Y,X)
    #     eta = np.empty( [3,len(X),len(Y)], order='F')

    #     if mat_shape=='gaussian1dx':
    #         u_x_e = x - rip_vx_e*t - rip_xoff_e
    #         u_x_m = x - rip_vx_m*t - rip_xoff_m
    #         u_y_e = y - rip_vy_e*t - rip_yoff_e
    #         u_y_m = y - rip_vy_m*t - rip_yoff_m

    #         u_e = (u_x_e/rip_xsig_e)**2 + (u_y_e/rip_ysig_e)**2
    #         u_m = (u_x_m/rip_xsig_m)**2 + (u_y_m/rip_ysig_m)**2

    #         eta[0,:,:] = d_e*np.exp(-u_e) + bkg_er
    #         eta[1,:,:] = d_e*np.exp(-u_e) + bkg_er
    #         eta[2,:,:] = d_m*np.exp(-u_m) + bkg_mr
    #     elif mat_shape=='homogeneous':
    #         eta[0,:,:] = bkg_er
    #         eta[1,:,:] = bkg_er
    #         eta[2,:,:] = bkg_mr
    #     elif mat_shape=='vacuum':
    #         eta[0,:,:] = 1.0
    #         eta[1,:,:] = 1.0
    #         eta[2,:,:] = 1.0
    #     elif mat_shape=='interfacex':
    #         eta[0,:,:] = 1*(x<x_change) + 4*(x>=x_change)
    #         eta[1,:,:] = 1*(x<x_change) + 4*(x>=x_change)
    #         eta[2,:,:] = 1*(x<x_change) + 4*(x>=x_change)
    #     elif mat_shape=='interfacey':
    #         yy = y_upper-y_lower
    #         eta[0,:,:] = 1*(y<yy/2) + 4*(x>=yy/2)
    #         eta[1,:,:] = 1*(y<yy/2) + 4*(x>=yy/2)
    #         eta[2,:,:] = 1*(y<yy/2) + 4*(x>=yy/2)
    #     elif mat_shape=='multilayer':
    #         for n in range(0,N_layers):
    #             yi = n*tlp
    #             for m in range(0,n_layers):
    #                 if m==0:
    #                     eta[0,:,:] = layers[m,0]*(yi<y)*(y<=yi+layers[m,3])
    #                     eta[1,:,:] = layers[m,1]*(yi<y)*(y<=yi+layers[m,3])
    #                     eta[2,:,:] = layers[m,1]*(yi<y)*(y<=yi+layers[m,3])
    #                 else:
    #                     eta[0,:,:] = layers[m,0]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])
    #                     eta[1,:,:] = layers[m,1]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])
    #                     eta[2,:,:] = layers[m,1]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])


    #         eta[0,:,:] = layers[0,0]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,3])
    #         eta[1,:,:] = layers[0,1]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,3])
    #         eta[2,:,:] = layers[0,1]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,3])

    #     return eta

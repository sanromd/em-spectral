#!/usr/bin/env python
# encoding: utf-8
#
# One dimensional solver for Maxwell's equation via semi-spectral methods
#--------- non-linear em --------------------------

import numpy as np

# -------- GLOBAL SCALAR DEFINITIONS -----------------------------
# ======== all definitions are in m,s,g unit system.
x_lower = 0.
x_upper = 10e-6                 # lenght [m]
# ........ material properties ...................................

# vacuum
eo = 8.854187817e-12            # vacuum permittivity   - [F/m]
mo = 4e-7*np.pi                 # vacuum peremeability  - [V.s/A.m]
co = 1/np.sqrt(eo*mo)           # vacuum speed of light - [m/s]
zo = np.sqrt(mo/eo)

# material
mat_shape = 'homogeneous'           # material definition: homogeneous, interface, rip (moving perturbation), multilayered

# background refractive index 
bkg_er = 1.5
bkg_mr = 1.5
bkg_n  = np.sqrt(bkg_er*bkg_mr)
bkg_e  = eo*bkg_er
bkg_m  = mo*bkg_mr

# if interface declare position
x_change = x_upper/2

# set moving refractive index parameters
rip_vx_e    = 0.0*co    # replace here the value of x
rip_vx_m    = rip_vx_e

rip_xoff_e  = 10e-6
rip_xoff_m  = rip_xoff_e

rip_xsig_e  = 10.0e-6
rip_xsig_m  = rip_xsig_e
s_x_e       = rip_xsig_e**2
s_x_m       = rip_xsig_m**2

prip        = 0.1
deltan      = prip*(bkg_n) # assumes epsilon = mu
d_e         = deltan #*(2.0*1.5+deltan)
d_m         = deltan #*(2.0*1.5+deltan)

# set multilayer parameters

# multilayered definition
n_layers = 2
layers = np.zeros([n_layers,7]) # _layer:  eps mu N t chi2e chi2m chi3e chi3m
layers[0,0] = 1.5
layers[0,1] = 1.5
layers[0,2] = 10
layers[0,3] = 15e-9
layers[1,0] = 2.5
layers[1,1] = 2.5
layers[1,2] = layers[0,2] - 1
layers[1,3] = 50e-9
N_layers = 5
if mat_shape=='multilayer':
    x_upper = N_layers*np.sum(layers[:,3])+layers[0,3]
    tlp = np.sum(layers[:,3])
    mlp = np.floor(tlp/1e-9)

# set non-linear parameters of the material
chi2_e      = 0.0
chi3_e      = 0.0
chi2_m      = 0.0
chi3_m      = 0.0

# ........ excitation - initial conditoons .......................
ex_type  = 'off'
alambda  = 1e-6             # wavelength
ex_t_sig = 1.0*alambda          # width in time (pulse only)
ex_x_sig = 1.0*alambda          # width in the x-direction (pulse)
ex_toff  = 0.0                  # offset in time
ex_xoff  = 0.0                  # offset in the x-direction
omega    = 2.0*np.pi*co/alambda # frequency
k        = 2.0*np.pi/alambda
amp_Ey   = 1.
amp_Hz   = 1.

	# ........ pre-calculations for wave propagation .................
v_r = 1./bkg_n
v = co*v_r
ex_vx = v
ex_kx = k

# Grid - mesh settings
# Adapt the domain to go from -2pi to 2pi
sd = (x_upper-x_lower)
sf = 4.*np.pi/sd
sf_off = -sd/2. + np.abs(x_lower)

Lx = np.floor((x_upper-x_lower)/(2*np.pi))
if mat_shape=='multilayer':
    mx = np.floor((x_upper-x_lower)/1e-9)
else:
    mx = np.floor(20*4*np.pi/alambda)

ddx = (x_upper-x_lower)/mx
ddt = 0.90/(co*np.sqrt(1.0/(ddx**2)))
max_steps = 250000
t_final = (x_upper-x_lower)/v

tolerance=0.1^6;    # tolerance for fixed point iterations

# -------- GLOBAL FUNCTION DEFINITIONS --------------

# refractive index map definition function 
def etar(t,x):
    pass

def qinit(x):
    """
    Initial conditions in simulation grid for electromagnetic components q
    """
    
    if ex_type=='off':
        dd = x_upper-x_lower
        q[0,:] = 0.
        q[1,:] = np.exp(-(x-dd/2)**2/(dd/10)**2)
    else:
        q[0,:] = 0.0
        q[1,:] = 0.0
	
	return q

# Grid points & initial time
x = np.linspace(x_lower,x_upper,mx+1)
t = 0

# Calculate initial field
fields = qinit(x)
electric = field[0,:]
magnetic = field[1,:]

# Define temporary variables
electricold = electric
magneticold = magnetic

# Calculate refractive index at t0
aux = etar(t,x)
eps = aux[0,:]
mu = aux[1,:]

dt = ddt # Timestep size
 
# initialise variables

k = 1i*[0:N/2-1 0 -N/2+1:-1]/Lx; # wave vector
 
% Setting up Plot
tmax = 150; 
nplots = round(tmax/dt);
electricdata = [electric; zeros(nplots,N)]; tdata = t;
magneticdata = [magnetic; zeros(nplots,N)]; 
 
for i = 1:nplots
    t=t+dt
    muold=mu;
    electricold=electric;
    magneticold=magnetic;
    mu=A*exp((x-v*t)/sigma^2)+B;
    difference=1;
    while (difference>tolerance)
        electrictemp = electric;
        magnetictemp = magnetic;
        electricmean= 0.5*(electric + electricold);
        magneticmean = 0.5*(magnetic+magneticold);
        electric = (muold*electricold + dt*np.fft.ifft(k*np.fft.ifft(magneticmean)) + (2*chi2*electricmean + 3*chi3*electricmean**2)*electricold)/(mu+2*chi2*electricmean+3*chi3*electricmean.^2);
        magnetic = (muold.*magneticold + dt*ifft(k.*fft(electricmean)) +  (2*chi2*magneticmean+3*chi3*magneticmean**2).*magneticold)/(mu+2*chi2*magneticmean+3*chi3*magneticmean.^2);        
        difference = max(abs(electrictemp-electric))+max(abs(magnetictemp-magnetic));

    electricdata(i+1,:) = real(electric); # Records data
    magneticdata(i+1,:) = real(magnetic); # Records data
    tdata = [tdata; t];
end
 
%Plot using mesh
figure(1)
subplot(2,1,1);
mesh(x,tdata,electricdata), grid on, %axis([-1 2*pi 0 tmax -1 1]),
view(-60,55), xlabel x, ylabel t, zlabel electric, colorbar
subplot(2,1,2);
mesh(x,tdata,magneticdata), grid on, %axis([-1 2*pi 0 tmax -1 1]),
view(-60,55), xlabel x, ylabel t, zlabel magnetic, colorbar

figure(2)
subplot(2,1,1);
pcolor(x,tdata,electricdata), shading interp, %axis([-1 2*pi 0 tmax -1 1]),
xlabel x, ylabel t, title electric, colorbar
subplot(2,1,2);
pcolor(x,tdata,magneticdata), shading interp, %axis([-1 2*pi 0 tmax -1 1]),
xlabel x, ylabel t, title magnetic, colorbar

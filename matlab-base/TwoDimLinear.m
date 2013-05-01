% A program to solve a 2D reduction of the 3D Maxwell's equations 
% using the implicit midpoint rule
% E_t= \nabla\times H
% H_t=-\nabla\times E
% the reduction is 
% E=(E^x(x,y),E^y(x,y),0)
% H=(0,0,H^z(x,y))
clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')
 
% set up grid
tic
Lx = 4; % period 2*pi*L
Ly = 4; % period 2*pi*L
Nx = 64; % number of harmonics
Ny = 64; % number of harmonics
Nt = 100; % number of time slices
dt = 1.0/Nt; % time step
plotgap=1; % time steps between plots 
tolerance=0.1^6; % tolerance for fixed point iterations

% initialise variables
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)'*Lx; % x coordinate
kx = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx; % wave vector
y = (2*pi/Ny)*(-Ny/2:Ny/2 -1)'*Ly; % y coordinate
ky = 1i*[0:Ny/2-1 0 -Ny/2+1:-1]'/Ly; % wave vector
[xx,yy]=meshgrid(x,y);
[kxm,kym]=meshgrid(kx,ky);
 
% initial conditions
Ex = 0.0*fftn(exp(-5*(yy.^2+xx.^2)));
Ey = 0.0*Ex;
Hz = fftn(exp(-5*(xx.^2+yy.^2)));

figure(1); clf;
subplot(3,1,1);
pcolor(xx,yy,real(ifftn(Hz))), grid on, shading interp, 
xlabel x, ylabel t, title Magnetic^z, colorbar
subplot(3,1,2);
pcolor(xx,yy,real(ifftn(Ex))), grid on, shading interp, 
xlabel x, ylabel t, title Electric^x, colorbar
subplot(3,1,3);
pcolor(xx,yy,real(ifftn(Ey))), grid on, shading interp, 
xlabel x, ylabel t, title Electric^y, colorbar
t=0; tdata(1)=t;
 
 
% solve pde and plot results
 
for n =2:Nt+1
    Exold=Ex; 
    Eyold=Ey; 
    Hzold=Hz;
    difference=1;
    while (difference>tolerance)
        Extemp=Ex;
        Eytemp=Ey;
        Hztemp=Hz;
        
        Exmean=0.5*(Ex+Exold);
        Eymean=0.5*(Ey+Eyold);
        Hzmean=0.5*(Hz+Hzold);
        
        Ex = Exold + dt*(kym.*Hzmean); 
        Ey = Eyold - dt*(kxm.*Hzmean); 
        Hz = Hzold - dt*(kxm.*Eymean-kym.*Exmean);   
        
        difference= max(abs(Extemp-Ex))+max(abs(Eytemp-Ey))+max(abs(Hztemp-Hz));
    end

    if (mod(n,plotgap)==0)
        figure(1); clf;
        subplot(3,1,1);
        pcolor(xx,yy,real(ifftn(Hz))), grid on, shading interp,     
        xlabel x, ylabel t, title Magnetic^z, colorbar
        subplot(3,1,2);
        pcolor(xx,yy,real(ifftn(Ex))), grid on, shading interp, 
        xlabel x, ylabel t, title Electric^x, colorbar
        subplot(3,1,3);
        pcolor(xx,yy,real(ifftn(Ey))), grid on, shading interp, 
        xlabel x, ylabel t, title Electric^y, colorbar
     end
end
toc
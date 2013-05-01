% A program to solve the 3D Maxwell's equations using the implicit midpoint
% rule
% E_t= \nabla\times H
% H_t=-\nabla\times E
clear all; format compact; format short;
set(0,'defaultaxesfontsize',30,'defaultaxeslinewidth',.7,...
    'defaultlinelinewidth',6,'defaultpatchlinewidth',3.7,...
    'defaultaxesfontweight','bold')
 
% set up grid
tic
Lx = 4; % period 2*pi*L
Ly = 4; % period 2*pi*L
Lz = 4; % period 2*pi*L
Nx = 64; % number of harmonics
Ny = 64; % number of harmonics
Nz = 64; % number of harmonics
Nt = 100; % number of time slices
dt = 1.0/Nt; % time step
plotgap=1; % time steps between plots 
tolerance=0.1^6; % tolerance for fixed point iterations

% initialise variables
x = (2*pi/Nx)*(-Nx/2:Nx/2 -1)'*Lx; % x coordinate
kx = 1i*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx; % wave vector
y = (2*pi/Ny)*(-Ny/2:Ny/2 -1)'*Ly; % y coordinate
ky = 1i*[0:Ny/2-1 0 -Ny/2+1:-1]'/Ly; % wave vector
z = (2*pi/Nz)*(-Nz/2:Nz/2 -1)'*Lz; % y coordinate
kz = 1i*[0:Nz/2-1 0 -Nz/2+1:-1]'/Lz; % wave vector
[xx,yy,zz]=meshgrid(x,y,z);
[kxm,kym,kzm]=meshgrid(kx,ky,kz);
 
% initial conditions
Ex = fftn(exp(-5*(yy.^2+zz.^2)));
Ey = 0.0*Ex;
Ez = 0.0*Ex;
Hx = 0.0*Ex;
Hy = 0.0*Ex;
Hz = fftn(exp(-5*(xx.^2+yy.^2)));

figure(1); clf; UP = real(ifftn(Ex));
p1 = patch(isosurface(x,y,z,UP,.0025),...
    'FaceColor','yellow','EdgeColor','none');
p2 = patch(isocaps(x,y,z,UP,.0025),...
    'FaceColor','interp','EdgeColor','none');
isonormals(UP,p1); lighting phong;
xlabel('x'); ylabel('y'); zlabel('z');
axis equal; axis square; view(3); drawnow;
t=0; tdata(1)=t;
 
 
% solve pde and plot results
 
for n =2:Nt+1
    Exold=Ex; Eyold=Ey; Ezold=Ez;
    Hxold=Hx; Hyold=Hy; Hzold=Hz;
    difference=1;
    while (difference>tolerance)
        Extemp=Ex;Eytemp=Ey;Eztemp=Ez;
        Hxtemp=Hx;Hytemp=Hy;Hztemp=Hz;
        
        Exmean=0.5*(Ex+Exold);Eymean=0.5*(Ey+Eyold);Ezmean=0.5*(Ez+Ezold);
        Hxmean=0.5*(Hx+Hxold);Hymean=0.5*(Hy+Hyold);Hzmean=0.5*(Hz+Hzold);
        Ex = Exold + dt*(kym.*Hzmean-kzm.*Hymean); 
        Ey = Eyold + dt*(kzm.*Hxmean-kxm.*Hzmean); 
        Ez = Ezold + dt*(kxm.*Hymean-kym.*Hxmean);
        Hx = Hxold - dt*(kym.*Ezmean-kzm.*Eymean); 
        Hy = Hyold - dt*(kzm.*Exmean-kxm.*Ezmean); 
        Hz = Hzold - dt*(kxm.*Eymean-kym.*Exmean);
        
        
        difference= max(abs(Extemp-Ex))+max(abs(Eytemp-Ey))+max(abs(Eztemp-Ez))+...
                    max(abs(Hxtemp-Hy))+max(abs(Hytemp-Hy))+max(abs(Hztemp-Hz));
    end

    if (mod(n,plotgap)==0)
        figure(1); clf; UP = real(ifftn(Ex));
        p1 = patch(isosurface(x,y,z,UP,.0025),...
            'FaceColor','yellow','EdgeColor','none');
        p2 = patch(isocaps(x,y,z,UP,.0025),...
            'FaceColor','interp','EdgeColor','none');
        isonormals(UP,p1); lighting phong;
        xlabel('x'); ylabel('y'); zlabel('z');
        axis equal; axis square; view(3); drawnow;
     end
end
figure(4); clf; UP = real(ifftn(Ex));
p1 = patch(isosurface(x,y,z,UP,.0025),...
    'FaceColor','yellow','EdgeColor','none');
p2 = patch(isocaps(x,y,z,UP,.0025),...
    'FaceColor','interp','EdgeColor','none');
isonormals(UP,p1); lighting phong;
xlabel('x'); ylabel('y'); zlabel('z');
axis equal; axis square; view(3); drawnow;
toc
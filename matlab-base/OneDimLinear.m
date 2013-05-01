%Solving One Dimensional Maxwell's equations using the implicit midpoint
%rule
%magnetic_t= electric_x
%electric_x= magnetic_t
% This is a reduction of Maxwell's equations with:
% E=(0,0,e(x)) and
% M=(0,m(x),0)
%BC = (Periodic)
clear all; clc;
 
%Grid
N = 512;            % grid points
Lx=6;                % size of box in multiples of 2pi
x = (2*pi/N)*(-N/2:N/2 -1)*Lx; % x coordinate
tolerance=0.1^6;    % tolerance for fixed point iterations


% Initial conditions
t = 0;
electric = fft(10*exp(-6*x.^2).*sin(8*x));
magnetic = 0.0*x;
dt = .01; % Timestep size
 
% initialise variables

k = 1i*[0:N/2-1 0 -N/2+1:-1]/Lx; % wave vector
 
% Setting up Plot
tmax = 10; 
nplots = round(tmax/dt);
electricdata = [real(ifft(electric)); zeros(nplots,N)]; tdata = t;
magneticdata = [real(ifft(magnetic)); zeros(nplots,N)]; 
 
for i = 1:nplots
    t=t+dt
    electricold=electric;
    magneticold=magnetic;
    difference=1;
    while (difference>tolerance)
        electrictemp=electric;
        magnetictemp=magnetic;
        electricmean=0.5*(electric+electricold);
        magneticmean=0.5*(magnetic+magneticold);
        electric = electricold + dt*k.*magneticmean;
        magnetic = magneticold + dt*k.*electricmean;
        
        difference=max(abs(electrictemp-electric))+max(abs(magnetictemp-magnetic));
    end
    electricdata(i+1,:) = real(ifft(electric)); %Records data
    magneticdata(i+1,:) = real(ifft(magnetic)); %Records data
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
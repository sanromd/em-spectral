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
% Problem parameters
A=0.15;
B=0.6;
v=0.9;
sigma=5;
chi2=0.01;
chi3=0.00001;

% ---------- Initial conditions
t = 0;
electric = 10*exp(-6*x.^2).*sin(8*x);
magnetic = 0.0*x;%;
electricold=electric;
magneticold=magnetic;

mu=A*exp(-((x-v*t)/sigma).^2)+B;
dt = .01; % Timestep size
 
% ---------- Initialise variables
k = 1i*[0:N/2-1 0 -N/2+1:-1]/Lx; % wave vector

for i = 1:nplots
    t=t+dt
    muold=mu;
    electricold=electric;
    magneticold=magnetic;
    mu=A*exp((x-v*t)/sigma^2)+B;
    difference=1;
    while (difference>tolerance)
        electrictemp=electric;
        magnetictemp=magnetic;
        electricmean=0.5*(electric+electricold);
        magneticmean=0.5*(magnetic+magneticold);
        electric=(muold.*electricold + dt*ifft(k.*fft(magneticmean)) ...
            +(2*chi2*electricmean+3*chi3*electricmean.^2).*electricold)./(mu+2*chi2*electricmean+3*chi3*electricmean.^2);
        magnetic=(muold.*magneticold + dt*ifft(k.*fft(electricmean)) ...
            +(2*chi2*magneticmean+3*chi3*magneticmean.^2).*magneticold)./(mu+2*chi2*magneticmean+3*chi3*magneticmean.^2);
        
        difference=max(abs(electrictemp-electric))+max(abs(magnetictemp-magnetic));
    end
    electricdata(i+1,:) = real(electric); %Records data
    magneticdata(i+1,:) = real(magnetic); %Records data
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
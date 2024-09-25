%
%   Simple script examining stability and phase speed error in the upwind
%   solution approach to the advection equation.
%
%   Jim Powell, 9 Sept 2024
%
%       c   -   wave speed in advection equation
%       L   -   length of interval
%       J   -   number of sample points (in x)
%       dx  -   delta x =L/J
%       dt  -   discete time step delta t
%       phi -   Courant ration c*dt/dx
%
%   If we assume u_j^n = u( j dx, n dt) = A^n exp(i k*2*pi/L j*dx), the
%   upwind method gives us
%
%           A = 1 -phi + phi exp(-i k*2*pi/L dx)
%           lambda exp(-i omega dt ) = A = 1 - phi + phi (cos(k*2*pi/L dx) - i sin(k*2*pi/L dx) )
%
%   Stability and phase errors can be analyzed via lambda and omega.
%

% set up grids for k to use in plotting
J=100; L=10; dx=L/J; c=1;
% note actual wavenumber khat = 2*pi*k/L = 2*pi*k/J J/L  = 2*pi*k/J 1/dx, 
% where k is the mode number (up to J/2);
khat=2*pi/L*[1:J/2]

[1:J/2]

dx

phi=0.8;

figure
subplot(2,1,1)  % top subplot of two
lambda=abs(1-phi + phi*exp(-i*khat*dx));
plot(khat,lambda, khat, 1+0*khat,'--'), axis([0 max(khat) 0 1.25])
title(['Stability Plot for \phi = ' num2str(phi) ])
ylabel('\lambda')

% now the frequency
dt=phi/c*dx
omega=1/dt*atan2(phi*sin(khat*dx), 1 - phi + phi*cos(khat*dx) )
c_k=omega./khat;  % phase speed of each numerical wave
c_ratio= c_k/c;

subplot(2,1,2)  % bottom subplot of two
plot(khat,c_ratio, khat, 1+0*khat,'--'), axis([0 max(khat) 0 1.25])
title(['Phase speed ratio for \phi = ' num2str(phi) ])
xlabel('Wave number'), ylabel('c_k / c')


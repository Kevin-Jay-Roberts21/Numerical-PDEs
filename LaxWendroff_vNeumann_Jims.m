% 
% Script performing stability/phase error analysis for a third-order
% accurate approach to solving the advection equation, u_t + c u_x =0.  The
% third order scheme is
% 
%   u(j,n+1) = u(j,n) -1/2 phi ( u(j+1,n) - u(j-1,n) )
%                   +1/2 phi^2 (u(j+1,n) + u(j-1,n) -2 u(j,n) )
% 
% where:
%   phi = c delta t/delta x is the Courant number
% 
% Written by Jim Powell, 20 Sept. 2024
% 

% Just for definiteness, assume we are working on an interval of length 10
% with 100 subdivisions:
L=10; J=100; dx=L/J; c=1;

% define a vector of wave numbers
khat=2*pi*[1:J/2]/J;    % equivalent to 2*pi*k/L*dx

% define a vector of phi to test
nphi=50; phi=linspace(1/nphi,2,nphi); % we start at a small phi rather than zero 
                        % to avoid division by zero below

% define a 2D grid of k, phi for density plotting
[K,Phi]=meshgrid(khat,phi);

% calculate the (complex) von Neumann multiplier
A = 1 -.5*Phi.*(exp(i*K) -exp(-i*K)) ...
    +.5*Phi.^2.*(exp(i*K) +exp(-i*K) -2);
Lambda=abs(A);

% plot results as a hot density plot, with an additional (green) contour 
% indicating lambda = 1
pcolor(K/dx,Phi,Lambda), shading flat, colormap hot, colorbar
hold on, contour(K/dx,Phi,Lambda,[1 1],'g','LineWidth',2), hold off
xlabel('Wave number'), ylabel('\phi')
title('Stability of Lax-Wendroff')

% now for phase speed ratio
Dt=dx/c*Phi;   % matrix of delta ts
Omega=-atan2(imag(A),real(A))./Dt;
Phase_Ratio=dx/c*Omega./K;      % because the K has a factor of dx in it

% density plot
figure
pcolor(K/dx,Phi,Phase_Ratio), shading flat, colormap parula, colorbar
hold on, contour(K/dx,Phi,Phase_Ratio,[1 1],'k','LineWidth',2), hold off
xlabel('Wave number'), ylabel('\phi')
title('Phase Speed Ratio for Lax-Wendroff')


% calculate von Neumann stability:


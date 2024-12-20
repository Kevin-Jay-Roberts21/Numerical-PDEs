% Stability and Phase Analysis of 3rd Order Upwind Method
% Kevin Roberts
% Written Spetember 2024

clear all
close all
clc

L=10; J=100; dx=L/J; c=1;

pi*J/L

% define a vector of wave numbers
khat=2*pi*[1:J/2]/J;    % equivalent to 2*pi*k/L*dx

length(khat)

% define a vector of phi to test
nphi=50; phi=linspace(1/nphi,2,nphi); % we start at a small phi rather than zero 
                        
% define a 2D grid of k, phi for density plotting
[K,Phi]=meshgrid(khat,phi);

a_neg2 = 1/6*(Phi.^3 - Phi);
a_neg1 = 1/2*(2*Phi + Phi.^2 - Phi.^3);
a_0 = 1 - 1/2*Phi - Phi.^2 + 1/2*Phi.^3;
a_1 = 1/6*(-2*Phi + 3*Phi.^2 - Phi.^3);

% calculate the (complex) von Neumann multiplier 3rd Order Upwind
A = a_1.*exp(i*K) + a_0 + a_neg1.*exp(-i*K) + a_neg2.*exp(-2*i*K);

% % FOR LAX-WENDROFF
% % calculate the (complex) von Neumann multiplier Lax-Wendroff
%A = 1 -.5*Phi.*(exp(i*K) -exp(-i*K)) ...
    %+.5*Phi.^2.*(exp(i*K) +exp(-i*K) -2);

Lambda=abs(A);

% plot results as a hot density plot, with an additional (green) contour 
% indicating lambda = 1
pcolor(K/dx,Phi,Lambda), shading flat, colormap hot, colorbar
hold on, contour(K/dx,Phi,Lambda,[1 1],'g','LineWidth',2), hold off
xlabel('Wave number'), ylabel('\phi')
% title('Stability of Lax-Wendroff')
title('Stability of 3rd Order Upwind')

% now for phase speed ratio
Dt=dx/c*Phi;   % matrix of delta ts
Omega=-atan2(imag(A),real(A))./Dt;
Phase_Ratio=dx/c*Omega./K;      % because the K has a factor of dx in it

% density plot
figure
pcolor(K/dx,Phi,Phase_Ratio), shading flat, colormap parula, colorbar
hold on, contour(K/dx,Phi,Phase_Ratio,[1 1],'b','LineWidth',2), hold off
xlabel('Wave number'), ylabel('\phi')
% title('Phase Speed Ratio for Lax-Wendroff')
title('Phase Speed Ratio for 3rd Order Upwind')
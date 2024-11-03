% FTCS on the heat equation with Robin Boundary Conditions
% Kevin Roberts
% 10-13-2024
clear all
close all
clc

% Set the parameters for the solution
L = 20; % length of solution domain (from zero)
J = 200; % number of grid points (not counting zero)
x = linspace(0, L, J+1); % set up a vector of x locations
dx = L/J; % size of a grid cell, delta x
D = 1; % diffusion coefficient
dt = 0.0016; % size of time step, delta t
p = D*dt/(dx^2) % defining the dimensionless combination to be p (rho)
T = 15; % duraction of simulation
N = round(T/dt); % number of times to iterate simuation to get to time T
a = 0.5; 

% defining the initial condition
i_c_1 = 2*sin(x.*(pi/L)) + sin(x.*(2*pi/L));

% modifying the boundary conditions
bL = 0; % left boundary condition

FTCS0 = i_c_1;
FTCS = FTCS0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FTCS Method on the Diffusion Equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on, plot(x, FTCS0, 'b'), hold off % plot the IC

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['FTCS method on the diffusion equation, \rho = ' num2str(D*dt/(dx^2))])

for n = 2:N % beginning the for time loop
    
    FTCSn(1) = FTCS(1) + p*((FTCS(1)*(1 + a*dx) - a*dx*T) - 2*FTCS(1));
    
    for j = 2:J % begin spatial loop; note that j = 1 is x = 0    
        
        FTCSn(j) = FTCS(j) + p*(FTCS(j-1) - 2*FTCS(j) + FTCS(j+1));
        
    end % ending the spatial loop
    
    FTCSn(J+1) = bL;
    
    % to see results as an animation uncomment the line below and comment
    % the if - end statement
    % plot(x,u0,'r',x,crankn,'b'), pause(.1) % plots the next calculated time slice and the IC for reference. Makes this an animation
    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,FTCSn,'b'), hold off % if so, add a plot of current solution to existing plot
    end
    
    FTCS = FTCSn; % reset present solution for next time step
end % end time loop

% Diffusion Solutions
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

% modifying the boundary conditions
b0 = 0; % left boundary condition
bL = 0; % right boundary condition

% modifying the initial condition
i_c = 2*sin(x.*(pi/L)) + sin(x.*(2*pi/L));

% setting the initial conditions for the various methods
exact0 = i_c;
FTCS0 = i_c;
exact = exact0;
FTCS = FTCS0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FTCS Method on the Diffusion Equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(x, FTCS0, 'b') % plot the IC

% if we want to compare with the exact solution
plot(x, exact0, 'r')

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['FTCS method on the diffusion equation, \rho = ' num2str(D*dt/(dx^2))])

for n = 2:N % beginning the for time loop
    
    t = n * dt;    % Current time
    
    for j = 1:J+1 % begin spatial loop; note that j = 1 is x = 0    
        
        % the if and the else if statement handle the boundary conditions
        if j == 1
            FTCSn(j) = b0;
        elseif j == J+1
            % used for no flux at the boundary:
            % FTCSn(j) = FTCS(j) + p*(FTCS(j-1) - 2*FTCS(j) + FTCS(j));
            FTCSn(j) = bL;
        else
            FTCSn(j) = FTCS(j) + p*(FTCS(j-1) - 2*FTCS(j) + FTCS(j+1));
        end
        
        % evaluating the exact solution
        % Compute the exact solution using Fourier series
        % Reset for next time step
        exact(j) = 0;
        for k = 1:2  % Only compute for n=1 and n=2
            A_k = (k == 1) * 2 + (k == 2) * 1;  % Coefficients for n=1 and n=2
            contribution = A_k * exp(-((k * pi / L)^2) * D * t) * sin(k * pi * x(j) / L);
            exact(j) = exact(j) + contribution;
        end
        
        
    end % ending the spatial loop
    
    % evaluating the exact solution (making it satisfy the boundary conds.)
    
    
    
    % to see results as an animation uncomment the line below and comment
    % the if - end statement
    % plot(x,u0,'r',x,crankn,'b'), pause(.1) % plots the next calculated time slice and the IC for reference. Makes this an animation
    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,exact, 'r'), hold off % for plotting the exact
        hold on, plot(x,FTCSn,'b'), hold off % if so, add a plot of current solution to existing plot
    end
    
    FTCS = FTCSn; % reset present solution for next time step
end % end time loop

% Crank-Nicolson on the heat equation
% Kevin Roberts
% October 2024

clear all
close all
clc

% FOR HOMEOGENOUS DIRCHLET BC
%%
% Set the parameters for the solution
L = 20; % length of solution domain (from zero)
J = 200; % number of grid points (not counting zero)
x = linspace(0, L, J+1); % set up a vector of x locations
dx = L/J; % size of a grid cell, delta x
D = 1; % diffusion coefficient
dt = 0.0016; % size of time step, delta t
p = D*dt/(dx^2); % defining the dimensionless combination to be p (rho)
T = 15; % duraction of simulation
N = round(T/dt); % number of times to iterate simuation to get to time T

% modifying the boundary conditions
b0 = 0; % left boundary condition
bL = 0; % right boundary condition

% Used for homogeneous bc
i_c = 2*sin(x*(pi/L)) + sin(x*(2*pi/L));

% setting the initial conditions for the various methods
exact = i_c;
FTCS = i_c';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FTCS Method on the Diffusion Equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on, plot(x, FTCS, 'b--'), hold off % plot the IC

% if we want to compare with the exact solution
hold on, plot(x, exact, 'r'), hold off

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['Crank Nicolson method on the diffusion equation, \rho = ' num2str(D*dt/(dx^2))])

% Define matrices A and B
diagA = (1 + p)*ones(J-1, 1); % diagonal
off_diagA = (-p/2)*ones(J-2, 1); % off diagonals
A = diag(diagA) + diag(off_diagA, 1) + diag(off_diagA, -1);

a = diagA;
b = off_diagA;
c = off_diagA;

diagB = (1 - p)*ones(J-1, 1); % diagonal
off_diagB = (p/2)*ones(J-2, 1); % off diagonals
B = diag(diagB) + diag(off_diagB, 1) + diag(off_diagB, -1);

for n = 1:N
    
    t = n*dt;

    % Compute d using the tridiagonal structure (as opposed to d = B * FTCS(2:end-1))
    u_inner = FTCS(2:end-1);
    d = zeros(size(u_inner));
    d(1) = diagB(1)*u_inner(1) + off_diagB(1)*u_inner(2);
    for j = 2:J-2
        d(j) = off_diagB(j-1)*u_inner(j-1) + diagB(j)*u_inner(j) + off_diagB(j)*u_inner(j+1);
    end
    d(end) = off_diagB(end)*u_inner(end-1) + diagB(end)*u_inner(end);
    
    FTCS(2:end-1) = tri_diag_sol(a, b, c, d); % updating FTCSn (the solution at time step n)
    
    % apply the boundary conditions
    FTCS(1) = b0;
    FTCS(end) = bL;
    
    FTCSn = FTCS;
    
    % evaluating the exact solution (for b0 = 0 and bL = 0)
    exactn = 2*exp(-pi^2*D*n*dt/L^2)*sin(x*(pi/L)) + exp(-4*pi^2*D*n*dt/L^2)*sin(x*(2*pi/L));
    
    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,exactn, 'r'), hold off % for plotting the exact
        hold on, plot(x,FTCSn,'b--'), hold off % if so, add a plot of current solution to existing plot
    end

    % Update for the next time step
    FTCS = FTCSn;
end

% FOR BC u(0, t) = 0 and u(L, t) = 2
%%
% Set the parameters for the solution
L = 20; % length of solution domain (from zero)
J = 200; % number of grid points (not counting zero)
x = linspace(0, L, J+1); % set up a vector of x locations
dx = L/J; % size of a grid cell, delta x
D = 1; % diffusion coefficient
dt = 0.0016; % size of time step, delta t
p = D*dt/(dx^2); % defining the dimensionless combination to be p (rho)
T = 15; % duraction of simulation
N = round(T/dt); % number of times to iterate simuation to get to time T

% modifying the boundary conditions
b0 = 0; % left boundary condition
bL = 2; % right boundary condition

% Used for homogeneous bc
i_c = 2*sin(x*(pi/L)) + sin(x*(2*pi/L));

% setting the initial conditions for the various methods
exact = i_c;
FTCS = i_c';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FTCS Method on the Diffusion Equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on, plot(x, FTCS, 'b--'), hold off % plot the IC

% if we want to compare with the exact solution
hold on, plot(x, exact, 'r'), hold off

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['Crank Nicolson method on the diffusion equation, \rho = ' num2str(D*dt/(dx^2))])

% Define matrices A and B
diagA = (1 + p)*ones(J-1, 1); % diagonal
off_diagA = (-p/2)*ones(J-2, 1); % off diagonals
A = diag(diagA) + diag(off_diagA, 1) + diag(off_diagA, -1);

a = diagA;
b = off_diagA;
c = off_diagA;

diagB = (1 - p)*ones(J-1, 1); % diagonal
off_diagB = (p/2)*ones(J-2, 1); % off diagonals
B = diag(diagB) + diag(off_diagB, 1) + diag(off_diagB, -1);

for n = 1:N
    
    t = n*dt;

    % Compute d using the tridiagonal structure (as opposed to d = B * FTCS(2:end-1))
    u_inner = FTCS(2:end-1);
    d = zeros(size(u_inner));
    d(1) = diagB(1)*u_inner(1) + off_diagB(1)*u_inner(2);
    for j = 2:J-2
        d(j) = off_diagB(j-1)*u_inner(j-1) + diagB(j)*u_inner(j) + off_diagB(j)*u_inner(j+1);
    end
    d(end) = off_diagB(end)*u_inner(end-1) + diagB(end)*u_inner(end);
    
    FTCS(2:end-1) = tri_diag_sol(a, b, c, d); % updating FTCSn (the solution at time step n)
    
    % apply the boundary conditions
    FTCS(1) = b0;
    FTCS(end) = bL;
    
    FTCSn = FTCS;
    
    % evaluating the exact solution (for b0 = 0 and bL = 0)
    exactn = 2*exp(-pi^2*D*n*dt/L^2)*sin(x*(pi/L)) + exp(-4*pi^2*D*n*dt/L^2)*sin(x*(2*pi/L));
    
    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,exactn, 'r'), hold off % for plotting the exact
        hold on, plot(x,FTCSn,'b--'), hold off % if so, add a plot of current solution to existing plot
    end

    % Update for the next time step
    FTCS = FTCSn;
end



% FOR NO FLUX BC at x=L
%%



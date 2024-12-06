% Crank-Nicolson on the heat equation
% Kevin Roberts
% October 2024



% FOR HOMEOGENOUS DIRCHLET BC
%%

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
    
    FTCS(2:end-1) = tri_diag_sol(a, b, c, d); 

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
% Parameters
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
hold on, plot(x, FTCS, 'r'), hold off % plot the IC

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['Crank Nicolson method on the diffusion equation, \rho = ' num2str(D*dt/(dx^2))])

% Define matrices A and B
diagA = (1 + p)*ones(J+1, 1); % diagonal
off_diagA = (-p/2)*ones(J, 1); % off diagonals
A = diag(diagA) + diag(off_diagA, 1) + diag(off_diagA, -1);

diagB = (1 - p)*ones(J+1, 1); % diagonal
sup_diagB = (p/2)*ones(J, 1); % super off diagonals
sub_diagB = (p/2)*ones(J, 1); % sub off diagonals

% account for boundary conditions (need to change the top and bottom row of A and B)
A(1, 1) = 1;
A(1, 2) = 0;
A(J+1, J+1) = 1;
A(J+1, J) = 0;

diagB(1) = 1;
diagB(J+1) = p;
sup_diagB(1) = 0;
sub_diagB(J) = 0;

B = diag(diagB) + diag(sup_diagB, 1) + diag(sub_diagB, -1);

% defining vectors that will be plugged into the tridiagonal function
a = diagA;
b = off_diagA;
c = off_diagA;

for n = 1:N
    
    t = n*dt;

    % Compute d using the tridiagonal structure (as opposed to d = B * FTCS(2:end-1))
    
    d(1) = b0; % diagB(1)*FTCS(2) + sup_diagB(1)*FTCS(1);
    for j = 2:J
        d(j) = sup_diagB(j-1)*FTCS(j-1) + diagB(j)*FTCS(j) + sub_diagB(j-1)*FTCS(j+1);
    end
    d(J+1) = bL; % sub_diagB(J)*FTCS(J) + diagB(J+1)*FTCS(J+1);
    
    FTCSn = tri_diag_sol(a, b, c, d); 

    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,FTCSn,'b'), hold off % if so, add a plot of current solution to existing plot
    end

    % Update for the next time step
    FTCS = FTCSn;
end


% FOR NO FLUX BC at x=L
%%
% Parameters
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
p = D*dt/(dx^2); % defining the dimensionless combination to be p (rho)
T = 15; % duraction of simulation
N = round(T/dt); % number of times to iterate simuation to get to time T

% Used for homogeneous bc
i_c = 2*sin(x*(pi/(2*L))) + 3*sin(x*(3*pi/(2*L)));

% setting the initial conditions for the various methods
exact = i_c;
FTCS = i_c';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FTCS Method on the Diffusion Equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on, plot(x, FTCS, 'r'), hold off % plot the IC

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['Crank Nicolson method on the diffusion equation, \rho = ' num2str(D*dt/(dx^2))])

% Define matrices A and B
diagA = (1 + p)*ones(J+1, 1); % diagonal
off_diagA = (-p/2)*ones(J, 1); % off diagonals
A = diag(diagA) + diag(off_diagA, 1) + diag(off_diagA, -1);

diagB = (1 - p)*ones(J+1, 1); % diagonal
sup_diagB = (p/2)*ones(J, 1); % super off diagonals
sub_diagB = (p/2)*ones(J, 1); % sub off diagonals

% account for boundary conditions (need to change the top and bottom row of A and B)
% A(1, 1) = 1;
% A(1, 2) = 0;
%A(J+1, J+1) = 1;
A(J+1, J) = p;

% diagB(1) = 1;
% diagB(J+1) = 1;
% sup_diagB(1) = 0;
%sub_diagB(J) = p;

% defining vectors that will be plugged into the tridiagonal function
a = diagA;
b = off_diagA;
c = off_diagA;

for n = 1:N
    
    t = n*dt;

    % Compute d using the tridiagonal structure (as opposed to d = B * FTCS(2:end-1))
    
    d(1) = diagB(1)*FTCS(2) + sup_diagB(1)*FTCS(1);
    for j = 2:J
        d(j) = sup_diagB(j-1)*FTCS(j-1) + diagB(j)*FTCS(j) + sub_diagB(j-1)*FTCS(j+1);
    end
    d(J+1) = sub_diagB(J)*FTCS(J) + diagB(J+1)*FTCS(J+1);
    
    FTCSn = tri_diag_sol(a, b, c, d); 

    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,FTCSn,'b'), hold off % if so, add a plot of current solution to existing plot
    end

    % Update for the next time step
    FTCS = FTCSn;
end


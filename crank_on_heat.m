% Crank-Nicolson on the heat equation
% Kevin Roberts
% October 2024

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
i_c_1 = 2*sin(x.*(pi/L)) + sin(x.*(2*pi/L));

% % Used for b0 =0 and bL = 2
% i_c_2 = x.*2/L + (2 - 4/pi)*sin(x.*(pi/L)) + sin(x.*(2*pi/L)) + 4/(9*pi)*sin(x.*(3*pi/L));
% 
% % Used for no-flux boundary and different IC
% i_c_3 = 2*sin(x.*(pi/(2*L))) + sin(x.*(3*pi/(2*L)));

% setting the initial conditions for the various methods
exact0 = i_c_1;
FTCS0 = i_c_1;
exact = exact0;
FTCS = FTCS0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FTCS Method on the Diffusion Equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on, plot(x, FTCS0, 'b'), hold off % plot the IC

% if we want to compare with the exact solution
hold on, plot(x, exact0, 'r'), hold off

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['FTCS method on the diffusion equation, \rho = ' num2str(D*dt/(dx^2))])


% Define matrices A and B for interior points (199x199)
diagA = (1 + p)*ones(J+1, 1);
off_diagA = (-p/2)*ones(J, 1);
A = diag(diagA) + diag(off_diagA, 1) + diag(off_diagA, -1);

diagB = (1 - p)*ones(J+1, 1);
off_diagB = (p/2)*ones(J, 1);
B = diag(diagB) + diag(off_diagB, 1) + diag(off_diagB, -1);

% setting the boundary conditions
% NEED TO ADDRESS THIS
A(1,:) = 0;
A(1,1) = b0;
B(1,:) = 0;
B(1,1) = bL;

size(A)
size(B)

for n = 2:N
    
    t = n*dt;

    rhs = B * FTCS(2:J)';           % Transpose FTCS(2:J) to a column vector (199 elements)
    FTCSn = zeros(J+1, 1);          % Initialize FTCSn
    FTCSn(2:J) = A \ rhs;           % Solve for interior points (199 elements)
    
    % Boundary conditions
    FTCSn(1) = b0;
    FTCSn(end) = bL;
    
    % evaluating the exact solution (for b0 = 0 and bL = 0)
    exactn = 2*exp(-pi^2*D*n*dt/L^2)*sin(x.*(pi/L)) + exp(-4*pi^2*D*n*dt/L^2)*sin(x.*(2*pi/L));

    % % evaluating the exact solution (for b0 = 0 and bL = 0)
    % exactn = x.*(2/L) + (2-4\pi)*exp(-pi^2*D*n*dt/L^2)*sin(x.*(pi/L)) + exp(-4*pi^2*D*n*dt/L^2)*sin(x.*(2*pi/L)) + 4/(9*pi)*exp(-9*pi^2*D*n*dt/L^2)*sin(x.*(3*pi/L))

    % % evaluating the exact solution for no-flux
    % exactn = 2*exp(-pi^2*D*n*dt/(4*L^2))*sin(x.*(pi/(2*L))) + exp(-9*pi^2*D*n*dt/(4*L^2))*sin(x.*(3*pi/(2*L)))
    
    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,exactn, 'r'), hold off % for plotting the exact
        hold on, plot(x,FTCSn,'b'), hold off % if so, add a plot of current solution to existing plot
    end

    % Update for the next time step
    FTCS = FTCSn;

end






% ATTMEPT TO SOLVE WITHOUT DOING MATRIX INVERSE
% for n = 2:N % beginning the for time loop
% 
%     t = n * dt;    % Current time
% 
%     for j = 1:J+1 % begin spatial loop; note that j = 1 is x = 0    
% 
%         % the if and the else if statement handle the boundary conditions
%         if j == 1
%             FTCSn(j) = b0;
% 
%         elseif j == J+1
%             FTCSn(j) = bL;
%             % % for the no flux boundary
%             % FTCSn(j) = FTCS(j) + p*(FTCS(j-1) - 2*FTCS(j) + FTCS(j));
% 
%         elseif j == 2
%             FTCSn(j) = -p/2*((1 - p)*FTCS(j-1) + p/2*FTCS(j) - (1 + p)*FTCSn(j-1));
% 
%         elseif j == J
%             FTCSn(j) = 1/(1+p)*(p/2*FTCS(j-1) + (1-p)*FTCS(j) + p/2*FTCSn(j-1));
% 
%         else
%             FTCSn(j) = -p/2*(p/2*FTCS(j-2) + (1 - p)*FTCS(j-1) + p/2*FTCS(j) + p/2*FTCSn(j-2) - (1 + p)*FTCSn(j-1));
% 
%         end
% 
%     end % ending the spatial loop
% 
%     % % evaluating the exact solution (for b0 = 0 and bL = 0)
%     exactn = 2*exp(-pi^2*D*n*dt/L^2)*sin(x.*(pi/L)) + exp(-4*pi^2*D*n*dt/L^2)*sin(x.*(2*pi/L));
% 
%     % % evaluating the exact solution (for b0 = 0 and bL = 0)
%     % exactn = x.*(2/L) + (2-4\pi)*exp(-pi^2*D*n*dt/L^2)*sin(x.*(pi/L)) + exp(-4*pi^2*D*n*dt/L^2)*sin(x.*(2*pi/L)) + 4/(9*pi)*exp(-9*pi^2*D*n*dt/L^2)*sin(x.*(3*pi/L))
% 
%     % % evaluating the exact solution for no-flux
%     % exactn = 2*exp(-pi^2*D*n*dt/(4*L^2))*sin(x.*(pi/(2*L))) + exp(-9*pi^2*D*n*dt/(4*L^2))*sin(x.*(3*pi/(2*L)))
% 
%     % to see results as an animation uncomment the line below and comment
%     % the if - end statement
%     % plot(x,u0,'r',x,crankn,'b'), pause(.1) % plots the next calculated time slice and the IC for reference. Makes this an animation
%     if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
%         hold on, plot(x,exactn, 'r'), hold off % for plotting the exact
%         hold on, plot(x,FTCSn,'b'), hold off % if so, add a plot of current solution to existing plot
%     end
% 
%     FTCS = FTCSn; % reset present solution for next time step
% end % end time loop

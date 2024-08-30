clear all
close all
clc

%
% Introductory Matlab script for solving the advection equation: 
% u_t + c*u_x = 0
% using simple upwinding
% 

% Set the parameters for the solution
L = 20; % length of solution domain (from zero)
J = 200; % number of grid points (not counting zero)
x = linspace(0, L ,J+1); % set up a vector of x locations
dx = L/J; % size of a grid cell, delta x
c = 1; % advection speed
dt = 0.08; % size of time step, delta t

T = 15; % duraction of simulation
N = round(T/dt); % number of times to iterate simuation to get to time T

% set an initial condition
u0 = x./(x.^4 + 1);
d0 = x./(x.^4 + 1);
c0 = x./(x.^4 + 1);

up = u0; % up will denote u in the Present during the simulation
down = d0;
centered = c0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the Upwind Method on u_x + cu_t = 0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot(x, u0, 'r') % plot the IC
% 
% % add labels to the plot
% xlabel('x'), ylabel('u(x,t)')
% title(['Upwinding for the advection equation, \phi = ' num2str(c*dt/dx)])
% 
% for n = 1:N % beginning the for time loop
% 
%     up(1) = 0; % enforce the zero BC at x = 0
% 
%     for j = 2:J+1 % begin spatial loop; note that j = 1 is x = 0
%         un(j) = up(j) - c*dt/dx*(up(j)-up(j-1));
%     end % ending the spatial loop
% 
%     % to see results as an animation uncomment the line below and comment
%     % the if - end statement
% 
%     plot(x,u0,'r',x,un,'b'), pause(.1) % plots the next calculated time slice and the IC for reference. Makes this an animation
% 
%     % if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
%     %     hold on, plot(x,un,'b'), hold off % if so, add a plot of current solution to existing plot
%     % end
% 
%     up = un; % reset present solution for next time step
% end % end time loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the Downwind Method on u_x + cu_t = 0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot(x, d0, 'r') % plot the IC
% 
% % add labels to the plot
% xlabel('x'), ylabel('u(x,t)')
% title(['Downwinding for the advection equation, \phi = ' num2str(c*dt/dx)])
% 
% for n = 1:N % beginning the for time loop
% 
%     down(1) = 0; % enforce the zero BC at x = 0
% 
%     for j = 2:J+1 % begin spatial loop; note that j = 1 is x = 0
%         dn(j) = down(j) - c*dt/dx*(down(j-1)-down(j));
%     end % ending the spatial loop
% 
%     % to see results as an animation uncomment the line below and comment
%     % the if - end statement
% 
%     plot(x,d0,'r',x,dn,'b'), pause(.1) % plots the next calculated time slice and the IC for reference. Makes this an animation
% 
%     % if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
%     %     hold on, plot(x,dn,'b'), hold off % if so, add a plot of current solution to existing plot
%     % end
% 
%     down = dn; % reset present solution for next time step
% end % end time loop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the Centered Method on u_x + cu_t = 0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(x, c0, 'r') % plot the IC

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['Centering Method for the advection equation, \phi = ' num2str(c*dt/dx)])

for n = 1:N % beginning the for time loop

    centered(1) = 0; % enforce the zero BC at x = 0
    for j = 2:J % begin spatial loop; note that j = 1 is x = 0
        cn(j) = centered(j) - c*dt/(2*dx)*(centered(j+1)-centered(j-1));
    end % ending the spatial loop

    centered_of_J_plus_1 = 20.1/(20.1^4 + 1); % maybe this is just 0?
    cn(201) = centered(J) - c*dt/(2*dx)*(centered_of_J_plus_1 - centered(J-1)); % setting anther boundary condition

    % to see results as an animation uncomment the line below and comment
    % the if - end statement

    plot(x,c0,'r',x,cn,'b'), pause(.1) % plots the next calculated time slice and the IC for reference. Makes this an animation

    % if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
    %     hold on, plot(x,cn,'b'), hold off % if so, add a plot of current solution to existing plot
    % end
    % 
    centered = cn; % reset present solution for next time step
end % end time loop
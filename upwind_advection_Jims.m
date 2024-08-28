% 
% Introductory Matlab script for a solving the advection equation 
%     u_t + c u_x =0
% using simple upwinding
%
% Jim Powell, 22 August 2024
%

% Set parameters for the solution

L=20;   % length of solution domain (from zero)
J=200;  % number of grid points (not counting zero)
x=linspace(0,L,J+1);    % set up a vector of x locations
dx=L/J; % size of a grid cell, delta x
c=1;    % advection speed
dt=.08; % size of time step, delta t

T=15;   % duration of simulation
N=round(T/dt);  % number of times to iterate simulation to get to time T

% set an initial condition
u0=x./(x.^4+1);   
up=u0;  % up will denote u in the Present during the simulation

plot(x,u0,'r')  % plot the IC
% add labels to plot
xlabel('x'), ylabel('u(x,t)')
title(['Upwinding for the advection equation, \phi = ' num2str(c*dt/dx) ])

for n=1:N       % begin the loop for time
    up(1)=0;    % enforce the zero BC at x=0
    for j=2:J+1     % begin spatial loop; note that j=1 is x=0
        un(j)=up(j)-c*dt/dx*(up(j)-up(j-1));
    end             % end spatial loop
    
    % to see results as an animation uncomment the line below and comment
    %   the if - end statement
    % plot(x,u0,'r',x,un,'b'), pause(.1)  % plots the next calculated time slice
                    % and the IC for reference.  Makes this an animation
                    
    if  mod(n*dt,2) < dt  % check if current time is close to a multiple of 2
        hold on, plot(x,un,'b'), hold off  % if so, add a plot of current solution to existing plot
    end
    
    up=un;      % reset present solution for next time step
end             % end time loop


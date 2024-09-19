% A MATLAB code that includes multiple methods of solving the advection equation.
% Written by Kevin Roberts
% Last updated on September 2024

clear all
close all
clc

% Set the parameters for the solution
L = 20; % length of solution domain (from zero)
J = 200; % number of grid points (not counting zero)
x = linspace(0, L, J+1); % set up a vector of x locations
dx = L/J; % size of a grid cell, delta x
c = 1; % advection speed (downward is stable is this is less than 0)
dt = 0.08; % size of time step, delta t
phi = c*dt/dx % phase speed

T = 15; % duraction of simulation
N = round(T/dt); % number of times to iterate simuation to get to time T

% setting the initial conditions for the various methods
u0 = x./(x.^4 + 1);
d0 = x./(x.^4 + 1);
c0 = x./(x.^4 + 1);
leap0 = x./(x.^4 + 1);
crank0 = x./(x.^4 + 1);

up = u0; 
down = d0;
centered = c0;
leapfrog = leap0;
crank = crank0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the Upwind Method on u_x + cu_t = 0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(x, u0, 'r') % plot the IC

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['Upwinding for the advection equation, \phi = ' num2str(c*dt/dx)])

for n = 1:N % beginning the for time loop

    up(1) = 0; % enforce the zero BC at x = 0

    for j = 2:J+1 % begin spatial loop; note that j = 1 is x = 0
        un(j) = up(j) - c*dt/dx*(up(j)-up(j-1));
    end % ending the spatial loop

    % to see results as an animation uncomment the line below and comment
    % the if - end statement

    % plot(x,u0,'r',x,un,'b'), pause(.1) % plots the next calculated time slice and the IC for reference. Makes this an animation

    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,un,'b'), hold off % if so, add a plot of current solution to existing plot
    end

    up = un; % reset present solution for next time step
end % end time loop

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the Downwind Method on u_x + cu_t = 0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(x, d0, 'r') % plot the IC

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['Downwinding for the advection equation, \phi = ' num2str(c*dt/dx)])

for n = 1:N % beginning the for time loop

    down(1) = 0; % enforce the zero BC at x = 0

    for j = 2:J % begin spatial loop; note that j = 1 is x = 0
        dn(j) = down(j) - c*dt/dx*(down(j+1)-down(j));
    end % ending the spatial loop

    dn(J+1) = 0; % setting another boundary condition

    % to see results as an animation uncomment the line below and comment
    % the if - end statement

    % plot(x,d0,'r',x,dn,'b'), pause(.1) % plots the next calculated time slice and the IC for reference. Makes this an animation

    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,dn,'b'), hold off % if so, add a plot of current solution to existing plot
    end

    down = dn; % reset present solution for next time step
end % end time loop

%%
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

    cn(J+1) = 0; % setting another boundary condition

    % to see results as an animation uncomment the line below and comment
    % the if - end statement

    % plot(x,c0,'r',x,cn,'b'), pause(.1) % plots the next calculated time slice and the IC for reference. Makes this an animation

    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,cn,'b'), hold off % if so, add a plot of current solution to existing plot
    end
    % 
    centered = cn; % reset present solution for next time step
end % end time loop

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the LeapFrog Method on u_x + cu_t = 0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(x, leap0, 'r') % plot the IC

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['Leap Frog Method for the advection equation, \phi = ' num2str(c*dt/dx)])

leap_all = leapfrog;

for n = 1:N % beginning the for time loop

    leapfrog(1) = 0; % enforce the zero BC at x = 0

    % using the upwind method for the first iteration
    if n == 1
        for j = 2:J+1 % begin spatial loop; note that j = 1 is x = 0
            leapn(j) = leapfrog(j) - c*dt/dx*(leapfrog(j)-leapfrog(j-1));
        end % ending the spatial loop
    else
        for j = 2:J % begin spatial loop; note that j = 1 is x = 0
            leapn(j) = leap_all((n+1)-2, j) - c*dt/(dx)*(leapfrog(j+1)-leapfrog(j-1));
        end % ending the spatial loop
    end

    ln(J+1) = 0; % setting another boundary condition

    % to see results as an animation uncomment the line below and comment
    % the if - end statement

    % plot(x,leap0,'r',x,leapn,'b'), pause(.01) % plots the next calculated time slice and the IC for reference. Makes this an animation

    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,leapn,'b'), hold off % if so, add a plot of current solution to existing plot
    end
    % 
    leap_all = [leap_all; leapn];
    leapfrog = leapn; % reset present solution for next time step
end % end time loop

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crank-Nicolson Method on the Advection Equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(x, crank0, 'r') % plot the IC
f = 0;

% add labels to the plot
xlabel('x'), ylabel('u(x,t)')
title(['Crank Nicolson for the advection equation, \phi = ' num2str(c*dt/dx)])

for n = 1:N % beginning the for time loop

    for j = 1:J+1 % begin spatial loop; note that j = 1 is x = 0
        
        if j == 1
            crankn(j) = (1-phi)/(1+phi)*crank(j) + f/(1+phi);
        else
            crankn(j) = -(1-phi)/(1+phi)*crankn(j-1) + (1-phi)/(1+phi)*crank(j) + crank(j-1)      
        end
    
    end % ending the spatial loop

    % to see results as an animation uncomment the line below and comment
    % the if - end statement

    % plot(x,u0,'r',x,crankn,'b'), pause(.1) % plots the next calculated time slice and the IC for reference. Makes this an animation

    if mod(n*dt,2) < dt % check if current time is close to a multiple of 2
        hold on, plot(x,crankn,'b'), hold off % if so, add a plot of current solution to existing plot
    end

    crank = crankn; % reset present solution for next time step
end % end time loop

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lax-Wendroff Method on the Advection Equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

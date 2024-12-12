% Solving the 2d Heat Eqaution with a no flux boundary conditon at y = L_y
% Kevin Roberts (Adapted from Jim Powell's code)
% December 2024

clear all
close all
clc

% set up spatial grid and discretization
Lx = 5; J = 50; dx = Lx/J; x = linspace(0, Lx, J+1);
Ly = 10; K = 100; dy = Ly/K; y = linspace(0, Ly, K+1);
[X,Y] = meshgrid(x,y); % this command sets up matrices of coordinate values
D = 1; dt = 0.24*dx^2/D; % dt just below stability boundary
T = 20; N = ceil(T/dt); % how long to run simulation
r1 = D*dt/dx^2; r2 = D*dt/dy^2; % dimensionless rhos

g1 = Ly/2 - abs(Ly/2 - y); % BC at x=0
g2 = sin(pi*y/Ly); % BC at x=Lx

% assume u=0 at y=0 and y=Ly

up = 0*X; % zero IC
un = 0*X; % just to initialize u^{n+1}
            
% Note - since the first and last rows are never updated below, this automatically sets zero boundary conditions

up(:, 1) = g1; % BC at x=0 (all time-independent)
up(:, J+1) = g2; % BC at x=Lx
un(:, 1) = g1; % BC at x=0     -- necessary since we will use up=un to update below
un(:, J+1) = g2; % BC at x=Lx    -- necessary since we will use up=un to update below

for n = 1:N % begin loop in time 
    for j = 2:J
        for k = 2:K
            un(k, j) = up(k, j) + r1*(up(k, j+1) + up(k, j-1) - 2*up(k, j))+...
                r2*(up(k+1, j) + up(k-1, j) - 2*up(k, j));
        end     
    end
    
    % applying the no-flux boundary condition at y = L_y
    un(K+1, :) = up(K+1, j) + r1*(up(K+1, j+1) + up(K+1, j-1) - 2*up(K+1, j))+...
                r2*(2*up(K, j) - 2*up(K+1, j));
    
    up = un;  % update for next time step

    if (mod(n*dt, 0.1) < dt) % plot every .1 sec
        pcolor(X, Y, up), colormap hot, shading flat, caxis([0 Ly/2]), axis image
            % set color axis so that it maxes at max of BC
            pause(0.1) % pause for rendering the plot
    end

end

disp("Code Finished")





% Using the ADI method to solve the 2D diffusion equation
% Kevin Roberts
% December 2024

clear all
close all
clc

% set up spatial grid and discretization
Lx = 10; J = 100; dx = Lx/J; x = linspace(0, Lx, J+1);
Ly = 5; K = 50; dy = Ly/K; y = linspace(0, Ly, K+1);
[X,Y] = meshgrid(x, y); % this command sets up matrices of coordinate values
D = 1; dt = 0.01; % dt just below stability boundary
T = 4; N = ceil(T/dt); % how long to run simulation
r1 = D*dt/dx^2; r2 = D*dt/dy^2; % dimensionless rhos

g1 = Ly/2 - abs(Ly/2 - y); % BC at x=0
g2 = sin(pi*y/Ly); % BC at x=Lx

up = (abs(X - 5) <= 2).*(abs(Y - 2.5) <= 1); % IC % zero IC

un = 0*X; % just to initialize u^{n+1}
rhs1 = 0*X;
rhs2 = 0*X;
u_star = 0*X;

% creating the matrices
diagA1 = (1 + r1)*ones(J-1, 1); % diagonal of A
sup_diagA1 = (-r1/2)*ones(J-2, 1); % super diagonal of A
sub_diagA1 = (-r1/2)*ones(J-2, 1); % sub diagonals of A
A1 = diag(diagA1, 0) + diag(-sup_diagA1, 1) + diag(-sub_diagA1, -1);

diagA2 = (1 + r2)*ones(K-1, 1); % diagonal of A
sup_diagA2 = (-r2/2)*ones(K-2, 1); % super diagonal of A
sub_diagA2 = (-r2/2)*ones(K-2, 1); % sub diagonals of A
A2 = diag(diagA2, 0) + diag(-sup_diagA2, 1) + diag(-sub_diagA2, -1);

m = 0;

for n = 1:N 
    % Find the RHS with loops for j and k
    for j = 2:J
        for k = 2:K
            rhs1(k, j) = (1 - r2)*up(k, j) + 1/2*r2*(up(k+1, j) + up(k-1, j));
        end     
    end
    
    % tridiagonalize with d = RHS and with a fixed k
    for k = 2:K
        d1 = rhs1(k, 2:J)';
        solution = tri_diag_sol(diagA1, sup_diagA1, sub_diagA1, d1);
        u_star(k, 2:J) = solution';
    end
    
    
    
    % Find the RHS with 2 loops
    for j = 2:J
        for k = 2:K
            rhs2(k, j) = (1 - r1)*u_star(k, j) + 1/2*r1*(u_star(k, j+1) + u_star(k, j-1));
        end     
    end
    
    % tridiagonalize with d = RHS and with a fixed k
    for j = 2:J
        d2 = rhs2(2:K, j);
        solution = tri_diag_sol(diagA2, sup_diagA2, sub_diagA2, d2);
        un(2:K, j) = solution;
    end
    
    up = un;
    
    if (mod(n*dt, 1) < dt)        % plot every .1 sec
        m = m + 1; figure(m)
        surf(X, Y, up), colormap hot, shading flat, axis image
            % set color axis so that it maxes at max of BC  %caxis([0 1]), 
            pause(.1)  % pause for rendering the plot
    end     % of plotting if
    
end







disp("Code Finished")





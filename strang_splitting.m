% Using Strang Splitting to solve PDEs
% Kevin Roberts
% November 2024



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STRANG-SPLITTING ON ADVECTION DIFFUSION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Parameters
L = 10;                % Domain length
J = 256;               % Number of spatial points
dx = 2*L/J;              % Spatial step size
x = linspace(-L, L-dx, J); % Spatial grid
c = 1;                 % Advection speed
D = 1;                 % Diffusion coefficient
dt = dx/c;             % Time step (CFL condition)
T = 1;               % Final time
steps = floor(T / dt); % Number of time steps

% Initial condition: block wave
u_initial = 10*(abs(x-3) < 2);

u = u_initial;

% plotting the initial condition
plot(x, u_initial, 'r'); hold on;

% Wavenumbers for Fourier transform
k = pi/L * [0:J/2 1-J/2:-1]; % FFT wavenumbers

% Time-stepping loop
for n = 1:steps
    % First diffusion half-step
    u_hat = fft(u);
    u_hat = u_hat .* exp(-D * (k.^2) * (dt/2));
    u = real(ifft(u_hat));
    
    % Advection full step
    u_p = u;
    u_new(1) = u_p(J); % for the periodic boundary conditions
    for j = 2:J
        u_new(j) = u_p(j) - c*dt/dx*(u_p(j) - u_p(j-1));
    end
    
    % Second diffusion half-step
    u_new_fft = fft(u_new);
    u_k = u_new_fft .* exp(-D * (k.^2) * (dt/2));
    u = real(ifft(u_k));
    
    
    plot(x, u, 'b'); 
    pause(0.05); % Pause to animate
end

% Plot the final solution
plot(x, u, 'LineWidth', 1.5);
xlabel('x');
ylabel('u');
title('Advection-Diffusion Solution (Square Wave)');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STRANG-SPLITTING ON FKPP EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Parameters
L = 2*100;                  % Domain length
J = 1000;                   % Number of spatial points
x = linspace(-L/2, L/2, J); % Spatial grid
dx = L/J;                   % Spatial step size
K = 1;                      
m = 1;                      
D = 1;
c = 1;
dt = dx/c;                  % Time step (CFL condition)
T = 20;                      % Final time
steps = floor(T / dt); % Number of time steps

% Initial condition: square wave
u_initial = exp(-x.^2); 

u = u_initial; % defining an initial condition

% plotting the initial condition
plot(x, u_initial, 'r'); hold on;

% Wavenumbers for Fourier transform
k = 2 * pi * [0:J/2 1-J/2:-1]/L; % FFT wavenumbers

% Time-stepping loop
for n = 1:steps
    % First diffusion half-step
    u_hat = fft(u);
    u_hat = u_hat .* exp(-D * (k.^2) * (dt/2));
    u = real(ifft(u_hat));
    
    % getting the analytical solution of logistic eqn
    v = f(u, K, m, dt);
    
    % Second diffusion half-step
    u_hat = fft(v);
    u_hat = u_hat .* exp(-D * (k.^2) * (dt/2));
    u = real(ifft(u_hat));
    if  mod(n*dt,2) < dt  % check if current time is close to a multiple of 2
        plot(x, u,'b')  % if so, add a plot of current solution to existing plot
    end
    
end

% Plot the final solution
plot(x, u, 'LineWidth', 1.5);
xlabel('x');
ylabel('u');
title('FKPP Solution');

% defining the f function as f(u) = mu(1 - u)
function v = f(u, K, m, dt)

    v = (K*u) ./ (((K - u) .* exp(-m * dt)) + u);

end


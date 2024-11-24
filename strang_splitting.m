% using Strang Splitting to sovle PDE
% Kevin Roberts
% November 2024

clear all
close all
clc

% Parameters
L = 10;                % Domain length
J = 256;               % Number of spatial points
x = linspace(0, L, J); % Spatial grid
dx = L/J;              % Spatial step size
c = 1;                 % Advection speed
D = 1;                 % Diffusion coefficient
dt = dx/c;             % Time step (CFL condition)
T = 2;               % Final time
steps = floor(T / dt); % Number of time steps

% Initial condition: square wave
u = zeros(size(x));
u(x < L/2) = 1.0; % Square wave

u_initial = u; % defining an initial condition

% plotting the initial condition
plot(x, u_initial, 'r'); hold on;

% Wavenumbers for Fourier transform
k = 2 * pi * [0:J/2-1 -J/2:-1] / L; % FFT wavenumbers

% Time-stepping loop
for n = 1:steps
    % First diffusion half-step
    u_hat = fft(u);
    u_hat = u_hat .* exp(-D * (k.^2) * (dt/2));
    u = real(ifft(u_hat));
    
    % Advection full step
    shift = round(c * dt / dx); % Number of grid points to shift
    u = circshift(u, -shift); % Shift left for periodic boundary

    % Second diffusion half-step
    u_hat = fft(u);
    u_hat = u_hat .* exp(-D * (k.^2) * (dt/2));
    u = real(ifft(u_hat));
    
    plot(x, u, 'b'); 
    pause(0.05); % Pause to animate
end

% Plot the final solution
plot(x, u, 'LineWidth', 1.5);
xlabel('x');
ylabel('u');
title('Advection-Diffusion Solution (Square Wave)');




clear all
close all
clc

J = 256;      % use a power of two for maximum efficiency
L = 10;       % half width of domain
dx = 2 * L / J;   % spatial step size
x = linspace(0, L - dx, J);  % FFT assumes f(L)=f(-L), so last point unnecessary
D = 1;        % diffusion constant

% set initial condition and take FFT
f = 10 * (abs(x - 3) < 2);     % block wave, width 4, centered at x=3

% Homogenous Dirichlet boundary conditions (u = 0 at x = -L and x = L)
% Reflect f with inverted sign at boundaries
% NOTE: the -fliplr(f) function takes the list, flips it, and makes it
% negative, for example if f = [1, 2, 3, 4], the -fliplr(f) = [-4, -3, -2, -1]
f_reflected_dirichlet = [-fliplr(f), f];

% Homogenous Neumann boundary conditions (du/dx = 0 at x = -L and x = L)
% Reflect f without inverting sign at boundaries
f_reflected_neumann = [fliplr(f), f];

% Mixed boundary conditions (Dirichlet at x = -L, Neumann at x = L)
f_mixed = [fliplr(f), f];  % Dirichlet reflection at -L
f_mixed(end/2+1:end) = fliplr(f);  % Neumann reflection at L

% Update ks to match the length of the expanded domain
ks_reflected = pi / L * [0:(2*J)/2 1-(2*J)/2:-1];  % Adjust for size 2*J

% Homogenous Dirichlet
figure;
fk = fft(f_reflected_dirichlet);  % transform!
plot(x, f, 'r'), hold on  % IC in red (only plot original domain)
for t = 1:5   % times for the various solutions
    uk = exp(-D * t * ks_reflected.^2).*fk;     % solution in wave space
    u = real(ifft(uk)); % inverse FFT, ignore small imaginary bits
    u = u(1:J); % extract original domain
    plot(x, u, 'b')
end
title('Homogenous Dirichlet Boundary Conditions')
xlabel('x'), ylabel('u(x,t)')
hold off

% Homogenous Neumann
figure;
fk = fft(f_reflected_neumann);  % transform!
plot(x, f, 'r'), hold on  % IC in red (only plot original domain)
for t = 1:5   % times for the various solutions
    uk = exp(-D * t * ks_reflected.^2) .* fk;     % solution in wave space
    u = real(ifft(uk)); % inverse FFT, ignore small imaginary bits
    u = u(1:J); % extract original domain
    plot(x, u, 'b')
end
title('Homogenous Neumann Boundary Conditions')
xlabel('x'), ylabel('u(x,t)')
hold off

% Mixed boundary conditions
figure;
fk = fft(f_mixed);  % transform!
plot(x, f_mixed, 'r'), hold on  % IC in red (only plot original domain)
for t = 1:5   % times for the various solutions
    uk = exp(-D * t * ks_reflected.^2) .* fk;     % solution in wave space
    u = real(ifft(uk)); % inverse FFT, ignore small imaginary bits
    u = u(1:J); % extract original domain
    plot(x, u, 'b')
end
title('Mixed Boundary Conditions')
xlabel('x'), ylabel('u(x,t)')
hold off
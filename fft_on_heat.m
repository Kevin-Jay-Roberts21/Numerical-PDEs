clear all
close all
clc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Homogeneous Dirichlet %
%%%%%%%%%%%%%%%%%%%%%%%%%

J = 256;      % use a power of two for maximum efficiency
L = 10;       % half width of domain
dx = 2 * L / J;   % spatial step size
x = linspace(-L, L - dx, J);  % FFT assumes f(L)=f(-L), so last point unnecessary
D = 1;        % diffusion constant

f = zeros(size(x));

% set initial condition and take FFT
f(x >= 0) = 10 * (abs(x(x >= 0) - 3) < 2); % f for x >= 0
f(x < 0) = -10 * (abs(x(x < 0) + 3) < 2); % f for x < 0

% Update ks to match the length of the expanded domain
ks = pi / L * [0:J/2 1-J/2:-1];  % Adjust for size 2*J

% Homogenous Dirichlet
figure;
fk = fft(f);  % transform!
plot(x, f, 'r'), hold on  % IC in red (only plot original domain)
for t = 1:5   % times for the various solutions
    uk = exp(-D * t * ks.^2).*fk;     % solution in wave space
    u = real(ifft(uk)); % inverse FFT, ignore small imaginary bits
    u = u(1:J); % extract original domain
    plot(x, u, 'b')
end
title('Homogenous Dirichlet Boundary Conditions')
xlabel('x'), ylabel('u(x,t)')
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%
% Homogeneous Neumann %
%%%%%%%%%%%%%%%%%%%%%%%

J = 256;      % use a power of two for maximum efficiency
L = 10;       % half width of domain
dx = 2 * L / J;   % spatial step size
x = linspace(-L, L - dx, J);  % FFT assumes f(L)=f(-L), so last point unnecessary
D = 1;        % diffusion constant

f = zeros(size(x));

% set initial condition and take FFT
f(x >= 0) = 10 * (abs(x(x >= 0) - 3) < 2); % f for x >= 0
f(x < 0) = 10 * (abs(x(x < 0) + 3) < 2); % f for x < 0

% Update ks to match the length of the expanded domain
ks = pi / L * [0:J/2 1-J/2:-1];  % Adjust for size 2*J

% Homogenous Neumann
figure;
fk = fft(f);  % transform!
plot(x, f, 'r'), hold on  % IC in red (only plot original domain)
for t = 1:5   % times for the various solutions
    uk = exp(-D * t * ks.^2).*fk;     % solution in wave space
    u = real(ifft(uk)); % inverse FFT, ignore small imaginary bits
    u = u(1:J); % extract original domain
    plot(x, u, 'b')
end
title('Homogenous Neumann Boundary Conditions')
xlabel('x'), ylabel('u(x,t)')
hold off

%%
%%%%%%%%%
% Mixed %
%%%%%%%%%

J = 2*256;      % use a power of two for maximum efficiency
L = 10;       % half width of domain
dx = 4 * L / J;   % spatial step size
x = linspace(-2*L, 2*L - dx, J);  % FFT assumes f(L)=f(-L), so last point unnecessary
D = 1;        % diffusion constant

f = zeros(size(x));

% set initial condition and take FFT
f(abs(x - 17) < 2 & x >= 0) = -10;
f(abs(x - 3) < 2 & x > 0) = 10;
f(abs(x + 3) < 2 & x < 0) = 10;
f(abs(x + 17) < 2 & x < 0) = -10;

% Update ks to match the length of the expanded domain
ks = pi / L * [0:J/2 1-J/2:-1];  % Adjust for size 2*J

% Mixed boundary conditions
figure;
fk = fft(f);  % Transform!

length(ks)
length(fk)
plot(x, f, 'r'), hold on  % Plot IC in red (only plot original domain)
for t = 1:5   % Times for the various solutions
    uk = exp(-D * t * ks.^2).*fk;  % Solution in wave space
    u = real(ifft(uk));  % Inverse FFT, ignore small imaginary bits
    u = u(1:J);  % Extract original domain
    plot(x, u, 'b')
end
title('Mixed Boundary Conditions')
xlabel('x'), ylabel('u(x,t)')
hold off
% Using FFT to solve the Advection Eqn
% Kevin Roberts
% November

clear all
close all
clc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE FOR u(x,0) = 10*exp(-(x-3)^2/8) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = 256;      % use a power of two for maximum efficiency
L = 10;       % half width of domain
dx = 2 * L/J;   % spatial step size
x = linspace(-L, L-dx, J);  % FFT assumes f(L)=f(-L) so last point unnecessary
c = 1;        % diffusion constant

% set initial condition and take FFT
f = 10 * exp(-(x-3).^2/8);     % block wave, width 4, centered at x=3
fk = fft(f);                 % transform!

% set up vector of wave numbers
ks = pi/L * [0:J/2 1-J/2:-1]; % where k = [0:J/2 1-J/2:-1]

% calculate and plot spectral solution at different times
plot(x, f, 'r')   % IC in red

for t=1:5   % times for the various solutions
    uk=exp(-i * c * t * ks).*fk;     % solution in wave space
    u=real(ifft(uk)); % inverse FFT, ignore small imaginary bits
    hold on, plot(x, u, 'b'), hold off
end

xlabel('x'), ylabel('u(x,t)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE FOR u(x,0) = 10*(H(x-1) - H(x-5)) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = 256;      % use a power of two for maximum efficiency
L = 10;       % half width of domain
dx = 2 * L/J;   % spatial step size
x = linspace(-L, L-dx, J);  % FFT assumes f(L)=f(-L) so last point unnecessary
c = 1;        % diffusion constant

% set initial condition and take FFT
f = 10*(heaviside(x-1) - heaviside(x-5));     % block wave, width 4, centered at x=3
fk = fft(f);                 % transform!

% set up vector of wave numbers
ks = pi/L * [0:J/2 1-J/2:-1]; % where k = [0:J/2 1-J/2:-1]

% calculate and plot spectral solution at different times
plot(x, f, 'r')   % IC in red

for t=1:5   % times for the various solutions
    uk=exp(-i * c * t * ks).*fk;     % solution in wave space
    u=real(ifft(uk)); % inverse FFT, ignore small imaginary bits
    hold on, plot(x, u, 'b--'), hold off
end

xlabel('x'), ylabel('u(x,t)')
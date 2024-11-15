% Heat Equation Convolution
% Kevin Roberts
% November 2024

clear all
close all
clc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Homogeneous Dirichlet %
%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters and define domain
J = 256;      % use a power of two for maximum efficiency
L = 10;       % half width of domain
dx = 2*L/J;   % spatial step size
x = linspace(-L, L-dx, J);  % FFT assumes f(L)=f(-L) so last point unnecessary
D = 1;        % diffusion constant

% set initial condition and take FFT
f(x >= 0) = 10 * (abs(x(x >= 0) - 3) < 2); % f for x >= 0
f(x < 0) = -10 * (abs(x(x < 0) + 3) < 2); % f for x < 0

fk = fft(f);  % transform!

% calculate and plot spectral solution at different times using convolution
plot(x, f, 'r')   % IC in red

for t=1:5   % times for the various solutions
    K = exp(-.25/(D*t)*x.^2)/sqrt(4*pi*D*t);     % fundamental soln at time t
    Kk = dx*fft(K);       % scaled FFT of K for convolution
    uk = Kk.*fk;          % implement convolution in wve space
    u = real(fftshift(ifft(uk))); % implement shift, ignore small imaginary bits
    hold on, plot(x, u, 'b'), hold off
end
xlabel('x'), ylabel('u(x,t)')

%%
%%%%%%%%%%%%%%%%%%%%%%%
% Homogeneous Neumann %
%%%%%%%%%%%%%%%%%%%%%%%

% set parameters and define domain
J = 256;      % use a power of two for maximum efficiency
L = 10;       % half width of domain
dx = 2*L/J;   % spatial step size
x = linspace(-L, L-dx, J);  % FFT assumes f(L)=f(-L) so last point unnecessary
D = 1;        % diffusion constant

% set initial condition and take FFT
f(x >= 0) = 10 * (abs(x(x >= 0) - 3) < 2); % f for x >= 0
f(x < 0) = 10 * (abs(x(x < 0) + 3) < 2); % f for x < 0

fk = fft(f);  % transform!

% calculate and plot spectral solution at different times using convolution
plot(x, f, 'r')   % IC in red

for t=1:5   % times for the various solutions
    K = exp(-.25/(D*t)*x.^2)/sqrt(4*pi*D*t);     % fundamental soln at time t
    Kk = dx*fft(K);       % scaled FFT of K for convolution
    uk = Kk.*fk;          % implement convolution in wve space
    u = real(fftshift(ifft(uk))); % implement shift, ignore small imaginary bits
    hold on, plot(x, u, 'b'), hold off
end
xlabel('x'), ylabel('u(x,t)')
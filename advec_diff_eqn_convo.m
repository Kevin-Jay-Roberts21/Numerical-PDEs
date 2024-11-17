% Solving the Advection heat equation using convolution and FFT
% Kevin Roberts
% November 2024

clear all
close all
clc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advection-diffusion eqn %%%%%%%%%%%%%%%%
% CODE FOR u(x,0) = 10*(H(x-1) - H(x-3)) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = 256;      % use a power of two for maximum efficiency
L = 10;       % half width of domain
dx = 2 * L/J;   % spatial step size
x = linspace(-L, L-dx, J);  % FFT assumes f(L)=f(-L) so last point unnecessary
c = 1;        % speed constant
D = 1;        % diffusion constant

% set initial condition and take FFT
f = heaviside(x-1) - heaviside(x-3);     % block wave, width 4, centered at x=3
fk = fft(f);                 % transform!

% set up vector of wave numbers
ks = pi/L * [0:J/2 1-J/2:-1]; % where k = [0:J/2 1-J/2:-1]

% calculate and plot spectral solution at different times
plot(x, f, 'r')   % IC in red

for t=1:5   % times for the various solutions
    K = exp(-.25/(D*t)*(x-c*t).^2)/sqrt(4*pi*D*t);     % fundamental soln at time t
    Kk = dx*fft(K);       % scaled FFT of K for convolution
    uk = Kk.*fk;          % implement convolution in wve space
    u = real(fftshift(ifft(uk))); % implement shift, ignore small imaginary bits
    hold on, plot(x, u, 'b'), hold off
end

xlabel('x'), ylabel('u(x,t)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Model for Cache Valley Radiant Heating %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

J = 256;      % use a power of two for maximum efficiency
L = 5;       % half width of domain
dx = 2 * L/J;   % spatial step size
x = linspace(-L, L-dx, J);  % FFT assumes f(L)=f(-L) so last point unnecessary
hr = 12; % total time in hours
km = 10; % total distance in kilimeters
D = 0.5*km^2/hr;        % diffusion constant

% set initial condition and take FFT
u_0 = zeros(size(x));     
uk = fft(u_0);                 % transform!
uk_p = fft(u_0);

% set up vector of wave numbers
ks = pi/L * [0:J/2 1-J/2:-1]; % where k = [0:J/2 1-J/2:-1]

% calculate and plot spectral solution at different times
plot(x, u_0, 'r')   % IC in red

for t=1:hr   % times for the various solutions (in hours)
    K = exp(-.25/(D*t)*x.^2)/sqrt(4*pi*D*t);     % fundamental soln at time t
    Kk = dx*fft(K);       % scaled FFT of K for convolution
    
    uk = Kk.*uk;          % implement convolution in wve space
    
    g = max(0, sin(pi/12*t).*cos(3*pi/(4*L)*(x + L*(t-6)/12)));
    gk = dx*fft(g); % Fourier transform of the source term
    
    % don't need to add the dtau here since dtau is just 1
    uk_p = 0.5*gk + 0.5*Kk.*uk_p

    % updated solution by adding the effect of g (null + particular)
    uk = uk.*Kk + gk;
    
    u = real(fftshift(ifft(uk))); % implement shift, ignore small imaginary bits
    
    hold on, plot(x, u, 'b'), hold off
end

xlabel('x'), ylabel('u(x,t)')
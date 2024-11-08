% 
% Introductory Matlab script for a solving the heat equation 
%     u_t = D u_xx
% using the fast Fourier Transform
%
% Jim Powell, 4 November 2024
%

% set parameters and define domain
J=256;      % use a power of two for maximum efficiency
L=10;       % half width of domain
dx=2*L/J;   % spatial step size
x=linspace(-L,L-dx,J);  % FFT assumes f(L)=f(-L) so last point unnecessary
D=1;        % diffusion constant

% set initial condition and take FFT
f=10*( abs(x-3)<2);     % block wave, width 4, centered at x=3
fk=fft(f);  % transform!

% set up vector of wave numbers
ks=pi/L*[0:J/2 1-J/2:-1];

% calculate and plot spectral solution at different times

plot(x,f,'r')   % IC in red

for t=1:5   % times for the various solutions
    uk=exp(-D*t*ks.^2).*fk;     % solution in wave space
    u=real( ifft(uk) ); % inverse FFT, ignore small imaginary bits
    hold on, plot(x,u,'b'), hold off
end
xlabel('x'), ylabel('u(x,t)')


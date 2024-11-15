% 
% Introductory Matlab script for a solving the heat equation 
%     u_t = D u_xx
% using the fast Fourier Transform and convolutions
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

% calculate and plot spectral solution at different times using convolution

plot(x,f,'r')   % IC in red

for t=1:5   % times for the various solutions
    K=exp(-.25/(D*t)*x.^2)/sqrt(4*pi*D*t);     % fundamental soln at time t
    Kk=dx*fft(K);      % scaled FFT of K for convolution
    uk=Kk.*fk;          % implement convolution in wve space
    u=real( fftshift( ifft(uk) ) ); % implement shift, ignore small imaginary bits
    hold on, plot(x,u,'b'), hold off
end
xlabel('x'), ylabel('u(x,t)')


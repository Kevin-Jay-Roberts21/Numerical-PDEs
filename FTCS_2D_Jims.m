%
%   FTCS method in two spatial dimensions for heat equation in rectangular
%   domain with Dirichlet boundary conditions.
%
%   written by Jim Powell for Math 6620
%   December 2, 2024
%

%   set up spatial grid and discretization
Lx=5; J=50; dx=Lx/J; x=linspace(0,Lx,J+1);
Ly=10; K=100; dy=Ly/K; y=linspace(0,Ly,K+1);
[X,Y]=meshgrid(x,y);    % this command sets up matrices of coordinate values
D=1;  dt=.24*dx^2/D;    % dt just below stability boundary
T=10; N=ceil(T/dt);     % how long to run simulation
r1=D*dt/dx^2; r2=D*dt/dy^2; % dimensionless rhos

g1=Ly/2-abs(Ly/2-y);    % BC at x=0
g2=sin(pi*y/Ly);        % BC at x=Lx
% assume u=0 at y=0 and y=Ly

up=0*X;     % zero IC
un=0*X;     % just to initialize u^{n+1}
            % Note - since the first and last rows are never updated below,
            %   this automatically sets zero boundary conditions
up(:,1)=g1;     % BC at x=0 (all time-independent)
up(:,J+1)=g2;   % BC at x=Lx
un(:,1)=g1;     % BC at x=0     -- necessary since we will use up=un to update below
un(:,J+1)=g2;   % BC at x=Lx    -- necessary since we will use up=un to update below

for n=1:N   % begin loop in time 

    for j=2:J
        for k=2:K
            un(k,j)=up(k,j)+r1*(up(k,j+1)+up(k,j-1)-2*up(k,j))+...
                r2*(up(k+1,j)+up(k-1,j)-2*up(k,j));
        end     % of k loop
    end     % of j loop
    up=un;  % update for next time step

    if (mod(n*dt,.1)<dt)        % plot every .1 sec
        pcolor(X,Y,up), colormap hot, shading flat, caxis([0 Ly/2]), axis image
            % set color axis so that it maxes at max of BC
            pause(.1)  % pause for rendering the plot
    end     % of plotting if

end         % end of time loop

disp("Code Finished")




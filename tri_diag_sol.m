function [x] = tri_diag_sol(a,b,c,d,b0,bL)

    J = length(a); % could also be length(d)

    % Adjust the first and last entries of d based on boundary conditions
    d(1) = d(1) - b(1) * b0;  % Include left boundary condition
    d(J) = d(J) - c(end) * bL;  % Include right boundary condition
    
    % Initialize alpha and delta
    alpha = zeros(J, 1);
    delta = zeros(J, 1);
    alpha(1) = a(1);
    delta(1) = d(1);

    % defining the entire alpha and delta vector
    for j = 2:J
        alpha(j) = a(j) - b(j-1)*c(j-1)/alpha(j-1);
        delta(j) = d(j) - delta(j-1)*c(j-1)/alpha(j-1);
    end
    
    x = zeros(1, J); % defining the solution
    x(J) = delta(J)/alpha(J); % last x element

    % solving for the x vector (the rest of the elements)
    for j = J-1:-1:1
        x(j) = (delta(j) - b(j)*x(j+1))/alpha(j);
    end
end
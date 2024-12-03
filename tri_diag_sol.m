function [x] = tri_diag_sol(a,b,c,d)

    J = length(a); % could also be length(d)

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
    for j = J:-1:2
        x(j-1) = (delta(j-1) - b(j-1)*x(j))/alpha(j-1);
    end
end
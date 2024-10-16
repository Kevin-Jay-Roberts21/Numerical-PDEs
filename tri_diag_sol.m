function [x] = tri_diag_sol(a,b,c,d)

    J = length(d) % could also be length(a)

    alpha(1) = a(1);
    delta(1) = d(1);

    % defining the entire alpha and delta vector
    for i = 2:J
        alpha(i) = a(i) - b(i-1)*c(i-1)/alpha(i-1);
        delta(i) = d(i) - delta(i-1)*c(i-1)/alpha(i-1);
    end

    x(1) = delta(J)/alpha(J) % first x element

    % solving for the x vector (middle elements)
    for i = 3:J
        x(j-1) = 1/alpha(j-1)*(delta(j-1) - b(j)*alpha(j)) 
    end

    x(J) =  (delta(1) - b(2)*alpha(2))/alpha(1) % last x element
    
end
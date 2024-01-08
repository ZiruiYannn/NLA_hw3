function [lambda, x, y] = DPA(E, A, b, c, s0, tol)
    % DPA: Dominant Pole Approximation
    % Inputs:
    %   E, A - System matrices
    %   b, c - Input and output vectors
    %   s0 - Initial pole estimate
    %   tol - Tolerance for convergence
    % Outputs:
    %   lambda - Approximate dominant pole
    %   x, y - Corresponding eigenpair
    
    % Initialization
    k = 0;
    err = 2*tol;
    s = s0;
    
    while err > tol
        % Solve linear systems
        v = (s * E - A) \ b;
        w = (s * E - A)' \ c;
        
        % Compute new pole estimate
        s_new = s - (c' * v) / (w' * E * v);
        
        % Normalize vectors
        x = v / norm(v);
        y = w / norm(w);
        
        % Compute error
        err = norm(A * x - s_new * E * x, 2);
        
        % Update values for the next iteration
        s = s_new;
        k = k + 1;
    end
    
    % Output the results
    lambda = s;
end

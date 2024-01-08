function [lambda, x, y] = QDPA(M, D, K, b, c, s0, tol)
    % QDPA: Quadratic Dominant Pole Algorithm
    
    k = 0;
    err = inf;
    
    % Initialize variables
    s = s0;

    s2M_minus_K = s^2*M - K;
    sM_plus_sM_D = 2*s*M + D;

    while err > tol
        disp(strcat('Err: ', num2str(err)))
        % Solve linear systems
        MDK = s^2*M + s*D + K;
        v = MDK\b;
        w = MDK'\c;
        
        % Compute new pole estimate
        s_new = (w' * s2M_minus_K * v) / (w' * sM_plus_sM_D * v);
        
        % Normalize vectors
        x = v / norm(v);
        y = w / norm(w);
        
        % Compute error
        s2M_minus_K = s_new^2*M - K;
        sM_plus_sM_D = 2*s_new*M + D;
        err = norm( s2M_minus_K* x - sM_plus_sM_D*x);
        
        % Update variables for the next iteration
        s = s_new;
        k = k + 1;
    end
    
    % Output results
    lambda = s;
end

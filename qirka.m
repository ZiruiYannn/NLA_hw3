function [Mhat, Dhat, Khat, bhat, chat, V, Ehat, Ahat] = irka(M, D, K, b, c, sigma0, tol)
    err = Inf;
    d = numel(sigma0);
    n = size(M,1);
    V = zeros(n,d);
    W = zeros(n,d);
    
    sigma = sigma0;

    while err > tol
        disp(strcat('Err: ', num2str(err)))
        for i=1:d
            MDK = sigma(i)^2*M + sigma(i)*D + K;
            Vi = MDK\ b;
            Wi = MDK' \ c;
            
            Vi = Vi - V(:,1:i-1)*(V(:,1:i-1)'*Vi);
            Vi = Vi - V(:,1:i-1)*(V(:,1:i-1)'*Vi);
            Wi = Wi - W(:,1:i-1)*(W(:,1:i-1)'*Wi);
            Wi = Wi - W(:,1:i-1)*(W(:,1:i-1)'*Wi);

            V(:,i) = Vi / norm(Vi);
            W(:,i) = Wi / norm(Wi);
        end
        Mhat = W'*M*V;
        Dhat = W'*D*V;
        Khat = W'*K*V;
        Ehat = [Mhat zeros(size(Mhat)); Dhat eye(size(Mhat))];
        Ahat = [zeros(size(Mhat)) eye(size(Mhat)); -Khat zeros(size(Mhat))];       
        lambda = sort(eig(Ahat, Ehat));
        diffs = zeros(d,1);
        for i=1:d
            diffs(i) = min(abs(lambda(i) + sigma));
        end
        
        err = norm(diffs) / norm(lambda);
        sigma = -lambda;
    end
    
    bhat = W'*b;
    chat = V'*c;
end

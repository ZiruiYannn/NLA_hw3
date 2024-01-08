function [Mhat, Dhat, Khat, bhat, chat, i] = qgrka(M,D,K, b, c, sigma1, smin, smax, scount, tol)
    V = zeros(size(M,1),0);
    W = zeros(size(M,1),0);
    Mhat = zeros(0);
    Dhat = zeros(0);
    Khat = zeros(0);
    bhat = zeros(0,1);
    
    ss = linspace(smin,smax,scount);
    i = 1;
    while true
        if i == 1
            sigma = sigma1;
        elseif mod(i,2)
            sigma = -sigma;
        else
            MV = M*V(:,1:i-1);
            DV = D*V(:,1:i-1);
            KV = K*V(:,1:i-1);
            ress = arrayfun(@(s) res(s,MV,DV,KV,b,Mhat(1:i-1,1:i-1),Dhat(1:i-1,1:i-1),Khat(1:i-1,1:i-1),bhat(1:i-1)),ss);
            [mx, imax] = max(ress);
            if mx / norm(bhat(1:i-1)) < tol
                chat = V'*c;
                i = i - 1;
                return
            end
            disp(strcat('Res: ', num2str(mx / norm(bhat(1:i-1)))))
            sigma = ss(imax);
        end
        Vi = (sigma^2*M +sigma*D + K) \ b;
        Wi = (sigma^2*M +sigma*D + K)' \ c;
        
        Vi = Vi - V(:,1:i-1)*(V(:,1:i-1)'*Vi);
        Vi = Vi - V(:,1:i-1)*(V(:,1:i-1)'*Vi);
        Wi = Wi - W(:,1:i-1)*(W(:,1:i-1)'*Wi);
        Wi = Wi - W(:,1:i-1)*(W(:,1:i-1)'*Wi);
        V(:,i) = Vi / norm(Vi); W(:,i) = Wi / norm(Wi);
        
        Mhat(i,1:i) = W(:,i)'*M*V(:,1:i);
        Mhat(1:i-1,i) = W(:,1:i-1)'*M*V(:,i);
        Dhat(i,1:i) = W(:,i)'*D*V(:,1:i);
        Dhat(1:i-1,i) = W(:,1:i-1)'*D*V(:,i);
        Khat(i,1:i) = W(:,i)'*K*V(:,1:i);
        Khat(1:i-1,i) = W(:,1:i-1)'*K*V(:,i);
        bhat(i,1) = W(:,i)'*b;
        i = i + 1;
    end
end

function res = res(s,MV,DV,KV,b,Mhat,Dhat,Khat,bhat)
    xhat = (s^2*Mhat +s * Dhat + Khat) \ bhat;
    res = norm((s^2*MV + s * DV + KV)*xhat - b);
end

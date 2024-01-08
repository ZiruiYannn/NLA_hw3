
function [resis_spiral_peec, induc_spiral_peec] = res_ind(E,A,B)
    L_spiral_peec = E;
    R_spiral_peec = -A;
    B_spiral_peec = B;
    
    freq = 10.^[1:0.1:10];
    
    RinvB_spiral_peec = R_spiral_peec\B_spiral_peec;
    RinvL_spiral_peec = R_spiral_peec\L_spiral_peec;
    
    [V_spiral_peec,D_spiral_peec] = eig(RinvL_spiral_peec);
    
    BtV_spiral_peec = transpose(B_spiral_peec)*V_spiral_peec;
    VinvRinvB_spiral_peec = V_spiral_peec\RinvB_spiral_peec;
    
    y_spiral_peec = zeros(size(freq));
    
    for j = 1:size(freq,2)
      y_spiral_peec(1,j) = sum(transpose(BtV_spiral_peec)./(1+sqrt(-1)*2*pi*freq(j)*diag(D_spiral_peec)).*VinvRinvB_spiral_peec);
    end
    
    z_spiral_peec = 1./y_spiral_peec;
    
    resis_spiral_peec = real(z_spiral_peec);
    induc_spiral_peec = imag(z_spiral_peec)./(2*pi*freq);
end
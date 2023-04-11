function [Ovec, hnvec, Pvec] = layerK(Ivec, hvec, theta_vec, K, dt, ltype, flag)
 
    
    
    Ovec = zeros(K,1);
    
    if mod(ltype,2) == 1
        hnvec = zeros(K/2,1);
        Pvec = zeros(K/2,1);
        for i = 1:K/2
        [Ovec(2*i-1), Ovec(2*i), hnvec(i), Pvec(i)] = mzi_cell(Ivec(2*i-1), Ivec(2*i), hvec(i), theta_vec(i), dt, flag(i));      
        end    
    else
        hnvec = zeros(K/2-1,1);
        Pvec = zeros(K/2-1,1);
        Ovec(1) = Ivec(1);
        Ovec(K) = Ivec(K);
        for i = 1:K/2-1
            [Ovec(2*i), Ovec(2*i+1), hnvec(i), Pvec(i)] = mzi_cell(Ivec(2*i), Ivec(2*i+1), hvec(i), theta_vec(i), dt, flag(i));      
        end       
    end
    
end
function [o_, hnvec, Pvec] = RC_KxC(inp_,theta_vec,hvec,dt,K,C)
    
    inpn = inp_;
    
    for i = 1:C 
        %disp(size(hvec((K-1)*(i-1)/2 + 3/2 - mod(i,2)/2:(K-1)*(i)/2 + mod(i,2)/2)));
        [outpn, hnvec((K-1)*(i-1)/2 + 3/2 - mod(i,2)/2:(K-1)*(i)/2 + mod(i,2)/2), Pvec((K-1)*(i-1)/2 + 3/2 - mod(i,2)/2:(K-1)*(i)/2 + mod(i,2)/2)] = layerK(inpn, hvec((K-1)*(i-1)/2 + 3/2 - mod(i,2)/2:(K-1)*(i)/2 + mod(i,2)/2), theta_vec((K-1)*(i-1)/2 + 3/2 - mod(i,2)/2:(K-1)*(i)/2 + mod(i,2)/2), K, dt,i);
        inpn = outpn;
    end
    
    o_ = inpn; 
    
end
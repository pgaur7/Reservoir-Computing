function [o1, o2, hnvec,P] = RC_1coup(i1,theta,hvec,dt)
    
    loss = 1/2;
    % Layer 1 %
    [out1,out2] = node_coup(i1,0,theta(1));
    P1 = (abs(out1))^2;
  
    P = P1;
    [hnvec, phi1] = liq(hvec,P1,dt);
    %[hnvec(2), phi2] = liq(hvec(2),P2,dt);
    o11_1 = out1*exp(1i*phi1)*sqrt(loss);
    o11_2 = out2;
    
    % Layer 2 %
    
    [out1,out2] = node_coup(o11_1,o11_2,theta(2));
    o1 = out1;
    o2 = out2;
   
    
    
end


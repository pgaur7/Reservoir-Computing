function [o1, o2, hn, P] = mzi_cell(i1, i2, h, theta, dt, flag)
    
    loss = 1/2;
    [out1,out2] = node_coup(i1,i2,theta);   
    if flag == 0
        P = (abs(out1))^2;
        [hn, phi] = liq(h,P,dt);
        o1 = out1*exp(1i*phi)*sqrt(loss);
        o2 = out2;
    end
    if flag == 1
        P = (abs(out2))^2;
        [hn, phi] = liq(h,P,dt);
        o1 = out1;
        o2 = out2*exp(1i*phi)*sqrt(loss);
    end
    if flag == -1
        P = 0;
        hn = 0.8;
        o1 = out1;
        o2 = out2;
    end

end

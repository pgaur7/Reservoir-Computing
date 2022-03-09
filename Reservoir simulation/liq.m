function [hn, phi] = liq(h,P,dt)
    c11 = 0.4368;
    c12 = -0.0589;
    c21 = -0.5546;
    c22 = 0.0469;
    hn = h + dt*( (c11*P+c12)*h + (c21*P+c22) );
    phi = 1*2.786*exp(-7.082*hn);
    
%     if hn<0.3
%         hn = 0.3;
%         phi = 1*2.786*exp(-7.082*hn);
%         disp('Error: Height less than 0.3');
%     end
    
    
end

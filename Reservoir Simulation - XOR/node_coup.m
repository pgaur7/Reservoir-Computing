function [out1, out2] = node_coup(in1,in2,theta)
    M = [cos(theta) -1i*sin(theta); -1i*sin(theta) cos(theta)];
    coup = M*[in1;in2];
    m1 = coup(1);
    m2 = coup(2);

    
    out1 = m1;
    out2 = m2;
    
end

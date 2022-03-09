function out = arr_XOR(inp,d)
    N = length(inp);
    out = zeros(N-d,1);
    for i = 1:N-d
        out(i) = xor(inp(i),inp(i+d));
    end
end

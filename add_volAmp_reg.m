function [A, b] = add_volAmp_reg(A, b, b_vol, beta)

if beta ~= 0.0
    b = [b; sqrt(beta) * b_vol.^(-1)];
    A = [A; zeros(length(b_vol), size(A, 2))];
end

end
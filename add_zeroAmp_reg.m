function [A, b] = add_zeroAmp_reg(A, b, b_zero, beta)

if beta ~= 0.0
    b = [b; sqrt(beta) * b_zero];
    A = [A; zeros(length(b_zero), size(A, 2))];
end

end
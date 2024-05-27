function [A, b] = add_zeroAmp_reg(A, b, A_zero, beta)

if beta ~= 0.0
    b = [b; zeros(size(A_zero, 1), 1)];
    A = [A; sqrt(beta) * A_zero];
end

end
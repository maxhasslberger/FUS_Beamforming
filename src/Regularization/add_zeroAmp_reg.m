function [A, b, gamma] = add_zeroAmp_reg(A, b, gamma, A_zero, beta)

if beta ~= 0.0
    b = [b; zeros(size(A_zero, 1), 1)];
    A = [A; sqrt(beta) * A_zero];
    gamma = [gamma; ones(size(A_zero, 1), 1)];
end

end
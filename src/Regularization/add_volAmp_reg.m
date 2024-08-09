function [A, b, gamma] = add_volAmp_reg(A, b, gamma, A_vol, beta)

if beta ~= 0.0
    b = [b; zeros(size(A_vol, 1), 1)];
    A = [A; sqrt(beta)^(-1) * A_vol];
    gamma = [gamma; -1 * ones(size(A_vol, 1), 1)];
end

end
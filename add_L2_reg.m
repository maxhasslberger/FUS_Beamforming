function [A, b, gamma] = add_L2_reg(A, b, gamma, beta)

if beta ~= 0.0
    b = [b; zeros(size(A, 2), 1)];
    A = [A; sqrt(beta) * eye(size(A, 2))];
    gamma = [gamma; ones(size(A, 2), 1)];
end

end
function [A, b] = add_L2_reg(A, b, beta)

if beta ~= 0.0
    b = [b; zeros(size(A, 2), 1)];
    A = [A; sqrt(beta) * eye(size(A, 2))];
end

end
function [A, b] = add_ineq(A, b, n_elements, beta)

if beta ~= 0.0
    b = beta * b;
    A = beta * A;
else
    b = 0.0;
    A = zeros(1, n_elements);
end

end
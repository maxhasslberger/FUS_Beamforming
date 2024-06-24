function [A, b] = add_ineq(A, b, n_elements)

if isempty(A)
    b = 0.0;
    A = zeros(1, n_elements);
end

end
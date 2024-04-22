function [c,ceq] = unitdisk(p, A2_r, A2_i, b2)

% p_r = p(2:ceil(end/2));
% p_i = p(ceil(end/2) + 1:end);
% p_0 = p(1);
% 
% A2_r = A2(1:end/2, :);
% A2_i = A2(end/2 + 1:end, :);
% 
% x_r = p_r ./ sqrt(p_r.^2 + p_i.^2);
% x_i = p_i ./ sqrt(p_r.^2 + p_i.^2);

[x_r, x_i, p_0] = getElements(p);

c = (p_0 * sqrt((A2_r * x_r - A2_i * x_i).^2 + (A2_i * x_r + A2_r * x_i)) - b2);
ceq = [];

end


function [x_r, x_i, p_0] = getElements(p)

p_r = p(2:ceil(end/2));
p_i = p(ceil(end/2) + 1:end);
p_0 = p(1);

x_r = p_r ./ sqrt(p_r.^2 + p_i.^2);
x_i = p_i ./ sqrt(p_r.^2 + p_i.^2);

end

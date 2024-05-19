function p = solvePhases_Amp(A, b, p_init, beta_L2)

% ip.u = max(abs(ip.p)) * ones(size(ip.p)); % universal transducer amplitude
% 
% % Apply constraints (same transducer amplitude among elements)
% ip.A = [ip.A; ip.sq_beta * eye(length(ip.p))];
% b_ip_des = [b_ip_des; ip.sq_beta * ip.u];

% Solve phase retrieval problem
opts = struct;
opts.initMethod = 'custom';
opts.customx0 = p_init;

% Add L2 regularization
[A, b] = add_L2_reg(A, b, beta_L2);

[p, outs, opts] = solvePhaseRetrieval(A, A', b, [], opts);

end
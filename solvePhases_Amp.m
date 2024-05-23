function p = solvePhases_Amp(A, b, p_init, beta_L2)

% Solve phase retrieval problem
opts = struct;
opts.initMethod = 'custom';
opts.customx0 = p_init;

% Add L2 regularization
[A, b] = add_L2_reg(A, b, beta_L2);

[p, outs, opts] = solvePhaseRetrieval(A, A', b, [], opts);

end
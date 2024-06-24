function p = solvePhases_Amp_phasepack(A, b, opt_ids, obs_ids, p_init, init_ids, beta)

% Separate A and b
[A1, ~, b1, ~, A_zero, A_vol, p_init] = prepare_opt_vars(A, b, p_init, opt_ids, obs_ids, init_ids);
clear A;

% Add regularization
[A1, b1] = add_L2_reg(A1, b1, beta(1));
[A1, b1] = add_zeroAmp_reg(A1, b1, A_zero, beta(2));
% [A1, b1] = add_volAmp_reg(A1, b1, A_vol, beta(3));

% Solve phase retrieval problem
opts = struct;
opts.initMethod = 'custom';
opts.customx0 = p_init;
opts.algorithm = 'GerchbergSaxton'; % 'GerchbergSaxton', 'PhaseLift', 'PhaseMax'
opts.verbose = 1;

[p, outs, opts] = solvePhaseRetrieval(A1, A1', b1, [], opts); % var Amplitude

end


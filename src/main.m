clear;
close all;

%% Init
freq_spacing = 10; % kHz
upper_limit = 550; % grid res changes based on this
lower_limit = 450;
f0 = [lower_limit:freq_spacing:upper_limit] * 1e3; % Hz - transducer frequency
% f0 = 500 * 1e3;
n_dim = 2;
dx = 1e-3; % set dx manually since we are using 500 kHz
% dx = [];es
plot_dx_factor = 1;

t1w_filename = fullfile('..', 'Scans', 'dummy_t1w.nii');
ct_filename = fullfile('..', 'Scans', 'dummy_pseudoCT.nii');
% t1w_filename = [];
% ct_filename = [];

sidelobe_tol = 50; % percent of max amplitude

% Simulation config
only_focus_opt = false; % Optimize only for focal spots or entire observation domain
use_greens_fctn = true; % Use Green's function to obtain propagation matrix A (assuming point sources and a lossless homogeneous medium)

% get_current_A = "A_2D_2Trs_70mm_skull"
get_current_A = true; % Use precomputed propagation matrix - can be logical or a string containing the file name in Lin_Prop_Matrices
do_time_reversal = false; % Phase retrieval with time reversal as comparison
do_ground_truth = false; % Ground truth k-wave simulation -> plot_dx_factor
ineq_active = true; % Activate inequality constraints
save_results = false;
get_excitation_vec = false; % Use precomputed excitation vector

if isempty(dx)
    dx_factor = 1;
else
    dx_factor = (1500 / max(f0) / 3) / dx; % = (c0 / f0 / ppw) / dx
end

[kgrid, medium, sensor, sensor_mask, b_des, b_des_pl, b_mask, full_bmask, t_mask_ps, karray_t, only_focus_opt, force_pressures, ...
    active_ids, mask2el, el_per_t, t_pos, t_rot, plot_offset, point_pos, point_pos_m, grid_size, dx_factor1, preplot_arg, domain_ids, input_args] = ...
    init(max(f0), n_dim, dx_factor, ...
    'sidelobe_tol', sidelobe_tol, 't1_scan', t1w_filename, 'ct_scan', ct_filename, 'only_focus_opt', only_focus_opt);

use_greens_fctn = use_greens_fctn & max(medium.sound_speed(:)) == min(medium.sound_speed(:)); % Update green's fctn flag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time Reversal
if do_time_reversal
    tr.p = sim_exe(kgrid, medium, sensor, max(f0), b_des, full_bmask, t_mask_ps, false, input_args);
    % tr.p = max(abs(tr.p)) * exp(-1j * angle(tr.p)); % All elements with same amplitude
    tr.p = conj(tr.p); % Var amplitude
    tr.b = sim_exe(kgridP, mediumP, sensorP, max(f0), tr.p, t_mask_psP, sensor_maskP, true, input_argsP, 'karray_t', karray_tP);
else
    tr = [];
end

%% Obtain propagation operator -> acousticFieldPropagator (Green's functions)
% A = linearPropagator_vs_acousticFieldPropagator(t_mask, f0, medium.sound_speed, kgrid.dx);
if ~exist('A', 'var')
    if ~get_current_A
        % % Old code:
        % A = [];
        % % for f = f0
        % %     A = [A, obtain_linear_propagator(kgrid, medium, sensor, sensor_mask, input_args, t_mask_ps, karray_t, f0, get_current_A, use_greens_fctn, ...
        % %         'active_ids', active_ids)]; %change this bc A's must be diag
        % % end
        % A = [A, obtain_linear_propagator(kgrid, medium, sensor, sensor_mask, input_args, t_mask_ps, karray_t, f0(1), get_current_A, use_greens_fctn, ...
        %         'active_ids', active_ids)];
        % Find A matrix for each frequency
        A_cells = {};
        for f = f0
            if exist(fullfile("..", "Lin_Prop_Matrices","A_" + num2str(f) + "Hz.mat"), "file")
                A_temp = load(fullfile("..", "Lin_Prop_Matrices", "A_" + num2str(f) + "Hz.mat")).A_temp;
                disp('============================================================================')
                disp(['A_', num2str(f),'Hz loaded succesfully!'])
                disp('============================================================================')
            else
                disp('============================================================================')
                disp(['Computing A for ', num2str(f), ' Hz'])                
                disp('============================================================================')
                A_temp = obtain_linear_propagator(kgrid, medium, sensor, sensor_mask, input_args, t_mask_ps, karray_t, f, get_current_A, use_greens_fctn, ...
                    'active_ids', active_ids);
                save(fullfile("..", "Lin_Prop_Matrices", "A_" + num2str(f) + "Hz.mat"), "A_temp");
                disp('============================================================================')
                disp(['A obtained for ', num2str(f), ' Hz and saved to Lin_Prop_Matrices'])
                disp('============================================================================')
            end
            A_cells{end + 1} = A_temp;
        end
        % Block diagonal matrix with A matrices for each frequency
        % A = blkdiag(A_cells{:});
        % save(fullfile("..", "Lin_Prop_Matrices", "A_current.mat"), "A", "-v7.3")
        save(fullfile("..", "Lin_Prop_Matrices", "A_cells_current.mat"), "A_cells", "-v7.3")
        disp("============================================================================")
        disp("Obtained full A matrix and saved to Lin_Prop_Matrices as A_cells_current.mat")
        disp("============================================================================")
    else
        % disp("Loading precomputed Propagation Matrix...")
        % A = load(fullfile("..", "Lin_Prop_Matrices", "A_current.mat")).A;
        % disp("Propagation Matrix loaded successfully!")
                       
        A_cells = load(fullfile("..", "Lin_Prop_Matrices", "A_cells_current.mat")).A_cells;                            
        disp("A_cells (contains individual propagation matrices for each frequency) loaded successfully!")
    end
end

% if kgrid.dim == 3
% A = ones(numel(b_mask), 1);
% end

%% Solve Inverse Problem
tic

obs_ids = reshape(logical(full_bmask), numel(full_bmask), 1); % indices of target region

if only_focus_opt
    domain_ids = 1; % To be discarded -> off-target pressures could be too large!
    skull_ids = 1;

    % Preserve sonicated points only
    ip.A = A(obs_ids, :); % Question --------> What is ip?
    % -> ip stands for inverse problem, the name of the method we're using
    % to compute array phases (and amplitudes). We used to compare it with
    % time reversal (tr) which is another method with different properties.
    % The struct variable ip contains all the data associated with the
    % method (ip.A -> propagation matrix; ip.p -> excitation vector that
    % contains all transducer phases and amplitudes (as complex numbers);
    % ip.b -> solution vector = acoustic pressure distribution in the
    % entire observation domain; ip.*_gt -> ground truth variables to
    % confirm the results in a higher spatial resolution. If ground_truth =
    % false we use the _gt variables to compare different inverse problems
    % against each other. But for you solvePhasesOnly is not really
    % relevant. Same for only_focus_opt, just assume it to be false for
    % now.

    b_ip_des = b_des; % = b_des_pl(obs_ids)

    vol_ids = true(size(ip.A, 1), 1);
    init_ids = vol_ids;
    ip.beta = 0.0;
else
%     [domain_ids, skull_ids] = limit_space(medium.sound_speed); % Indices considered in optimization (intracranial and skull)
    skullMask = medium.sound_speed > min(medium.sound_speed(:));
    skull_ids = reshape(skullMask, [], 1);

    % Take entire observation grid into account
    ip.A = A_cells; % Use cell array of A matrices rather than large block diagonal A matrix
    % Question ---> I am still a little confused about what information the A matrix
    % contains. I understand that this is the propagation matrix and that
    % this is multiplied with the excitation vector to get the target
    % pressures, but I am not sure about the contents of the A matrix.

    % Answer: The contents are a bunch of complex values z * exp(j * phi). 
    % Each one describes A) the attenuation (z) and B) the phase shift 
    % (phi) from ONE transducer element to ONE point in the observation 
    % domain. Pls see my example in solvePhasesAmp and refer to my slides 
    % where I describe how to obtain A. Let me know if it's clear.
    b_ip_des = b_des_pl;

    vol_ids = obs_ids; % Indices that correspond to the target volume(s)
    [init_ids, ~, b_mask_plot] = get_init_ids(kgrid, min(medium.sound_speed(:)) / mean(f0), b_mask, force_pressures); % Indices where pressure values given
    
    ip.beta = 0.0;

    % New preplot with init point arg
    preplot_arg = preplot_arg + b_mask_plot;
    preplot_arg = preplot_arg / max(preplot_arg(:));
    plot_results(kgrid, [], preplot_arg, 'Plot Preview 2', mask2el, t1w_filename, plot_offset, grid_size, dx_factor1, false, [], 'slice', point_pos.slice, ...
        'colorbar', false, 'cmap', hot());
end

% Compute pseudoinverse for each frequency, then stack the resulting
% excitation vectors (each column is one excitation vector) to obtain the initial solution
if ~get_excitation_vec
    m = size(A_cells, 2);
    p_init = [];
    
    for i = 1:m
        A_temp = A_cells{i};
        p_init = [p_init, pinv(A_temp(init_ids, :)) * b_ip_des(init_ids, :)];
    end       

% % ---------------------------------------------
% % >>>>>>>>>>>> Test Cost Function <<<<<<<<<<<<<
% [val,num_failed_constraints] = testCostFctn(ip.A, b_ip_des, domain_ids, skull_ids, vol_ids, p_init, init_ids, ip.beta, f0);
% val % 4.2264e+05
% num_failed_constraints % 3512
% 
% % For initial solution, cost func value is 1.0e+05 * [3.0000, -3.0000], with a norm of 4.2426e+05
% % 0 indices failed the inequality constraint. Question --> Is this surprising? I thought it was because we are
% % using an initial excitation vector that simply combines 2 initial
% % solutions, which are each obtained by solving the one frequency problem
% % for each frequency
% 
% % ---------------------------------------------

% Obtain optimal p for multiple frequencies

    ip.p = solvePhasesAmpMultiFreq(ip.A, b_ip_des, domain_ids, skull_ids, vol_ids, p_init, init_ids, ip.beta, f0);
else
    ip.p = load(fullfile("..", "Excitation_Vecs", "p_470khz_5f.mat")).p; % Change excitation vectors file name manually
end
% % Single freq:
% p_init = pinv(ip.A(init_ids, :)) * b_ip_des(init_ids);
% ip.p = solvePhasesAmp(ip.A, b_ip_des, domain_ids, skull_ids, vol_ids, p_init, init_ids, ip.beta, ineq_active); % var Amp

% ip.p_gt = solvePhasesAmpMultiFreq(ip.A, b_ip_des, domain_ids, skull_ids, vol_ids, p_init, init_ids, ip.beta); % var Amp
% ip.p = p_init;

ip.t_solve = toc;

%% Obtain Acoustic profile
nfreq = size(A_cells,2);
% Single freq
% ip.b = A * ip.p;
% Multifreq
ip.b = timeDomainSum(f0, A_cells, ip.p);
% plot_focalpt_waveform(f0,A_cells,ip.p,init_ids);

disp("pressures at the foci:")
disp(double(ip.b(domain_ids &  init_ids)))

ip.b = reshape(ip.b, size(kgrid.k));

% Comment this in for intracranial pressures
% ip.b(~domain_ids) = 0.0;
% ip.b(~domain_ids & ~skull_ids) = 0.0;

if do_ground_truth % For different resolution: Only supported in 3D at the moment
    [kgridP, mediumP, sensorP, sensor_maskP, ~, ~, ~, ~, t_mask_psP, karray_tP, ~, ~, ...
    ~, ~, ~, ~, ~, plot_offsetP, point_posP, ~, grid_sizeP, dx_factorP, ~, ~, input_argsP] = ...
    init(max(f0), n_dim, dx_factor * plot_dx_factor, 't1_scan', t1w_filename, 'ct_scan', ct_filename, 'only_focus_opt', only_focus_opt);

    % Need to fix handling of f0. Might need to adjust for using multiple frequencies   
    if use_greens_fctn
        [amp_in, phase_in] = get_amp_phase_mask(kgrid, f0, ip.p, t_mask_psP, karray_tP);
        ip.b_gt = acousticFieldPropagator(amp_in, phase_in, kgrid.dx, f0, medium.sound_speed);
    else
        ip.b_gt = sim_exe(kgridP, mediumP, sensorP, f0, ip.p, t_mask_psP, sensor_maskP, true, input_argsP, 'karray_t', karray_tP);
        ip.b_gt = reshape(ip.b_gt, size(kgridP.k));
    end

    % [domain_ids_gt, skull_ids_gt] = limit_space(mediumP.sound_speed);
    % ip.b_gt(~domain_ids_gt) = 0.0;
else
    % ip.p_gt = solvePhasesOnly(ip.A, b_ip_des, domain_ids, skull_ids, vol_ids, p_init, init_ids, ip.beta, mask2el, el_per_t, true); % Amp fixed
%     ip.p_gt = solvePhases_Amp_phasepack(ip.A, b_ip_des, domain_ids, vol_ids, p_init, init_ids, beta); % var Amp
%     ip.p_gt = p_init;
    % ip.b_gt = A * ip.p_gt;
    % ip.b_gt = reshape(ip.b_gt, size(kgrid.k));

    % ip.b_gt(~domain_ids) = 0.0;

    % For multifreq without ground truth, we set ground truth values to same values as inverse problem for the sake of simplicity so that
    % comp plots don't error out
    ip.p_gt = ip.p;
    ip.b_gt = ip.b;
end

plot_thr = min(b_ip_des) +  1e3;

%% Save Results in mat-file
current_datetime = string(datestr(now, 'yyyymmddHHMMSS'));
if save_results
    res_filename = "results";
    if ~only_focus_opt
        ip.A = []; % A might be very large...
    end
    save(fullfile("..", "Results", current_datetime + "_" + res_filename + ".mat"), ...
        "f0", "kgrid", "b_mask", "t_mask_ps", "active_ids", "init_ids", "mask2el", "t_pos", "t_rot", "tr", "ip", "b_ip_des", "point_pos", "point_pos_m", ...
        "only_focus_opt", "input_args");
end

%% TR Results
if do_time_reversal
    % tr.b(~domain_ids) = 0.0;
    plot_results(kgrid, tr.p, tr.b, 'Time Reversal', mask2el, t1w_filename, plot_offset, grid_size, dx_factor1, save_results, current_datetime, 'slice', point_pos.slice);
    plot_results(kgrid, [], abs(tr.b) > plot_thr, 'Time Reversal', mask2el, t1w_filename, plot_offset, grid_size, dx_factor1, save_results, current_datetime, ...
        'slice', point_pos.slice, 'cmap', gray());
end

%% IP Results
plot_results(kgrid, ip.p, ip.b, 'Inverse Problem', mask2el, t1w_filename, plot_offset, grid_size, dx_factor1, save_results, current_datetime, 'slice', point_pos.slice);
plot_results(kgrid, [], abs(ip.b) > plot_thr, 'Inverse Problem', mask2el, t1w_filename, plot_offset, grid_size, dx_factor1, save_results, current_datetime, ...
    'slice', point_pos.slice, 'colorbar', true, 'cmap', gray()); % plot mask with pressure above off-target limit

if do_ground_truth % For different resolution: Only supported in 3D at the moment
    plot_results(kgridP, [], ip.b_gt, 'Ground Truth', mask2el, t1w_filename, plot_offsetP, grid_sizeP, dx_factorP, save_results, ...
        current_datetime, 'slice', point_posP.slice);
    plot_results(kgridP, [], abs(ip.b_gt) > plot_thr, 'Ground Truth', mask2el, t1w_filename, plot_offsetP, grid_sizeP, dx_factorP, save_results, ...
        current_datetime, 'slice', point_posP.slice, 'colorbar', true, 'cmap', gray());
    
%     if plot_dx_factor == 1
%         err = abs(ip.b) - abs(ip.b_gt);
%         plot_results(kgrid, [], err, 'Difference', mask2el, t1w_filename, plot_offset, grid_size, dx_factor1, save_results, current_datetime, 'slice', point_pos.slice);
%         figure
%         histogram(err(:))
%         xlabel("Pressure Deviation (Pa)")
%     end
% else
    % plot_results(kgrid, ip.p_gt, ip.b_gt, 'Inverse Problem Comp', mask2el, t1w_filename, plot_offset, grid_size, dx_factor1, save_results, ...
        % current_datetime, 'slice', point_pos.slice);
    % plot_results(kgrid, [], abs(ip.b_gt) > max(b_des) / 2, 'Inverse Problem Comp', mask2el, t1w_filename, plot_offset, grid_size, dx_factor1, save_results, ...
        % current_datetime, 'slice', point_pos.slice, 'colorbar', true, 'cmap', gray());
end

%% Metrics evaluation
disp("Time until solver converged: " + string(ip.t_solve) + " s")

real_ip = abs(reshape(ip.b, [], 1));
real_ip_gt = abs(reshape(ip.b_gt, [], 1));

if only_focus_opt
    real_ip = real_ip(obs_ids);
    real_ip_gt = real_ip_gt(obs_ids);
    
    err_ip = real_ip - b_ip_des;
    err_ip_gt = real_ip_gt - b_ip_des;
    
    std_ip = sqrt(sum(err_ip.^2) / length(err_ip));
    std_ip_gt = sqrt(sum(err_ip_gt.^2) / length(err_ip_gt));
    
    maxDev_ip = max(abs(err_ip));
    maxDev_ip_gt = max(abs(err_ip_gt));
    
    % fprintf("\nDesired Amplitudes (kPa):\n")
    % disp(b_ip_des' * 1e-3)
    % 
    % fprintf("\nReal Amplitudes (kPa):\n")
    % disp(real_ip' * 1e-3)
    % disp(real_ip_gt' * 1e-3)
    
    fprintf("\nInverse Problem STD (kPa):\n")
    disp(std_ip * 1e-3)
    disp(std_ip_gt * 1e-3)
    
    fprintf("\nInverse Problem Max Dev (kPa):\n")
    disp(maxDev_ip * 1e-3)
    disp(maxDev_ip_gt * 1e-3)
else
    offTar_real_ip = real_ip;
    offTar_real_ip(vol_ids | ~domain_ids) = [];
    offTar_real_ip_gt = real_ip_gt;
    offTar_real_ip_gt(vol_ids | ~domain_ids) = [];

    tar_real_ip = real_ip(vol_ids);
    tar_real_ip_gt = real_ip_gt(vol_ids);

    init_real_ip = real_ip(init_ids);
    init_real_ip_gt = real_ip_gt(init_ids);

    fprintf("\nInverse Problem Init Points (kPa):\n")
    disp(init_real_ip' * 1e-3)
    disp(init_real_ip_gt' * 1e-3)

    fprintf("\nInverse Problem Min Target Area (kPa):\n")
    disp(min(tar_real_ip) * 1e-3)
    disp(min(tar_real_ip_gt) * 1e-3)

    fprintf("\nInverse Problem Max Off-Target (kPa):\n")
    disp(max(offTar_real_ip) * 1e-3)
    disp(max(offTar_real_ip_gt) * 1e-3)
end

% point_coord = [-18, -27];
% point_val_ip = abs(ip.b(plot_offset(1) + point_coord(1), plot_offset(3) + point_coord(2)));
% point_val_ip_gt = abs(ip.b_gt(plot_offset(1) + point_coord(1), plot_offset(3) + point_coord(2)));
% 
% fprintf("\nInverse Problem sel. Point Pressure (kPa):\n")
% disp(point_val_ip * 1e-3)
% disp(point_val_ip_gt * 1e-3)

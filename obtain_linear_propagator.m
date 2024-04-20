function A = obtain_linear_propagator(kgrid, medium, sensor, sensor_mask, input_args, t_mask, karray_t, f0, get_current_A, use_greens_fctn, varargin)

active_ids = [];

if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'active_ids'
                active_ids = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

if islogical(get_current_A)

    if ~get_current_A
    
%         use_greens_fctn = use_greens_fctn & max(medium.sound_speed(:)) == min(medium.sound_speed);
        if use_greens_fctn
        
            el_ids = find(t_mask);
%             el2mask_ids = sort(el2mask_ids); % Mapping mask -> element index
%             el_ids = el_ids(el2mask_ids);

            phase_in = 0;
            input = 1.0;
            A = single(zeros(kgrid.total_grid_points, numel(el_ids)));

            input = zeros(length(el_ids), 1);
            input(1) = 1;

            % Excite one element at a time and obtain one column (observation) after the other
            for i = 1:length(el_ids)
                disp("Offline Simulation " + string(i) + "/" + string(size(A, 2)))

                amp_in = zeros(size(t_mask));
                amp_in(el_ids(i)) = 1;

%                 a_coli = single(acousticFieldPropagator(amp_in, phase_in, kgrid.dx, f0, medium.sound_speed));
                a_coli = single(sim_exe(kgrid, medium, sensor, f0, input, t_mask, sensor_mask, true, input_args)); % TODO: Test and confirm
                a_coli = reshape(a_coli, [], 1);
                A(:, i) = a_coli;

                input = circshift(input, 1); % Move unitary excitation by 1 element
            end
        
        else
            A = single(zeros(kgrid.total_grid_points, karray_t.number_elements));

            input = 1.0;
            karray_tmp = kWaveArray();

            % Excite one element at a time and obtain one column (observation) after the other
            for i = 1:karray_t.number_elements
                disp("Offline Simulation " + string(i) + "/" + string(size(A, 2)))

                el_x = karray_t.elements{i};
                karray_tmp.addRectElement(el_x.position, el_x.length, el_x.width, el_x.orientation);

                a_coli = single(sim_exe(kgrid, medium, sensor, f0, input, [], sensor_mask, true, input_args, 'karray_t', karray_tmp));
                a_coli = reshape(a_coli, [], 1);
                A(:, i) = a_coli;

                karray_tmp.removeElement(1);
%                 input = circshift(input, 1); % Move unitary excitation by 1 element
            end
        end
        
        save(fullfile("Lin_Prop_Matrices", "A_current.mat"), "A", "-v7.3")
        
    else
        disp("Loading precomputed Propagation Matrix...")
        A = load(fullfile("Lin_Prop_Matrices", "A_current.mat")).A;
        disp("Propagation Matrix loaded successfully!")
    end

else
    disp("Loading precomputed Propagation Matrix...")
    A = load(fullfile("Lin_Prop_Matrices", string(get_current_A) + ".mat")).A;
    disp("Propagation Matrix loaded successfully!")
end

% Only use subset of available transducers
if ~isempty(active_ids)
    if length(active_ids) < size(A, 2)
        
        % Update A
        A = A(:, active_ids);
    end
end

end
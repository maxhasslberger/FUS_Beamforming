function A = obtain_linear_propagator(kgrid, medium, sensor, sensor_mask, input_args, t_mask, karray_t, f0, get_current_A, use_greens_fctn, varargin)

active_ids = [];
save_iter = 10;
std_filename = "A_current";
tmp_filename = "A_current_tmp";

if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'active_ids'
                active_ids = varargin{arg_idx+1};
            case 'save_iter'
                save_iter = varargin{arg_idx+1};
            case 'std_filename'
                std_filename = varargin{arg_idx+1};
            case 'tmp_filename'
                tmp_filename = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

if islogical(get_current_A)

    if ~get_current_A
    
        if exist(fullfile("..", "Lin_Prop_Matrices", tmp_filename + ".mat"), "file") % Start with tmp matrix if prior run interrupted
            disp("Loading precomputed Tmp Matrix...")
            A = load(fullfile("..", "Lin_Prop_Matrices", tmp_filename + ".mat")).A;
            disp("Tmp Matrix loaded successfully!")

            i_start = find(all(A==0), 1);
        else
            i_start = 1;
        end
    
        phase_in = 0;

        if isempty(karray_t) % Use point sources as transducers
        
            el_ids = find(t_mask);
%             el2mask_ids = sort(el2mask_ids); % Mapping mask -> element index
%             el_ids = el_ids(el2mask_ids);

            if ~exist("A", "var")
                A = zeros(kgrid.total_grid_points, numel(el_ids), 'single');
            end

            input = zeros(length(el_ids), 1);
            input(i_start) = 1;

            % Excite one element at a time and obtain one column (observation) after the other
            for i = i_start:length(el_ids)
                disp("Offline Simulation " + string(i) + "/" + string(size(A, 2)))
                
                if use_greens_fctn
                    amp_in = zeros(size(t_mask));
                    amp_in(el_ids(i)) = 1;
                    a_coli = single(acousticFieldPropagator(amp_in, phase_in, kgrid.dx, f0, medium.sound_speed));
                else
                    a_coli = single(sim_exe(kgrid, medium, sensor, f0, input, t_mask, sensor_mask, true, input_args));
                    input = circshift(input, 1); % Move unitary excitation by 1 element
                end

%                 a_col_plot = a_coli;
%                 max_val = 0.25 * max(abs(a_col_plot(:)));
%                 a_col_plot(abs(a_col_plot) > max_val) = max_val * exp(1j * angle(a_col_plot(abs(a_col_plot) > max_val)));
%                 a_col_plot = abs(a_col_plot) .* cos(angle(a_col_plot));
%                 plot_results(kgrid, [], a_col_plot, '', [], '..\Scans\dummy_t1w.nii', [97 128 127], [0.1920 0.2560], 1, false, []);

                a_coli = reshape(a_coli, [], 1);
                A(:, i) = a_coli;

                if mod(i, save_iter) == 0
                    save(fullfile("..", "Lin_Prop_Matrices", tmp_filename + ".mat"), "A", "-v7.3")
                end
            end
        
        else % Use karray class
            if ~exist("A", "var")
                A = zeros(kgrid.total_grid_points, karray_t.number_elements, 'single');
            end

            input = 1.0;
            karray_tmp = kWaveArray();

            % Excite one element at a time and obtain one column (observation) after the other
            for i = i_start:karray_t.number_elements
                disp("Offline Simulation " + string(i) + "/" + string(size(A, 2)))

                el_x = karray_t.elements{i};

                if strcmp(el_x.type, 'disc')
                    karray_tmp.addDiscElement(el_x.position, el_x.diameter, el_x.focus_position);
                else
                    karray_tmp.addRectElement(el_x.position, el_x.length, el_x.width, el_x.orientation);
                end

                if use_greens_fctn
                    amp_in = karray_tmp.getArrayGridWeights(kgrid);
                    [amp_in, phase_in] = get_amp_phase_mask(kgrid, f0, input, [], karray_tmp);
                    a_coli = single(acousticFieldPropagator(amp_in, phase_in, kgrid.dx, f0, medium.sound_speed));
                else
                    a_coli = single(sim_exe(kgrid, medium, sensor, f0, input, [], sensor_mask, true, input_args, 'karray_t', karray_tmp));
                end

%                 a_col_plot = a_coli;
%                 max_val = 0.25 * max(abs(a_col_plot(:)));
%                 a_col_plot(abs(a_col_plot) > max_val) = max_val * exp(1j * angle(a_col_plot(abs(a_col_plot) > max_val)));
%                 a_col_plot = abs(a_col_plot) .* cos(angle(a_col_plot));
%                 plot_results(kgrid, [], a_col_plot, '', [], '..\Scans\dummy_t1w.nii', [97 128 127], [0.1920 0.2560 0.2560], 1, false, []);

                a_coli = reshape(a_coli, [], 1);
                A(:, i) = a_coli;

                karray_tmp.removeElement(1);
                if mod(i, save_iter) == 0
                    save(fullfile("..", "Lin_Prop_Matrices", tmp_filename + ".mat"), "A", "-v7.3")
                end
            end
        end
        
        delete(fullfile("..", "Lin_Prop_Matrices", tmp_filename + ".mat")) % delete tmp file
        save(fullfile("..", "Lin_Prop_Matrices", std_filename + ".mat"), "A", "-v7.3")
        
    else
        disp("Loading precomputed Propagation Matrix...")
        A = load(fullfile("..", "Lin_Prop_Matrices", std_filename + ".mat")).A;
        disp("Propagation Matrix loaded successfully!")
    end

else
    disp("Loading precomputed Propagation Matrix...")
    A = load(fullfile("..", "Lin_Prop_Matrices", string(get_current_A) + ".mat")).A;
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
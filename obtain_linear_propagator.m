function A = obtain_linear_propagator(t_mask, f0, sound_speed, dx, set_current_A, varargin)

mask2el_ids = [];
active_tr_ids = [];

if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'mask2el_ids'
                mask2el_ids = varargin{arg_idx+1};
            case 'tr_ids'
                active_tr_ids = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

if islogical(set_current_A)

    if ~set_current_A
    
        if min(sound_speed(:)) == max(sound_speed(:)) % homogeneous medium
        
            el_ids = find(t_mask);
%             el2mask_ids = sort(el2mask_ids); % Mapping mask -> element index
%             el_ids = el_ids(el2mask_ids);

            phase_in = 0;
            A = single(zeros(numel(t_mask), numel(el_ids)));

            % Excite one element at a time and obtain one column (observation) after the other
            for i = 1:length(el_ids)
                amp_in = zeros(size(t_mask));
                amp_in(el_ids(i)) = 1;
                a_coli = single(acousticFieldPropagator(amp_in, phase_in, dx, f0, sound_speed));
                a_coli = reshape(a_coli, [], 1);
                A(:, i) = a_coli;
            end
        
        else
            error("Only linear media supported at the moment!")
        end
        
        save(fullfile("Lin_Prop_Matrices", "A_current.mat"), "A", "-v7.3")
        
    else
        A = load(fullfile("Lin_Prop_Matrices", "A_current.mat")).A;
    end

else
    A = load(fullfile("Lin_Prop_Matrices", string(set_current_A) + ".mat")).A;
end

% Only use subset of available transducers
if ~isempty(mask2el_ids) && ~isempty(active_tr_ids)
    if sum(t_mask(:)) < size(A, 2)

        % Obtain transducer element ids
        mask2el_ids = mask2el_ids(:, active_tr_ids);
        mask2el_ids = sort(mask2el_ids(:));
        
        % Update A
        A = A(:, mask2el_ids);
    end
end

end
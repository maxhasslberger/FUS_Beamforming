function A = obtain_linear_propagator(t_mask, b_mask, f0, sound_speed, dx, only_focus_opt, set_current_A)

if islogical(set_current_A)

    if ~set_current_A
    
        if min(sound_speed(:)) == max(sound_speed(:)) % homogeneous medium
        
            % Excite one element after the other
            phase_in = 0;
            el_ids = find(t_mask);
            
            if only_focus_opt
                obs_ids = find(b_mask);
            else
                obs_ids = 1:numel(t_mask);
            end
            
            A = zeros(numel(obs_ids), numel(el_ids));
            for i = 1:length(el_ids)
                amp_in = zeros(size(t_mask));
                amp_in(el_ids(i)) = 1;
                a_coli = acousticFieldPropagator(amp_in, phase_in, dx, f0, sound_speed);
                a_coli = reshape(a_coli, [], 1);
                A(:, i) = a_coli(obs_ids);
            end
        
        else
            error("Only linear media supported at the moment!")
        end
        
        save(fullfile("Lin_Prop_Matrices\A_current.mat"), "A")
        
    else
        A = load(fullfile("Lin_Prop_Matrices\A_current.mat")).A;
    end

else
    A = load(fullfile("Lin_Prop_Matrices\" + string(set_current_A) + ".mat")).A;
end

% figure;
% imagesc(abs(A))
% title("A")
% colorbar;

end
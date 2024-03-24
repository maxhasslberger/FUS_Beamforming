function A = linearPropagator_vs_acousticFieldPropagator(t_mask, f0, sound_speed, dx) % for now entire plane instead of only b_mask

if min(sound_speed(:)) == max(sound_speed(:))
    amp_in = t_mask * 1.0;
    phase_in = 0;
    [all_elements_active, ~] = acousticFieldPropagator(amp_in, phase_in, dx, f0, sound_speed);

    % Excite one element after the other
    el_ids = find(t_mask);
    
    A = zeros(numel(t_mask), numel(el_ids));
    for i = 1:length(el_ids)
        amp_in = zeros(size(t_mask));
        amp_in(el_ids(i)) = 1;
        a_coli = acousticFieldPropagator(amp_in, phase_in, dx, f0, sound_speed);
        A(:, i) = reshape(a_coli, [], 1);
    end

else
    error("Only linear media supported at the moment!")
end

all_elements_active_2 = abs(A * ones(numel(el_ids), 1));
all_elements_active_2 = reshape(all_elements_active_2, size(t_mask, 1), size(t_mask, 2));

figure;
imagesc(all_elements_active_2)
title("Linear propagator solution")
colorbar;

figure;
imagesc(all_elements_active)
title("Exact Green's function solution")
colorbar;

figure;
imagesc(abs(all_elements_active_2 - all_elements_active))
title("Difference")
colorbar;

figure;
imagesc(abs(A))
title("Linear propagator entries")
colorbar;

end
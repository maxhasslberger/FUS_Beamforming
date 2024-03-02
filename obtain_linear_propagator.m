function A = obtain_linear_propagator(t_mask, b_mask, f0, sound_speed, dx) % for now entire plane instead of only b_mask

if min(sound_speed(:)) == max(sound_speed(:))
    amp_in = t_mask * 1.0;
    phase_in = 0;
    [all_elements_active, ~] = acousticFieldPropagator(amp_in, phase_in, dx, f0, sound_speed);

    % Excite one element after the other
    el_ids = find(amp_in);
    
    A = zeros(numel(t_mask), numel(el_ids));
    for i = length(el_ids)
        amp_in = zeros(size(t_mask));
        amp_in(el_ids(i)) = 1;
        [A(:, i), ~] = acousticFieldPropagator(amp_in, phase_in, dx, f0, sound_speed);
    end

else
    error("Only linear media supported right now!")
end

all_elements_active_2 = A * ones(numel(el_ids), 1);
all_elements_active_2 = reshape(all_elements_active_2, size(t_mask, 1), size(t_mask, 2));

figure;
imagesc(all_elements_active_2)
title("Approx")

figure;
imagesc(all_elements_active)
title("Exact")

figure;
imagesc(abs(all_elements_active_2 - all_elements_active))
title("Diff")

end
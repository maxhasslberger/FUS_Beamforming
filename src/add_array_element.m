function karray_t = add_array_element(karray_t, el_pos, el_sz, el_rot)

if length(el_sz) == 1 % Only diameter -> disc elements
    % Convert rotation into vector on transducer axis
    R = rotation_matrixXYZ(el_rot);
    tr_axis = el_pos + R * [0; 0; 1] * 1e-3; % Make sure that transducer is not placed too close to the grid boundaries
    
    karray_t.addDiscElement(el_pos, el_sz, tr_axis);
else
    karray_t.addRectElement(el_pos, el_sz(1), el_sz(2), el_rot);
end

end


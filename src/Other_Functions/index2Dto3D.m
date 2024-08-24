function mat_2D = index2Dto3D(mat_3D, dim, idx)

dim = dim2num(dim);

if dim == 1
    mat_2D = squeeze(mat_3D(idx, :, :));
elseif dim == 2
    mat_2D = squeeze(mat_3D(:, idx, :));
else
    mat_2D = squeeze(mat_3D(:, :, idx));
end

end


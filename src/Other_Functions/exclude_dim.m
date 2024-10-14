function rem_dims = exclude_dim(dim)

dim = dim2num(dim);

rem_dims = [1, 2, 3];
rem_dims(rem_dims == dim) = [];

end
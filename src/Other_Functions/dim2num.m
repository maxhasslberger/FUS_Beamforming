function dim = dim2num(dim)

if ischar(dim)
    dim = double(dim) - double('X') + 1;
end

end
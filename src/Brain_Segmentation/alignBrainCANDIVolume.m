function alignData = alignBrainCANDIVolume(data,aff,affRef)
% Aligns a volume to a reference orientation (axis and direction) 
% specified by an affine matrix

% Copyright 2022 The MathWorks, Inc.

dim = ndims(data);
rasAff = getRAS(aff,dim);
rasRef =  getRAS(affRef,dim);

if ~issorted(rasAff)
    alignAx = rasAff;
else 
    alignAx = rasRef;
end

% Align axes
aff(:,[rasRef,size(aff,1)]) = aff(:,[rasAff,size(aff,1)]);
alignData= permute(data,alignAx);
dot_products = sum(aff(1:dim,1:dim).*affRef(1:dim,1:dim));
for i = 1:dim
    if dot_products(i) <0
        alignData = flip(alignData,i);
    end
end
end


function affRAS =  getRAS(aff,dim)
affInv = inv(aff);
[~,affRAS] = max(abs(affInv(1:dim,1:dim)));
end
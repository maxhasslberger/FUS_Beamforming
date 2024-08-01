function [B,cropIdx,volSize]  = preProcessBrainCANDIData(data,metaData,cropSize,resample)
% Pre-process the data before passing it to the pre-trained model.

% Copyright 2022 The MathWorks, Inc.

% Perform resampling to make data isotropic
if resample
    aff = metaData.Transform.T';
    voxelSize = metaData.PixelDimensions;
    newVoxelSize = [1 1 1];
    [data,aff]  = resampleVolume(data,aff,voxelSize,newVoxelSize); %#ok<*UNRCH> 
else
    aff = metaData.Transform.T';
end

% Align volume to ref
affRef = eye(4);
data = alignBrainCANDIVolume(data,aff,affRef);

% Cropping
if isscalar(cropSize)
    cropSize= repmat(cropSize,[1 3]);
end
volSize = size(data);

% Check that padded image size divisible by 2^numLevels if not find closest
% number divisable by N
numLevels = 5;
N = 2^numLevels;
if ~rem(cropSize,N)==0
    cropSize = floor(cropSize/N)*N;
end

% Center the data
minCropIdx = max((volSize-cropSize)/2,[0 0 0]);
maxCropIdx = min(minCropIdx+cropSize,volSize);
cropIdx = [minCropIdx maxCropIdx];
B = data(cropIdx(1)+1:cropIdx(4), cropIdx(2)+1:cropIdx(5), cropIdx(3)+1:cropIdx(6));

% Normalize
inMax = prctile(B(:),99.5);
inMin = prctile(B(:),0.5);
B = rescale(B,0,1,InputMax=inMax,InputMin=inMin);
end

function [resampleData,aff] = resampleVolume(data,aff,voxelSize,newVoxelSize)
% Resample data to make it isotropic of voxel size 1mmx1mmx1mm.

% Copyright 2022 The MathWorks, Inc.

% Make isotropic volume with new voxel size
scale = voxelSize./newVoxelSize;
outputSz = [ceil(size(data,1)*scale(1)),ceil(size(data,2)*scale(2)),...
    ceil(size(data,3)*scale(3))];
resampleData = imresize3(data,outputSz,method="linear",AntiAliasing=false);

% Update affine transform
aff(1:3,1:3)= aff(1:3,1:3)./scale;
aff(1:3,end) = aff(1:3,end)-(aff(1:3,1:3)*0.5*(scale'-1));
end
function seg = postProcessBrainCANDIData(postIn,postFlipIn,imSize,cropIdx,metaData,classNames,labelIDs)
% Post process the data to get the segmentation maps.

% Copyright 2022 The MathWorks, Inc.

% Apply 3D Gaussian smoothing
postIn = squeeze(extractdata(postIn));
sigma = 0.5;
for n = 1:size(postIn,4)
    postIn(:,:,:,n) = imgaussfilt3(postIn(:,:,:,n),sigma,...
        FilterSize=3,Padding=0);   
end

if ~isempty(postFlipIn)
    postFlipIn = squeeze(extractdata(postFlipIn));
    % Apply 3D Gaussian smoothing on fliped data
    for n = 1:size(postFlipIn,4)
        postFlipIn(:,:,:,n) = imgaussfilt3(postFlipIn(:,:,:,n),...
            sigma,FilterSize=3,Padding=0);
    end
    % Get posteriors and segmentation
    postFlipIn = flip(postFlipIn,1);
    postFlipIn = flip(postFlipIn,2);
    postFlipIn = fliplr(postFlipIn);
    lrIndices = [2 3 4 5 6 7 8 9 10 11 15 16 17 18; 19 20 21 22 ...
        23 24 25 26 27 28 29 30 31 32];
    rlIndices = flip(lrIndices);
    postFlipIn(:,:,:,reshape(lrIndices',1,[])) = postFlipIn(:,:,:,reshape(rlIndices',1,[]));
    postIn = 0.5 * (postIn + postFlipIn);
end

% Use largest connected component
th = 0.25;
temp = postIn(:,:,:,2:end);
postInMask = sum(temp,4)>th;
largeComp = findLargestComponent(postInMask);
S = repmat(largeComp,1,1,1,size(temp,4));
temp(~S)=0;
postIn(:,:,:,2:end) = temp;

% Make posteriors to zero outside the largest connected component 
% of each topological class
postInMask = postIn> th;
topology_classes = [0  4  4  5  5  6  6  7  8  9 10  1  2  3  5 11 ...
    12 13 14 14 15 15 16 16 17 18 19 20 15 21 22 23];
topology_classesUnique = unique(topology_classes);
for topology_class = topology_classesUnique(2):topology_classesUnique(end)
    [~,tmp_topology_indices] = find(topology_classes==topology_class);
    tmp_mask = postInMask(:,:,:,tmp_topology_indices);
    tmp_mask =  findLargestComponent(tmp_mask);
    postIn(:,:,:,tmp_topology_indices) = ...
        postIn(:,:,:,tmp_topology_indices).*tmp_mask;
end

% Renormalize posteriors and get hard segmentation
postIn = postIn./sum(postIn,4);
[~,segPatch] = max(postIn,[],4);

% Make segmentation maps to original image size
seg = ones(imSize);
seg(cropIdx(1)+1:cropIdx(4),cropIdx(2)+1:cropIdx(5), ...
    cropIdx(3)+1:cropIdx(6)) = segPatch;
seg = labelIDs(seg); 

% Align prediction back to the input orientation
aff = eye(4);
affRef = metaData.Transform.T';
seg = alignBrainCANDIVolume(seg,aff,affRef);

% Convert segmentation result to categorical type
seg = categorical(seg,labelIDs,classNames);
end


function largeCompFinal =  findLargestComponent(postInMask)
largeCompFinal = zeros(size(postInMask));
for ii = 1:size(postInMask,4)
    tmp =  postInMask(:,:,:,ii);
    % get_largest_connected_component
    CC = bwconncomp(postInMask(:,:,:,ii),6);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    largeComp = false(size(tmp));
    largeComp(CC.PixelIdxList{idx}) = true;
    largeCompFinal(:,:,:,ii) = largeComp;
end
end
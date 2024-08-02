function [predictedSegMaps] = segment_space(imFile, dx_scan)

addpath('Brain_Segmentation')

str_sep = strfind(imFile, '.nii');
mat_file = strcat(imFile(1:str_sep(end)-1), "_seg.mat");

if ~exist(mat_file, "file")
    % Download label data
    zipFile = matlab.internal.examples.downloadSupportFile("image","data/brainSegData.zip");
    filepath = fileparts(zipFile);
    unzip(zipFile,filepath);
    
    dataDir = fullfile(filepath,"brainSegData");
    
    if ~exist(fullfile(dataDir, "trainedBrainSynthSegNetwork.h5"), "file")
        % Load pretrained network
        trainedBrainCANDINetwork_url = "https://www.mathworks.com/supportfiles/"+ ...
        "image/data/trainedBrainSynthSegNetwork.h5";
        downloadTrainedNetwork(trainedBrainCANDINetwork_url,dataDir);
    end
    
    % Load t1w data and labels
    metaData = niftiinfo(imFile);
    vol = niftiread(metaData);
    
    [classNames,labelIDs] = getBrainCANDISegmentationLabels;
    
    % Preprocessing
    des_dx = 1e-3; % assume 1 mm voxel resolution
    resample = dx_scan ~= des_dx;
    cropSize = 192;
    [volProc,cropIdx,imSize] = preProcessBrainCANDIData(vol,metaData,cropSize,resample);
    inputSize = size(volProc);
    
    volDL = dlarray(volProc,"SSSCB");
    
    % Define network architecture
    modelFile = fullfile(dataDir,"trainedBrainSynthSegNetwork.h5");
    lgraph = importKerasLayers(modelFile,ImportWeights=true,ImageInputSize=inputSize);
    
    % placeholderLayers = findPlaceholderLayers(lgraph);
    
    sf = softmaxLayer;
    lgraph = replaceLayer(lgraph,"unet_prediction",sf);
    net = dlnetwork(lgraph);
    layerGraph(net);
    
    % Predict
    predictIm = predict(net,volDL);

    flipVal = false;
    if flipVal
        flippedData = fliplr(volProc);  
        flippedData = flip(flippedData,2);
        flippedData = flip(flippedData,1);
        flippedData = dlarray(flippedData,"SSSCB");
        flipPredictIm = predict(net,flippedData);
    else
        flipPredictIm = [];  
    end

    predictedSegMaps = postProcessBrainCANDIData(predictIm,flipPredictIm,imSize,...
        cropIdx,metaData,classNames,labelIDs);
    
    % Save as mat file
    save(mat_file, "predictedSegMaps");
else
    % Load mat file
    predictedSegMaps = load(mat_file).predictedSegMaps;

end


clear all

%% -- Analysis parameters

componentConnectivity = 18;

binning = 3;

% Assign channels for different purposes
segChannel = 1; % Segmentation of nuclei


% Intensity threshold for cytoplasm segmentation
cyto_thresh = 300;

% Use parallel processing?
parallel_switch = true; % allowd values: true, false

ps = parallel.Settings;
ps.Pool.AutoCreate = parallel_switch;


%% --- pick directory and find .czi files

[sourceFile,sourceDir] = uigetfile('*.*');

thisPath = fullfile(sourceDir,sourceFile);

% --- Make a reader instance
reader = bfGetReader(thisPath);

% --- extract stack and microscope info from meta data
omeMeta = reader.getMetadataStore();
reader.close()

numChannels = omeMeta.getChannelCount(0);
numImages = omeMeta.getImageCount();
numFrames = omeMeta.getPixelsSizeT(0).getValue();

% --- get time stamps of the individual frames, in seconds
tt_vector = zeros(1,numFrames);
for kk = 1:numFrames
    
    tt_vector(kk) = omeMeta.getPlaneDeltaT(0,kk-1).value(ome.units.UNITS.SECOND);
    
end

tt_start_recording = tt_vector(1);

tt_vector = tt_vector-tt_vector(1);

% --- get the voxel edge sizes

voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0);
voxelSizeX = voxelSizeX.value(ome.units.UNITS.MICROM);
rawVoxelSizeX = voxelSizeX.doubleValue();
voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0);
voxelSizeY = voxelSizeY.value(ome.units.UNITS.MICROM);
rawVoxelSizeY = voxelSizeY.doubleValue();
voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0);
voxelSizeZ = voxelSizeZ.value(ome.units.UNITS.MICROM);
rawVoxelSizeZ = voxelSizeZ.doubleValue();

voxelSizeY = rawVoxelSizeY.*binning;
voxelSizeX = rawVoxelSizeY.*binning;
voxelSizeZ = rawVoxelSizeZ.*binning;

% --- get the spatial stack dimensions
rawStackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
rawStackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
rawStackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices

stackSizeY = floor(rawStackSizeY./binning);
stackSizeX = floor(rawStackSizeX./binning);
stackSizeZ = floor(rawStackSizeZ./binning);

fprintf('Raw stack dimensions: %d by %d by %d voxels.\n',...
    rawStackSizeY,rawStackSizeX,rawStackSizeZ)



%% --- Create boolean masks from cytoplasm


if parallel_switch

    fprintf('Starting to process frames in parallel.\n')
    fprintf('Progress dealing with individual frames not displayed.\n')
    fprintf('Showing only percentage of frames fully processed.\n')
    
    parfor_progress(numFrames);
    
end

errorFlagVec = false(1,numFrames);

segmentationMasks = cell(1,numFrames);

parfor ff = 1:numFrames

    if ~parallel_switch
        
        fprintf('Accessing frame %d of %d.\n',...
            ff,numFrames)
        
    end
    
    rawStack = cell(1,numChannels);
    
    if ~parallel_switch
        fprintf('Reading in stack.\n')
        
        parfor_progress(numChannels.*rawStackSizeZ);
    end
    
    % --- Make a reader instance
    reader = bfGetReader(thisPath);

    for cc = 1:numChannels
        
        rawStack{cc} = ...
            zeros(rawStackSizeY,rawStackSizeX,rawStackSizeZ,'uint16');
        
        for ll = 1:rawStackSizeZ
            
            % Direct assigment
            planeInd = reader.getIndex(ll-1,cc-1,ff-1)+1;
            rawStack{cc}(:,:,ll) = bfGetPlane(reader,planeInd);
            
            if ~parallel_switch
                parfor_progress;
            end
            
        end
        
    end

    reader.close();
    
    if ~parallel_switch
        parfor_progress(0);
    end
    
    
        
    
    
    % --- Make binned stack
    if binning > 1
        
        if ~parallel_switch
            fprintf('Binning stack, %d binning size.\n',binning)
        end
        
        binnedStack = cell(1,numChannels);
        
        [m,n,o]=size(rawStack{1});
        
        m = binning.*floor(m/binning);
        n = binning.*floor(n/binning);
        o = binning.*floor(o/binning);
        
        if ~parallel_switch
            parfor_progress(numChannels.*o);
        end
        
        for cc = 1:numChannels
            
            % dimension 1 binning
            
            xyBinStack = zeros(m/binning,n/binning,o/binning);
            
            for ll = 1:o
                
                container = squeeze(rawStack{cc}(1:m,1:n,ll));
                
                container = sum(reshape(container,binning,[]),1);
                container = reshape(container,m/binning,[]);
                container = permute(container,[2,1]);
                
                container = sum(reshape(container,binning,[]),1);
                container = reshape(container,n/binning,[]);
                
                xyBinStack(:,:,ll) = permute(container,[2,1]);
                
                if ~parallel_switch
                    parfor_progress;
                end
                
            end
            
            xyBinStack = permute(xyBinStack,[3,2,1]);
            xyBinStack = sum(reshape(xyBinStack,binning,[],m/binning),1);
            xyBinStack = reshape(xyBinStack,o/binning,[],m/binning);
            xyBinStack = permute(xyBinStack,[3,2,1]);
            
            binnedStack{cc} = xyBinStack./binning.^3;
            
        end
        
        if ~parallel_switch
            parfor_progress(0);
        end
        
    else
        
        binnedStack = cell(1,numChannels);

        for cc = 1:numChannels
        
            binnedStack{cc} = cast(rawStack{cc},'double');
        
        end
        
    end
    
    rawStack = 0;
        
    % --- cytoplasmic volume calculation
    
    % Take volume above Otsu threshold as cytoplasmic volume
    segmentationMasks{ff} = binnedStack{segChannel}>=cyto_thresh;
        
    if parallel_switch
        parfor_progress;
    end

end

parfor_progress(0);

reader.close()


%% --- summing up boolean masks over several frames

sumFrames = 30;

sumStartInds = [1:sumFrames:numFrames];
sumEndInds = sumStartInds+sumFrames-1;

sumEndInds(end) = min(sumEndInds(end),numFrames);

numSums = numel(sumStartInds);

summedMasks = cell(1,numSums);
cytoVol = zeros(1,numSums);
summed_tt = zeros(1,numSums);

maxXProj = cell(1,numSums);
maxYProj = cell(1,numSums);
maxZProj = cell(1,numSums);

for kk = 1:numSums
    
    summedMasks{kk} = segmentationMasks{sumStartInds(kk)};
    
    for nn = (1+sumStartInds(kk)):sumEndInds(kk)
       
        summedMasks{kk} = summedMasks{kk} | segmentationMasks{nn};
        
    end
    
    cytoVol(kk) = voxelSizeY.*voxelSizeX.*voxelSizeZ...
        .* sum(summedMasks{kk}(:));
    
    summed_tt(kk) = mean(tt_vector(sumStartInds(kk):sumEndInds(kk)));
    
    maxYProj{kk} = max(summedMasks{kk},[],1);
    maxXProj{kk} = max(summedMasks{kk},[],2);
    maxZProj{kk} = max(summedMasks{kk},[],3);
    
end






%% --- save analysis results

clear('omeMeta')
clear('reader')

save([thisPath,'_VolOutput.mat'])
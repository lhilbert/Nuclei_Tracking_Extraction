clear all

%% -- Analysis parameters

componentConnectivity = 18;

binning = 2;

% Do intensity evaluation only in xy, in the maximum contrast slice
onlyXY = true;

% Assign channels for different purposes
segChannel = 3; % Segmentation of nuclei

erodeSteps = 1;
dilateSpacer = 2;
dilateMeasure = 3;

% Intensity threshold for cytoplasm segmentation
cyto_thresh = 30;

% --- iterative thresholding parameters
adaptiveRangeFlag = false; % Switches adaptive threshold range on
minThresh = 55; % Lowest intensity threshold
maxThresh = 250; % Highest intensity threshold
threshSteps = 60; % How many threshold steps

distCutoff = 5; % How far can a centroid jump between two threshold values
%values before being considered a different object
volRatioCutoff = 0.5; % Change in volume to not trace object across
%thresholds, lower is more lenient, 1.0 the maximum stringency
segMinVol = 300; % minimum volume for nucleus to be recognized,
% unit: cubic micrometers
segMaxVol = 2.*1e5; % maximum volume for nucleus to be recognized,
% unit: cubic micrometers
segMinTrackLength = 3;

plot_nearest_neighbor = true;

% Use parallel processing?
parallel_switch = true; % allowd values: true, false

ps = parallel.Settings;
ps.Pool.AutoCreate = parallel_switch;

% Analyze Tokyo Tech data?
TokyoTechData = true;


%% --- pick directory and find .czi files

if ~TokyoTechData
    % in case of Zeiss .czi stacks
    
    [sourceFile,sourceDir] = uigetfile('*.*');
    
    thisPath = fullfile(sourceDir,sourceFile);
    
    % --- Make a reader instance
    reader = bfGetReader(thisPath);
    
    % --- extract stack and microscope info from meta data
    omeMeta = reader.getMetadataStore();
    reader.close()
    
    numChannels = omeMeta.getChannelCount(0);
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
    
    
    % --- get the spatial stack dimensions
    rawStackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
    rawStackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
    rawStackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices

else
    % in case of Tokyo Tech data set
    
    sourceDir = uigetdir;
    % --- Make a reader instance
    stackStruct =  read_TT_TiffStacks(sourceDir);
        
    numChannels = stackStruct.numChannels;
    numFrames = stackStruct.numTimes;
    
    % --- get time stamps of the individual frames, in seconds
    timeStep = 90; % time step between frames, in seconds
    % --- voxel edge sizes in micrometers
    rawVoxelSizeX = 0.828;
    rawVoxelSizeY = 0.828;
    rawVoxelSizeZ = 4;
    
    tt_vector = zeros(1,numFrames);
    for kk = 1:numFrames
        
        tt_vector(kk) = kk.*timeStep; % Time points in seconds
        
    end
    
    tt_start_recording = tt_vector(1);
    tt_vector = tt_vector-tt_vector(1);
    
    % --- get the spatial stack dimensions
    rawStackSizeX = stackStruct.rawStackSizeX; % image width, pixels
    rawStackSizeY = stackStruct.rawStackSizeY; % image height, pixels
    rawStackSizeZ = stackStruct.rawStackSizeZ; % number of Z slices
    
    [mainPath,DirName,~] = fileparts(sourceDir);
    thisPath = fullfile(mainPath,DirName);
    
end

voxelSizeY = rawVoxelSizeY.*binning;
voxelSizeX = rawVoxelSizeY.*binning;
voxelSizeZ = rawVoxelSizeZ.*binning;

stackSizeY = floor(rawStackSizeY./binning);
stackSizeX = floor(rawStackSizeX./binning);
stackSizeZ = floor(rawStackSizeZ./binning);

fprintf('Raw stack dimensions: %d by %d by %d voxels.\n',...
    rawStackSizeY,rawStackSizeX,rawStackSizeZ)



%% --- Segmentation and analysis within single frames

numNuc = zeros(1,numFrames);
nuc_vol_cell = cell(1,numFrames);
cyto_vol_vec = zeros(1,numFrames);
NN_median_vec = zeros(1,numFrames);
NN_distances_cell = cell(1,numFrames);

centroid_cell = cell(1,numFrames);

nuc_vxl_indices = cell(1,numFrames);

nucInt_cell = cell(1,numFrames);
cytoInt_cell = cell(1,numFrames);

ymaxProj_cell = cell(1,numFrames);
xmaxProj_cell = cell(1,numFrames);
zmaxProj_cell = cell(1,numFrames);

ysegProj_cell = cell(1,numFrames);
xsegProj_cell = cell(1,numFrames);
zsegProj_cell = cell(1,numFrames);

ycytoProj_cell = cell(1,numFrames);
xcytoProj_cell = cell(1,numFrames);
zcytoProj_cell = cell(1,numFrames);

thresh_vec = zeros(1,numFrames);

if parallel_switch
    
    fprintf('Starting to process frames in parallel.\n')
    fprintf('Progress dealing with individual frames not displayed.\n')
    fprintf('Showing only percentage of frames fully processed.\n')
    
    parfor_progress(numFrames);
    
end

errorFlagVec = false(1,numFrames);

for ff = 1:numFrames
    
    if ~parallel_switch
        
        fprintf('Accessing frame %d of %d.\n',...
            ff,numFrames)
        
    end
    
    try
        
        rawStack = cell(1,numChannels);
        
        if ~parallel_switch
            fprintf('Reading in stack.\n')
            
            parfor_progress(numChannels.*rawStackSizeZ);
        end
        
        if ~TokyoTechData
        
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
        
        elseif TokyoTechData
            
            for cc = 1:numChannels
                
                rawStack{cc} = ...
                    zeros(rawStackSizeY,rawStackSizeX,rawStackSizeZ,'uint16');
                
                for ll = 1:rawStackSizeZ
                    
                    % Direct assigment
                    rawStack{cc}(:,:,ll) = ...
                        read_TT_TiffStacks(stackStruct,[cc,ff,ll]);
                    
                    if ~parallel_switch
                        parfor_progress;
                    end
                    
                end
                
            end
            
        end
        
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
        
        voxelVol = voxelSizeY.*voxelSizeX.*voxelSizeZ;
        
        ymaxProj_cell{ff} = squeeze(max(binnedStack{segChannel},[],1));
        xmaxProj_cell{ff} = squeeze(max(binnedStack{segChannel},[],2));
        zmaxProj_cell{ff} = squeeze(max(binnedStack{segChannel},[],3));
        
        
        
        volCell = cell(1,threshSteps);
        centrCell = cell(1,threshSteps);
        vxlIdxCell = cell(1,threshSteps);
        coordsVxlIdxCell = cell(1,threshSteps);
        imagesCell = cell(1,threshSteps);
        
        numRegVec = zeros(1,threshSteps);
        
        featureVecCell = cell(1,threshSteps);
        
        if ~parallel_switch
            fprintf('Nuclei segmentation with iterative thresholds.\n')
        end
        
        % --- get threshold values
        
        jointInts = [ymaxProj_cell{ff}(:);...
            xmaxProj_cell{ff}(:);...
            zmaxProj_cell{ff}(:)];
        
        if adaptiveRangeFlag
            
            threshVals = ...
                linspace(cyto_thresh,max(jointInts),threshSteps+1);
            threshVals = threshVals(1:end-1);
            
        else
            
            threshVals = ...
                linspace(minThresh,maxThresh,threshSteps);
            
        end
        
        if ~parallel_switch
            parfor_progress(threshSteps);
        end
        
        % --- cytoplasmic volume calculation
        
        % Take volume above Otsu threshold as cytoplasmic volume
        cytoStack = binnedStack{segChannel}>=cyto_thresh;
        cyto_vol_vec(ff) = sum(cytoStack(:)).*voxelVol;
        ycytoProj_cell{ff} = squeeze(max(cytoStack,[],1));
        xcytoProj_cell{ff} = squeeze(max(cytoStack,[],2));
        zcytoProj_cell{ff} = squeeze(max(cytoStack,[],3));
        
        for tt = 1:threshSteps
            
            currThr = threshVals(tt);
            
            regions = bwconncomp(binnedStack{segChannel}>=currThr,...
                componentConnectivity);
            regionProps = regionprops(regions,'Area');
            regionVols = [regionProps(:).Area].*voxelVol;
            
            % Restrict to objects greater than minimum volume
            validVolInds = regionVols>segMinVol & regionVols<segMaxVol;
            nucleiRegions = struct;
            nucleiRegions.Connectivity = regions.Connectivity;
            nucleiRegions.ImageSize = regions.ImageSize;
            nucleiRegions.NumObjects = sum(validVolInds);
            nucleiRegions.PixelIdxList = regions.PixelIdxList(validVolInds);
            
            regionProps = ...
                regionprops(nucleiRegions,'Area','Centroid',...
                'PixelIdxList');
            volumes = [regionProps(:).Area].*voxelVol;
            centroids = {regionProps(:).Centroid};
            vxlIdx = {regionProps(:).PixelIdxList};
            
            yCoords = cellfun(@(elmt)elmt(1).*voxelSizeY,centroids);
            xCoords = cellfun(@(elmt)elmt(2).*voxelSizeX,centroids);
            zCoords = cellfun(@(elmt)elmt(3).*voxelSizeZ,centroids);
            
            featureVecCell{tt} = [volumes;xCoords;yCoords;zCoords].';
            
            numRegVec(tt) = numel(volumes);
            volCell{tt} = volumes;
            centrCell{tt} = centroids;
            vxlIdxCell{tt} = vxlIdx;
            
            if ~parallel_switch
                parfor_progress;
            end
            
        end
        
        if ~parallel_switch
            parfor_progress(0);
        end
        
        % --- Linking of objects
        
        if ~parallel_switch
            fprintf('Selection of individual nuclei thresholds.\n')
            
            parfor_progress(2.*threshSteps);
        end
        rvsLinkInds = cell(1,threshSteps-1);
        physicalDist = cell(1,threshSteps-1);
        volRatio = cell(1,threshSteps-1);
        
        
        for tt = 1:(threshSteps-1)
            
            if ~isempty(featureVecCell{tt+1}) ...
                    && ~isempty(featureVecCell{tt})
                
                distMatr = ...
                    pdist2(featureVecCell{tt+1}(:,2:end),...
                    featureVecCell{tt}(:,2:end),'euclidean');
                
                [rvsLinkDist,rvsLinkInds{tt}] = min(distMatr,[],2);
                
                physicalDist{tt} = sqrt(sum((...
                    featureVecCell{tt}(rvsLinkInds{tt},2:end) ...
                    - featureVecCell{tt+1}(:,2:end)).^2,2));
                
                volRatio{tt} = ...
                    featureVecCell{tt+1}(:,1)...
                    ./featureVecCell{tt}(rvsLinkInds{tt},1);
                volRatio{tt}(volRatio{tt}>1) = 1./volRatio{tt}(volRatio{tt}>1);
                
            else
                rvsLinkInds{tt} = [];
                physicalDist{tt} = [];
                volRatio{tt} = [];
            end
            
            if ~parallel_switch
                parfor_progress;
            end
            
        end
        
        objects = {};
        if numRegVec(end)>0
            for rr = 1:numRegVec(end)
                
                thisRegion = struct;
                
                thisRegion.origObjIndTrack = rr;
                thisRegion.volTrack = volCell{end}(rr);
                thisRegion.threshTrack = threshSteps;
                thisRegion.centrTrack = centrCell{end}(rr);
                thisRegion.truncated = false;
                
                objects = [objects,{thisRegion}];
                
            end
        end
        
        %%% --- beginning bug hunt iterative thresholding gap
        
        % Here, we start make objects with traces through the different
        % segmentation intensities. In each step, we check for every region
        % that has been segmented at this intensity, if a region from the
        % previous intensity maps to it. If so, we assign it to the object that
        % was with this region at the previous intensity. If not, we start a
        % new object from this region. This way we will have all regions at all
        % threshold values contained within objects.
        
        for tt = (threshSteps-1):-1:1
            
            nonTruncInds = ...
                find(cellfun(@(elmt)~elmt.truncated,objects));
            
            usedInds = [];
            
            if ~isempty(rvsLinkInds{tt})
                
                % make all objects truncated, and only bring those back that
                % get connected
                
                for oo = 1:numel(objects)
                    objects{oo}.truncated = true;
                end
                
                useFlags = ...
                    (volRatio{tt}>=volRatioCutoff ...
                    & physicalDist{tt}<=distCutoff);
                
                %                 disp('+++')
                %                 disp(numRegVec(tt+1))
                %                 disp(numel(nonTruncInds))
                
                for rr = 1:numRegVec(tt+1)
                    
                    objInd = nonTruncInds(rr);
                    origObjInd = objects{objInd}.origObjIndTrack(end);
                    rvsInd = rvsLinkInds{tt}(origObjInd);
                    
                    if useFlags(origObjInd)
                        
                        thisRegion = struct;
                        
                        objects{objInd}.origObjIndTrack = ...
                            [objects{objInd}.origObjIndTrack,rvsInd];
                        objects{objInd}.volTrack = ...
                            [objects{objInd}.volTrack,volCell{tt}(rvsInd)];
                        objects{objInd}.threshTrack = ...
                            [objects{objInd}.threshTrack,tt];
                        objects{objInd}.centrTrack = ...
                            [objects{objInd}.centrTrack,centrCell{tt}(rvsInd)];
                        objects{objInd}.truncated = false;
                        
                        usedInds = [usedInds,rvsInd];
                        
                    end
                    
                end
                
            end
            
            notUsedInds = setdiff(1:numRegVec(tt),usedInds);
            
            
            for rr = 1:numel(notUsedInds)
                
                regInd = notUsedInds(rr);
                
                thisRegion = struct;
                
                thisRegion.origObjIndTrack = regInd;
                thisRegion.volTrack = volCell{tt}(regInd);
                thisRegion.threshTrack = tt;
                thisRegion.centrTrack = centrCell{tt}(regInd);
                thisRegion.truncated = false;
                
                objects = [objects,{thisRegion}];
                
            end
            
            if ~parallel_switch
                parfor_progress;
            end
            
        end
        
        
        % -- Find individual objects' threshold by minimum volume change
        
        keepFlags = ...
            cellfun(@(elmt)numel(elmt.threshTrack)>segMinTrackLength,objects);
        keepObjects = objects(keepFlags);
        
        for oo = 1:numel(keepObjects)
            
            % normalized volume difference between different thresholds
            volDiff = diff(keepObjects{oo}.volTrack) ...
                ./keepObjects{oo}.volTrack(1:end-1);
            
            minDiff = min(volDiff);
            minDiffInd = find(volDiff==minDiff,1,'last');
            
            keepObjects{oo}.threshInd = keepObjects{oo}.threshTrack(minDiffInd);
            keepObjects{oo}.centr = keepObjects{oo}.centrTrack{minDiffInd};
            keepObjects{oo}.vol = keepObjects{oo}.volTrack(minDiffInd);
            keepObjects{oo}.origObjInd = ...
                keepObjects{oo}.origObjIndTrack(minDiffInd);
            keepObjects{oo}.vxlIdx = ...
                vxlIdxCell{keepObjects{oo}.threshInd}{keepObjects{oo}.origObjInd};
            
        end
        
        if ~parallel_switch
            parfor_progress;
        end
        
        % -- rejection of overlapping objects
        
        numObjs = numel(keepObjects);
        keepFlags = true(1,numObjs);
        
        yCoords = cellfun(@(elmt)elmt.centr(1).*voxelSizeY,keepObjects);
        xCoords = cellfun(@(elmt)elmt.centr(2).*voxelSizeX,keepObjects);
        zCoords = cellfun(@(elmt)elmt.centr(3).*voxelSizeZ,keepObjects);
        
        coordVector = [xCoords.',yCoords.',zCoords.'];
        objectDist = pdist2(coordVector,coordVector);
        
        for rr = (numRegVec(end)+1):numObjs
            
            %             numObjs - rr
            
            checkInds = 1:(rr-1);
            [dists,sortInds] = sort(objectDist(rr,checkInds),'ascend');
            
            sortInds = sortInds(dists<50); % limit to a sphere of 50 microns radius
            
            checkInds = checkInds(sortInds);
            
            for qq = checkInds;
                
                if ~isempty(intersect(keepObjects{rr}.vxlIdx,keepObjects{qq}.vxlIdx))
                    
                    keepFlags(rr) = false;
                    break;
                    
                end
                
            end
            
        end
        
        
        %%%%% --- end of bug hunt "iterative thresholding gap"
        
        
        validObjs = keepObjects(keepFlags);
        
        % -- Make binary image of valid objects
        % Add all objects to binary 3D image
        segArray = false(stackSizeY,stackSizeX,stackSizeZ);
        numValObj = numel(validObjs);
        for oo = 1:numValObj
            segArray(validObjs{oo}.vxlIdx) = true;
        end
        
        % Store maximum z projection and central x section
        ysegProj_cell{ff} = max(segArray,[],1);
        xsegProj_cell{ff} = max(segArray,[],2);
        zsegProj_cell{ff} = max(segArray,[],3);
        
        if ~parallel_switch
            parfor_progress;
            parfor_progress(0);
        end
        
        % -- Quantification based on extracted regions
        
        if ~parallel_switch
            fprintf('Quantification of nuclear and cytoplasmic intensities.\n')
        end
        
        % Construct region structure for further intensity feature
        % extraction
        nucleiRegions = struct;
        nucleiRegions.Connectivity = componentConnectivity;
        nucleiRegions.ImageSize = [stackSizeY,stackSizeX,stackSizeZ];
        nucleiRegions.NumObjects = numel(validObjs);
        nucleiRegions.PixelIdxList = ...
            cellfun(@(elmt)elmt.vxlIdx,validObjs,'UniformOutput',false);
        
        nucleiPropsBinary = regionprops(nucleiRegions,...
            'Area','Image','BoundingBox');
        
        nucleiVol = [nucleiPropsBinary(:).Area].*voxelVol;
        nucleiImage = {nucleiPropsBinary(:).Image};
        nucleiBoundingBox = {nucleiPropsBinary(:).BoundingBox};
        nucleiVxlIdx = nucleiRegions.PixelIdxList;
        
        nucleiPropsSegmentation = ...
            regionprops(nucleiRegions,binnedStack{segChannel},...
            'WeightedCentroid');
        
        nucleiCentroid = ...
            {nucleiPropsSegmentation(:).WeightedCentroid};
        
        numNuc(ff) = numel(nucleiVol);
        nuc_vol_cell{ff} = nucleiVol;
        centroid_cell{ff} = nucleiCentroid;
        nuc_vxl_indices{ff} = nucleiVxlIdx;
        for ll = 1:numNuc(ff)
            centroid_cell{ff}{ll} = ...
                centroid_cell{ff}{ll}.*[voxelSizeY,voxelSizeX,voxelSizeZ];
        end
        
        nucInt = cell(1,numChannels);
        cytoInt = cell(1,numChannels);
        
        for cc = 1:numChannels
            
            nucInt{cc} = zeros(1,numNuc(ff));
            cytoInt{cc} = zeros(1,numNuc(ff));
            
        end
        
        
        % Make masks to determine cytoplasmic intensity
        
        totalDil = dilateSpacer+dilateMeasure;
        
        if ~parallel_switch
            parfor_progress(numNuc(ff));
        end
        
        if onlyXY
            
            % --- 2D treatment in maximum contrast slice
            
            for nn = 1:numNuc(ff)
                
                % Determine small bounding box (without dilate)
                smallMinY = nucleiBoundingBox{nn}(2)+0.5;
                smallMinX = nucleiBoundingBox{nn}(1)+0.5;
                MinZ = nucleiBoundingBox{nn}(3)+0.5;
                smallMaxY = nucleiBoundingBox{nn}(2)+nucleiBoundingBox{nn}(5)-0.5;
                smallMaxX = nucleiBoundingBox{nn}(1)+nucleiBoundingBox{nn}(4)-0.5;
                MaxZ = nucleiBoundingBox{nn}(3)+nucleiBoundingBox{nn}(6)-0.5;
                
                % Determine extended bounding box (after dilate)
                fullExtMinY = smallMinY-totalDil;
                fullExtMinX = smallMinX-totalDil;
                fullExtMaxY = smallMaxY+totalDil;
                fullExtMaxX = smallMaxX+totalDil;
                
                % Limit extended bounding box to within image limits
                extMinY = max(1,fullExtMinY);
                yLoDiff = extMinY - fullExtMinY;
                extMinX = max(1,fullExtMinX);
                xLoDiff = extMinX - fullExtMinX;
                extMaxY = min(stackSizeY,fullExtMaxY);
                yHiDiff = fullExtMaxY - extMaxY;
                extMaxX = min(stackSizeX,fullExtMaxX);
                xHiDiff = fullExtMaxX - extMaxX;
                
                % Extended bounding box size
                extSizeY = extMaxY - extMinY + 1;
                extSizeX = extMaxX - extMinX + 1;
                SizeZ = MaxZ - MinZ + 1;
                
                cutoutImage = binnedStack{cc}(extMinY:extMaxY,...
                    extMinX:extMaxX,...
                    MinZ:MaxZ);
                
                % Find maximum contrast slice
                
                contrastVec = zeros(1,rawStackSizeZ);
                
                for ll = 1:SizeZ
                    
                    section = cutoutImage(:,:,ll);
                    
                    diffMatrix = ...
                        (section(1:end-1,1:end-1)-section(2:end,2:end)).^2;
                    
                    contrastVec(ll) = sqrt(mean(diffMatrix(:)));
                    
                end
                
                [~,maxContrastInd] = max(contrastVec);
                
                maxContrastImage = ...
                    squeeze(cutoutImage(:,:,maxContrastInd));
                
                % Inclusion mask
                inclMask = zeros(extSizeY,extSizeX);
                inclMask((1+totalDil-yLoDiff):(end-totalDil+yHiDiff),...
                    (1+totalDil-xLoDiff):(end-totalDil+xHiDiff))...
                    = squeeze(nucleiImage{nn}(:,:,maxContrastInd));
                
                % Nucleus mask
                nucMask = inclMask>0;
                
                % Exclusion mask
                exclMask = squeeze(segArray(extMinY:extMaxY,...
                    extMinX:extMaxX,...
                    MinZ+maxContrastInd-1));
                
                dilateNeighborhood = zeros(3,3);
                dilateNeighborhood(2,2) = 1;
                dilateNeighborhood(1,2) = 1;
                dilateNeighborhood(3,2) = 1;
                dilateNeighborhood(2,1) = 1;
                dilateNeighborhood(2,3) = 1;
                
                for ee = 1:erodeSteps
                    nucMask = imerode(nucMask,dilateNeighborhood);
                end
                
                for dd = 1:dilateSpacer
                    inclMask = imdilate(inclMask,dilateNeighborhood);
                    exclMask = imdilate(exclMask,dilateNeighborhood);
                end
                
                for dd = 1:dilateMeasure
                    inclMask = imdilate(inclMask,dilateNeighborhood);
                end
                
                measureMask = inclMask & ~exclMask;
                
                for cc = 1:numChannels
                    
                    % -- Determine nuclear and cytoplasmic intensities
                    
                    quantImage = binnedStack{cc}(...
                        extMinY:extMaxY,extMinX:extMaxX,...
                        MinZ+maxContrastInd-1);
                    
                    cytoInt{cc}(nn) = mean(quantImage(measureMask));
                    nucInt{cc}(nn) = mean(quantImage(nucMask));
                    
                end
                
                if ~parallel_switch
                    parfor_progress;
                end
                
            end
            
        else
            
            % Full 3D treatment
            
            for nn = 1:numNuc(ff)
                
                % Determine small bounding box (without dilate)
                smallMinY = nucleiBoundingBox{nn}(2)+0.5;
                smallMinX = nucleiBoundingBox{nn}(1)+0.5;
                smallMinZ = nucleiBoundingBox{nn}(3)+0.5;
                smallMaxY = nucleiBoundingBox{nn}(2)+nucleiBoundingBox{nn}(5)-0.5;
                smallMaxX = nucleiBoundingBox{nn}(1)+nucleiBoundingBox{nn}(4)-0.5;
                smallMaxZ = nucleiBoundingBox{nn}(3)+nucleiBoundingBox{nn}(6)-0.5;
                
                % Determine extended bounding box (after dilate)
                fullExtMinY = smallMinY-totalDil;
                fullExtMinX = smallMinX-totalDil;
                fullExtMinZ = smallMinZ-totalDil;
                fullExtMaxY = smallMaxY+totalDil;
                fullExtMaxX = smallMaxX+totalDil;
                fullExtMaxZ = smallMaxZ+totalDil;
                
                % Limit extended bounding box to within image limits
                extMinY = max(1,fullExtMinY);
                yLoDiff = extMinY - fullExtMinY;
                extMinX = max(1,fullExtMinX);
                xLoDiff = extMinX - fullExtMinX;
                extMinZ = max(1,fullExtMinZ);
                zLoDiff = extMinZ - fullExtMinZ;
                extMaxY = min(stackSizeY,fullExtMaxY);
                yHiDiff = fullExtMaxY - extMaxY;
                extMaxX = min(stackSizeX,fullExtMaxX);
                xHiDiff = fullExtMaxX - extMaxX;
                extMaxZ = min(stackSizeZ,fullExtMaxZ);
                zHiDiff = fullExtMaxZ - extMaxZ;
                
                % Extended bounding box size
                extSizeY = extMaxY - extMinY + 1;
                extSizeX = extMaxX - extMinX + 1;
                extSizeZ = extMaxZ - extMinZ + 1;
                
                % Inclusion mask
                inclMask = zeros(extSizeY,extSizeX,extSizeZ);
                inclMask((1+totalDil-yLoDiff):(end-totalDil+yHiDiff),...
                    (1+totalDil-xLoDiff):(end-totalDil+xHiDiff),...
                    (1+totalDil-zLoDiff):(end-totalDil+zHiDiff))...
                    = nucleiImage{nn};
                
                % Nucleus mask
                nucMask = inclMask>0;
                
                % Exclusion mask
                exclMask = segArray(extMinY:extMaxY,...
                    extMinX:extMaxX,...
                    extMinZ:extMaxZ);
                
                dilateNeighborhood = zeros(3,3,3);
                dilateNeighborhood(2,2,2) = 1;
                dilateNeighborhood(1,2,2) = 1;
                dilateNeighborhood(3,2,2) = 1;
                dilateNeighborhood(2,1,2) = 1;
                dilateNeighborhood(2,3,2) = 1;
                dilateNeighborhood(2,2,1) = 1;
                dilateNeighborhood(2,2,3) = 1;
                
                for ee = 1:erodeSteps
                    nucMask = imerode(nucMask,dilateNeighborhood);
                end
                
                for dd = 1:dilateSpacer
                    inclMask = imdilate(inclMask,dilateNeighborhood);
                    exclMask = imdilate(exclMask,dilateNeighborhood);
                end
                
                for dd = 1:dilateMeasure
                    inclMask = imdilate(inclMask,dilateNeighborhood);
                end
                
                measureMask = inclMask & ~exclMask;
                
                for cc = 1:numChannels
                    
                    % -- Determine nuclear and cytoplasmic intensities
                    
                    cutoutImage = binnedStack{cc}(extMinY:extMaxY,...
                        extMinX:extMaxX,...
                        extMinZ:extMaxZ);
                    
                    cytoInt{cc}(nn) = mean(cutoutImage(measureMask));
                    nucInt{cc}(nn) = mean(cutoutImage(nucMask));
                    
                end
                
                if ~parallel_switch
                    parfor_progress;
                end
                
            end
            
        end
        
        
        if ~parallel_switch
            parfor_progress(0);
        end
        
        % --- Nearest neighbor distances
        
        if numNuc(ff) <= 2
            
            distances = [];
            median_val = NaN;
            
        else
            
            if ~parallel_switch
                fprintf('Calculating nearest neighbor distances...')
            end
            
            yCoords = cellfun(@(elmt)elmt(1),centroid_cell{ff});
            xCoords = cellfun(@(elmt)elmt(2),centroid_cell{ff});
            zCoords = cellfun(@(elmt)elmt(3),centroid_cell{ff});
            
            coord_matrix = [yCoords;xCoords;zCoords].';
            
            % calculate median pairwise distance to nearest neighbor
            distMatrix = squareform(pdist(coord_matrix));
            distMatrix(distMatrix == 0) = Inf;
            [distances,NNinds] = min(distMatrix);
            
            median_val = median(distances);
            
            if ~parallel_switch
                fprintf('done.\n')
            end
            
        end
        
        NN_distances_cell{ff} = distances;
        NN_median_vec(ff) = median_val;
        
        
        % Upon successful completion, set error flag to false, no error!
        errorFlagVec(ff) = false;
        
    catch message
        
        errorFlagVec(ff) = true;
        errorMessages{ff} = message;
        
        fprintf('Image segmentation error in step %d.\n',kk)
        
    end
    
    
    % Assign storage variables based on this frames's results
    numNuc(ff) = numNuc(ff);
    nuc_vol_cell{ff} = nucleiVol;
    
    nucInt_cell{ff} = nucInt;
    cytoInt_cell{ff} = cytoInt;
    
    if parallel_switch
        parfor_progress;
    end
    
end

parfor_progress(0);

%% --- save analysis results

clear('omeMeta')
clear('reader')

save([thisPath,'_AnalysisOutput.mat'])
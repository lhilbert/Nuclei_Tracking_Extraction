function varargout = read_TT_TiffStacks(varargin)

if nargin == 1
    
    listing = dir(fullfile(varargin{1},'*.tif'));
    
    tiffStackAccessStruct = struct;
    
    tiffStackAccessStruct.sourceDir = varargin{1};
    tiffStackAccessStruct.listing = listing;
    
    numImages = numel(listing);
    
    channelInds = uint16(zeros(1,numImages));
    sectionInds = uint16(zeros(1,numImages));
    timeInds = uint16(zeros(1,numImages));
    
    for kk = 1:numImages
        
        indString = listing(kk).name(end-15:end-4);
        
        channelInds(kk) = str2double(indString(2:4));
        sectionInds(kk) = str2double(indString(6:8));
        timeInds(kk) = str2double(indString(10:12));
        
    end
    
    tiffStackAccessStruct.channelList = unique(channelInds);
    tiffStackAccessStruct.timeList = unique(timeInds);
    tiffStackAccessStruct.sectionList = unique(sectionInds);

    tiffStackAccessStruct.channelInds = channelInds;
    tiffStackAccessStruct.timeInds = timeInds;
    tiffStackAccessStruct.sectionInds = sectionInds;
    
    tiffStackAccessStruct.numChannels = ...
        numel(tiffStackAccessStruct.channelList);
    tiffStackAccessStruct.numTimes = ...
        numel(tiffStackAccessStruct.timeList);
    tiffStackAccessStruct.numSections = ...
        numel(tiffStackAccessStruct.sectionList);

    
    % Run through listing again and route
    
    tiffPathCell = cell(...
        tiffStackAccessStruct.numChannels,...
        tiffStackAccessStruct.numTimes,...
        tiffStackAccessStruct.numSections);
    
    for kk = 1:numImages
        
        tiffPathCell{channelInds(kk),...
            timeInds(kk),sectionInds(kk)} = ...
            fullfile(varargin{1},listing(kk).name);
        
    end
    
    tiffStackAccessStruct.tiffPathCell = tiffPathCell;
    
    intensityImage = ...
        imread(tiffStackAccessStruct.tiffPathCell{1,1,1});
    
    [ySize,xSize] = size(intensityImage);
    
    tiffStackAccessStruct.rawStackSizeX = xSize;
    tiffStackAccessStruct.rawStackSizeY = ySize;
    tiffStackAccessStruct.rawStackSizeZ = ...
        tiffStackAccessStruct.numSections;
    
    varargout{1} = tiffStackAccessStruct;
    
elseif nargin > 1
    
    tiffStackAccessStruct = varargin{1};
    accessInds = varargin{2};
    
    intensityImage = ...
        imread(tiffStackAccessStruct.tiffPathCell{...
        accessInds(1),accessInds(2),accessInds(3)});
            
    varargout{1} = intensityImage;
    
end
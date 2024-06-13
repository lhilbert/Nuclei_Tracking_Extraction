clear all

%% --- analysis parameters

intChannel = 1; % Channel to get intensity ratios from

%% --- Load raw analysis results
    
[sourceFile,sourceDir] = uigetfile('*.*');

thisPath = fullfile(sourceDir,sourceFile);

load(thisPath);


%% -- Get nearest neighbor distance only for sufficiently bright nuclei

minIntRatio = 2;
minNucCount = 8;

selectedNucleiCentroids = cell(1,numFrames);

numValidNuc = zeros(1,numFrames);

validDistCell = cell(1,numFrames);
validDistMedian = zeros(1,numFrames);

for ff = 1:numFrames
    
   % Indices of nuclei with high enough intensity 
   
   if numNuc(ff)>0
       
       nucInts = nucInt_cell{ff}{intChannel};
       cytoInts = cytoInt_cell{ff}{intChannel};
       
       intRatios = nucInts./cytoInts;
       
       validNucInds = find(intRatios>=minIntRatio);

       numValidNuc(ff) = numel(validNucInds);

       selectedNucleiCentroids{ff} = ...
           centroid_cell{ff}(validNucInds);
       
   else
       
       numValidNuc(ff) = 0;
       selectedNucleiCentroids{ff} = {};
       
   end
   
   
   
   
   % --- Nearest neighbor distances
    
    if numValidNuc(ff) <= 2
    
        distances = [];
        median_val = NaN;
        var_val = NaN;
        
    else
        
        yCoords = cellfun(@(elmt)elmt(1),selectedNucleiCentroids{ff});
        xCoords = cellfun(@(elmt)elmt(2),selectedNucleiCentroids{ff});
        zCoords = cellfun(@(elmt)elmt(3),selectedNucleiCentroids{ff});

        coord_matrix = [yCoords;xCoords;zCoords].';
        
        % calculate median pairwise distance to second-nearest neighbor
        distMatrix = squareform(pdist(coord_matrix));
        distMatrix(distMatrix == 0) = Inf;
        distMatrix = sort(distMatrix);        
        distances = squeeze(distMatrix(2,:));
        
        median_val = median(distances);
        var_val = var(distances);
    
    end
        
    validDistCell{ff} = distances;
    validDistMedian(ff) = median_val;
    validDistVar(ff) = var_val;
    
end




% --- plot overview of nuclei counts and nearest neighbor distances

upper_lim = cellfun(@(distr) prctile(distr,70),validDistCell);
lower_lim = cellfun(@(distr) prctile(distr,30),validDistCell);

subplot(1,2,inputFile)

errorbar(tt_vector(numValidNuc>minNucCount)./60,...
    validDistMedian(numValidNuc>minNucCount),...
    lower_lim(numValidNuc>minNucCount)...
    -validDistMedian(numValidNuc>minNucCount),...
    upper_lim(numValidNuc>minNucCount)...
    -validDistMedian(numValidNuc>minNucCount),...
    'k-o')

xlabel('Recording time [min]')
ylabel('Neighbor distance [\mum]')
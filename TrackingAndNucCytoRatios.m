clear all

%% --- tracking parameters

maxMoveDist = 25; % in unit of micrometers

minIntRatio = 1.15; % minimum intensity increase in nucleus
% vs. surrounding cytoplasm

minVol = 80; % minimum object volume in cubic microns
maxVol = 1800; % maximum volume

intChannel = 2; % Channel to get intensity ratios from

inputFile = 0; % 0 - new file, 1,2,{3,4} correspond to known data sets

deltat = 90; % time interval in seconds

%% --- Load raw analysis results

if inputFile == 0
    
    [load_sourceFile,load_sourceDir] = uigetfile('*.*');
    
elseif inputFile == 1
    
    load_sourceDir = '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/';
    load_sourceFile = ...
        '10_G1_Subset_fullTimeCourse.czi_AnalysisOutput.mat';

elseif inputFile == 2
    
    load_sourceDir = '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/';
    load_sourceFile = ...
        '10_G2_Subset.czi_AnalysisOutput.mat';

elseif inputFile == 3
    
    load_sourceDir = ...
        '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/PCNASubsets_5Feb2016/';
    load_sourceFile = ...
        'Image_16_Subset.czi_AnalysisOutput.mat';

elseif inputFile == 4
    
    load_sourceDir = ...
        '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/PCNASubsets_5Feb2016/';
    load_sourceFile = ...
        'Image_18_Subset.czi_AnalysisOutput.mat';

end

thisPath = fullfile(load_sourceDir,load_sourceFile);

load(thisPath);

tt_vector = 0:deltat:deltat.*(numFrames-1);

%% -- Initial plotting and evaluation

figure(1)

clf

maxInt = 0;

jointRatios = [];

for kk = 1:numFrames
    
    validInds = nuc_vol_cell{kk}>=minVol & nuc_vol_cell{kk}<=maxVol;
    
    subplot(1,2,1)
    
    plot(cytoInt_cell{kk}{intChannel}(validInds),...
        nucInt_cell{kk}{intChannel}(validInds)-cytoInt_cell{kk}{intChannel}(validInds),'ko',...
        'MarkerEdgeColor',[1,0.6,0.6].*(kk-1)./(numFrames-1)) 
    %     hold on
        
    set(gca,'XLim',[0,10000],...
        'YLim',[0,10000])

    subplot(1,2,2)
    
    theseRatios = nucInt_cell{kk}{intChannel}(validInds)...
        ./cytoInt_cell{kk}{intChannel}(validInds);
    
    plot(ones(1,numel(theseRatios)).*kk.*1.5,theseRatios,'ko')
    
    hold on
    
    jointRatios = [jointRatios,theseRatios];
    
    maxInt = max([maxInt,max(nucInt_cell{kk}{intChannel}(validInds))]);

    waitforbuttonpress;
    
end

subplot(1,2,1)

plot([0,maxInt],[0,maxInt],'b--')

plot([0,maxInt],[0,maxInt.*minIntRatio],'b-.')

xlabel('I_{cyto}')
ylabel('I_{nuc}')

hold off

axis equal

set(gca,'XLim',[0,maxInt],...
    'YLim',[0,maxInt])

subplot(1,2,2)


xlabel('Time')


%% --- Nucleus tracking

trackInds = cell(1,numFrames);
trackInds{1} = 1:numNuc(1);

rvsInds = cell(1,numFrames-1);
fwdInds = cell(1,numFrames-1);

for kk = 2:numFrames
    
    if numNuc(kk)>0 && numNuc(kk-1)>0
        
        xCoords = cellfun(@(elmt)elmt(1),centroid_cell{kk-1}).';
        yCoords = cellfun(@(elmt)elmt(2),centroid_cell{kk-1}).';
        zCoords = cellfun(@(elmt)elmt(3),centroid_cell{kk-1}).';
        
        propMatrix_last = [xCoords,yCoords,zCoords];
        
        xCoords = cellfun(@(elmt)elmt(1),centroid_cell{kk}).';
        yCoords = cellfun(@(elmt)elmt(2),centroid_cell{kk}).';
        zCoords = cellfun(@(elmt)elmt(3),centroid_cell{kk}).';
        
        propMatrix_this = [xCoords,yCoords,zCoords];
        
        dist = pdist2(propMatrix_last,propMatrix_this,'seuclidean');
        
        [fwdDist,fwdInds{kk-1}] = min(dist,[],2);
        [rvsDist,rvsInds{kk-1}] = min(dist,[],1);
        
        trackInds{kk} = trackInds{kk-1}(rvsInds{kk-1});
        
    elseif numNuc(kk)>0
        
        trackInds{kk} = 1:numNuc(kk);
        
    end
    
end

% Get tracking connection length statistics, to reject jumps
moveDist = cell(1,numFrames-1);
moveDistPrctl = zeros(1,numFrames-1);
moveDistMax = zeros(1,numFrames-1);

for kk = 1:(numFrames-1)
    
    if numNuc(kk)>0 && numNuc(kk+1)>0
        
        moveDist{kk} = zeros(1,numNuc(kk+1));
        
        for nn = 1:numNuc(kk+1)
            
            currentCent = centroid_cell{kk+1}{nn};
            lastCent = centroid_cell{kk}{rvsInds{kk}(nn)};
            
            moveDist{kk}(nn) = sqrt(sum((currentCent-lastCent).^2));
            
        end
        
        moveDistPrctl(kk) = prctile(moveDist{kk},[97]);
        moveDistMax(kk) = max(moveDist{kk});
        
    end
    
end

for kk = 2:numFrames
    
    if numNuc(kk-1)>0 && numNuc(kk)>0
        
        for nn = 1:numNuc(kk)
            
            currentCent = centroid_cell{kk}{nn};
            lastCent = centroid_cell{kk-1}{rvsInds{kk-1}(nn)};
            
            % target intensity ratio check
            currentIntRatio = ...
                nucInt_cell{kk}{intChannel}(nn) ...
                ./cytoInt_cell{kk}{intChannel}(nn);
            
            
            if sqrt(sum((currentCent-lastCent).^2))>maxMoveDist ...
                    || currentIntRatio<minIntRatio
                
                rvsInds{kk-1}(nn) = NaN;
                
            end
        end
        
    end
    
end


%% --- Construct tracks from tracking and trace interuptions

lastPoint = numFrames;

tracks = cell(1,numNuc(lastPoint));

traceAtObjectRegister = cell(1,numFrames);

traceAtObjectRegister{numFrames} = 1:numNuc(lastPoint);

for nn = 1:numel(tracks)
    
    tracks{nn} = struct();
    tracks{nn}.time = tt_vector(lastPoint);
    tracks{nn}.timeInd = lastPoint;
    tracks{nn}.origObjInd = nn;
    tracks{nn}.centroid = centroid_cell{lastPoint}(nn);
    tracks{nn}.volume = nuc_vol_cell{lastPoint}(nn);
    tracks{nn}.intensity = ...
        nucInt_cell{lastPoint}{intChannel}(nn);
    tracks{nn}.cytoIntensity = ...
        cytoInt_cell{lastPoint}{intChannel}(nn);
    
    tracks{nn}.truncated = false;
    tracks{nn}.joint = [];
    
    traceAtObjectRegister{numFrames}(nn) = nn;
    
end

truncatedVec = false(1,numNuc(lastPoint));
jointVec = false(1,numNuc(lastPoint));

for kk = (lastPoint-1):-1:1
    
    kk
    
    if numNuc(kk) == 0
        % Make sure all tracks are truncated
        
        for tt = 1:numel(tracks)
            
            if ~truncatedVec(tt) && ~jointVec(tt)
                
                tracks{tt}.truncated = true;
                truncatedVec(tt) = true;
                
            end
            
        end
        
    else
        
        % positions of original objects detected in this frame, while this
        % vector holds the index of the trace associated with this original
        % object
        usedOrigInds = zeros(1,numNuc(kk));
        
        % Register that points from every original object in every frame to
        % the trace that it has been placed into
        traceAtObjectRegister{kk} = zeros(1,numNuc(kk));
        
        for tt = find(~truncatedVec & ~jointVec)
            
            if ~isempty(rvsInds{kk})
                if isfinite(rvsInds{kk}(tracks{tt}.origObjInd(end)))

                    connectInd = rvsInds{kk}(tracks{tt}.origObjInd(end));

                    % Check if another trace already has taken the object to
                    % connect, if yes then join this trace to that trace
                    joinInd = find(usedOrigInds==connectInd,1,'first');

                    if ~isempty(joinInd)
                        % Save time point and target trace that this trace was
                        % joint with
                        tracks{tt}.joint = [kk,joinInd];
                        jointVec(tt) = true;
                    else

                        usedOrigInds(tt) = connectInd;

                        traceAtObjectRegister{kk}(connectInd) = tt;

                        tracks{tt}.time = [tracks{tt}.time,tt_vector(kk)];
                        tracks{tt}.timeInd = [tracks{tt}.timeInd,kk];
                        tracks{tt}.origObjInd = ...
                            [tracks{tt}.origObjInd,connectInd];
                        tracks{tt}.centroid = ...
                            [tracks{tt}.centroid,centroid_cell{kk}(connectInd)];
                        tracks{tt}.volume = ...
                            [tracks{tt}.volume,nuc_vol_cell{kk}(connectInd)];
                        tracks{tt}.intensity = ...
                            [tracks{tt}.intensity,...
                            nucInt_cell{kk}{intChannel}(connectInd)];
                        tracks{tt}.cytoIntensity = ...
                            [tracks{tt}.cytoIntensity,...
                            cytoInt_cell{kk}{intChannel}(connectInd)];
                    end
                
                else
                    
                    tracks{tt}.truncated = true;
                    truncatedVec(tt) = true;
                    
                end
            end
            
        end
                        
        % Make new tracks for objects whose indices were not used
        origInds = 1:numNuc(kk);
        unusedInds = setdiff(origInds,usedOrigInds);
        
        addTracksCell = cell(1,numel(unusedInds));
        truncatedVecAddition = false(1,numel(unusedInds));
        jointVecAddition = false(1,numel(unusedInds));
        
        for jj = 1:numel(unusedInds)
            
            addInd = unusedInds(jj);
            
            addTrack = struct;
            
            addTrack.time = tt_vector(kk);
            addTrack.timeInd = kk;
            addTrack.origObjInd = addInd;
            addTrack.centroid = centroid_cell{kk}(addInd);
            addTrack.volume = nuc_vol_cell{kk}(addInd);
            addTrack.intensity = ...
                nucInt_cell{kk}{intChannel}(addInd);
            addTrack.cytoIntensity = ...
                cytoInt_cell{kk}{intChannel}(addInd);
            
            addTrack.truncated = false;
            addTrack.joint = [];
            
            addTracksCell{jj} = addTrack;
                        
            traceAtObjectRegister{kk}(addInd) = numel(tracks);
            
        end
        
        tracks = [tracks,addTracksCell];
        truncatedVec = [truncatedVec,truncatedVecAddition];
        jointVec = [jointVec,jointVecAddition];
        
    end
    
end


%% -- save track file

save([load_sourceDir,'TrackSet_',load_sourceFile],'tracks',...
    'maxMoveDist','tt_vector','numFrames',...
    'xmaxProj_cell','ymaxProj_cell','zmaxProj_cell',...
    'voxelSizeX','voxelSizeY','voxelSizeZ',...
    'rawStackSizeX','rawStackSizeY','rawStackSizeZ','binning')


%% -- Overview of all tracks

minTimePoints = 4; % Minimum number of time points a trace must be present
minIntRatioPeak = minIntRatio; % Minimum peak intensity ratio

figure(2)

clf

validTrackInds = ...
    cellfun(@(elmt)numel(elmt.time)>=minTimePoints,tracks) ...
    & cellfun(@(elmt)max(elmt.intensity./elmt.cytoIntensity),tracks) ...
    >=minIntRatioPeak ...
    & cellfun(@(elmt)max(elmt.volume),tracks)<=maxVol ...
    & cellfun(@(elmt)min(elmt.volume),tracks)>=minVol;

validTracks = tracks(validTrackInds);

for kk = 1:numel(validTracks)
    
    kk
    
    plot(validTracks{kk}.time,...
        validTracks{kk}.intensity./validTracks{kk}.cytoIntensity)
    
    hold on
    
end

hold off

xlabel('Time [min]')
ylabel('V_{nuc} [\mum^3]')


%% Extract sum of nuclear volumes for different stages

% -- set parameters and context

if inputFile == 1
    
    timeWindows = ...
        {[4,14],[20,29],[32,45],[48,60],[65,77],[81,96],...
        [102,122],[124,145],[148,182]};
    
    window_labels = {'32','64','128','256','512','1K',...
        'high','oblong','sphere'};
    
    % cell cycles corresponding to stages
    windowCellCycles = [5,6,7,8,9,10,11,12,13];
    
elseif inputFile == 2
    
    timeWindows = ...
        {[3,14],[18,28],[33,43],[48,58],[64,76],...
        [81,100],[108,138]};
    
    window_labels = {'64','128','256','512','1K',...
        'high','oblong'};
    
    % cell cycles corresponding to stages
    windowCellCycles = [6,7,8,9,10,11,12];

elseif inputFile == 3
    % still needs adjustment
    
        timeWindows = ...
        {[3,8],[17,24],[30,39]};
    
    window_labels = {'16','32','64'};
    
    % cell cycles corresponding to stages
    windowCellCycles = [4,5,6];

elseif inputFile == 4
    
    timeWindows = ...
        {[0,7],[12,21],[28,39],...
        [45,58],[66,86]};
    
    window_labels = {'128','256','512',...
        '1K','high'};
    
    % cell cycles corresponding to stages
    windowCellCycles = [7,8,9,10,11];
    
end


% --- execute analysis

numWindows = numel(timeWindows);

windowNumNuc = zeros(1,numWindows);
windowTotalVol = zeros(1,numWindows);
windowIndivVol = zeros(1,numWindows);
windowCytoVol = zeros(1,numWindows);

for ww = 1:numWindows
    
    ww
    
    minTime = timeWindows{ww}(1);
    maxTime = timeWindows{ww}(2);
    
    inclInd = tt_vector./60>=minTime & tt_vector./60<=maxTime;
    
    windowCytoVol(ww) = max(cyto_vol_vec(inclInd));
    
    pickTrackMask = false(size(validTracks));
    peakVols = zeros(size(validTracks));
    
    for kk = 1:numel(validTracks)
        
        thisTrack = validTracks{kk};
        
        traceTimes = thisTrack.time;
        
        traceInts = thisTrack.intensity./thisTrack.cytoIntensity;
                
        firstIntTime = traceTimes(find(traceInts>=minIntRatio,1,'last'))./60;
        lastIntTime = traceTimes(find(traceInts>=minIntRatio,1,'first'))./60;
                
        traceVols = thisTrack.volume;
        peakVols(kk) = max(traceVols);
        
        pickTrackMask(kk) =  ~(lastIntTime<minTime) && ~(firstIntTime>=maxTime);
        
    end
    
    windowNumNuc(ww) = sum(pickTrackMask);
    windowTotalVol(ww) = sum(peakVols.*pickTrackMask);
    windowIndivVol(ww) = median(peakVols(pickTrackMask));
    
end



%% --- Total cytoplasmic volume

[cytoVolPeaks,peakInds] = findpeaks(cyto_vol_vec,'MinPeakDistance',6);

figure(3)

plot(tt_vector(peakInds)./60,cytoVolPeaks,'ko')
hold on
plot(tt_vector./60,cyto_vol_vec,'b-')
hold off

ylim([0,inf])

xlabel('Time [min]')
ylabel('V_{cyto} [\mum^3]')

cytoVol = mean(cytoVolPeaks(1:5));

%% --- save tracking and volume results

if inputFile == 1
    
    saveFile = 'LiveImagingVolumes_10_G1.mat';
    
elseif inputFile == 2

    saveFile = 'LiveImagingVolumes_10_G2.mat';

elseif inputFile == 3

    saveFile = 'LiveImagingVolumes_5Feb_Early.mat';

elseif inputFile == 4

    saveFile = 'LiveImagingVolumes_5Feb_Late.mat';
        
end

saveTarget = fullfile(load_sourceDir,saveFile);

save(saveTarget,'windowCellCycles','windowNumNuc',...
    'windowTotalVol','windowIndivVol','windowCytoVol')
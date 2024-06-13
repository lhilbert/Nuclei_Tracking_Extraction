clear all

%% --- tracking parameters

maxMoveDist = 20; % in unit of micrometers

minIntRatio = 1.025; % minimum intensity increase in nucleus
% vs. surrounding cytoplasm

intChannel = 3; % Channel to get intensity ratios from

inputFile = 0; % 0 - new file, 1,2,{3,4} correspond to known data sets

overrideTimstamps = false; % rplace timestamps saved during acquisition?
deltat = 90; % time interval in seconds

%% --- Load raw analysis results
    
[load_sourceFile,load_sourceDir] = uigetfile('*.*');

thisPath = fullfile(load_sourceDir,load_sourceFile);

load(thisPath);

if overrideTimstamps
    % Overrides time data of stack itself, comment out when time info not
    % corrupted
    tt_vector = 0:deltat:deltat.*(numFrames-1);
end


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
    tracks{nn}.nucIntensity = cell(1,numChannels);
    tracks{nn}.cytoIntensity = cell(1,numChannels);
    for cc = 1:numChannels
        tracks{nn}.nucIntensity{cc} = ...
            nucInt_cell{lastPoint}{cc}(nn);
        tracks{nn}.cytoIntensity{cc} = ...
            cytoInt_cell{lastPoint}{cc}(nn);
    end
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
                        for cc = 1:numChannels
                            tracks{tt}.nucIntensity{cc} = ...
                                [tracks{tt}.nucIntensity{cc},...
                                nucInt_cell{kk}{cc}(connectInd)];
                            tracks{tt}.cytoIntensity{cc} = ...
                                [tracks{tt}.cytoIntensity{cc},...
                                cytoInt_cell{kk}{cc}(connectInd)];
                        end
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
            
            addTrack.nucIntensity = cell(1,numChannels);
            addTrack.cytoIntensity = cell(1,numChannels);
            
            for cc = 1:numChannels
                addTrack.nucIntensity{cc} = ...
                    nucInt_cell{kk}{cc}(addInd);
                addTrack.cytoIntensity{cc} = ...
                    cytoInt_cell{kk}{cc}(addInd);
            end
    
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
    'maxMoveDist','tt_vector','numFrames','intChannel',...
    'xmaxProj_cell','ymaxProj_cell','zmaxProj_cell',...
    'voxelSizeX','voxelSizeY','voxelSizeZ',...
    'rawStackSizeX','rawStackSizeY','rawStackSizeZ','binning',...
    'cyto_vol_vec')
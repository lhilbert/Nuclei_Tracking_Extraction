clear all

%% --- prepare containers

fileTargets = [1,2,3];

numTimeLapses = numel(unique(fileTargets));

tt_vec_cell = cell(1,numTimeLapses);
numNuc_vec_cell = cell(1,numTimeLapses);

cellCycles_cell = cell(1,numTimeLapses);
numNuc_cell = cell(1,numTimeLapses);
totalVol_cell = cell(1,numTimeLapses);
indivVol_cell = cell(1,numTimeLapses);
cytoVol_cell = cell(1,numTimeLapses);

singleNucleiVolumes_cell = cell(1,numTimeLapses);

timeWindow_cell = cell(1,numTimeLapses);

%% --- prepare time windows

for kk = 1:numTimeLapses
    
   if kk == 1
        
        timeWindows = ...
            {[3,7.5],[17,22.5],[32,37],[45,51],[62,66],[75,79.5],...
            [90,94.5],[106,112],[126,134],[151,159],[173,181],[215,225]};
        
        % cell cycles corresponding to stages
        windowCellCycles = [3,4,5,6,7,8,9,10,11,12,13,14];
        
    elseif kk == 2
        
        timeWindows = ...
            {[6,9],[21,24],[34,39],[47,52],[64,67],[78,80],...
            [95,98],[109,111],[129,134],[153,161],[183,191],[223,233]};
        
        % cell cycles corresponding to stages
        windowCellCycles = [3,4,5,6,7,8,9,10,11,12,13,14];
        
    elseif kk == 3
        
        timeWindows = ...
            {[4,9],[18,23],[33,38],[49,52],[62,66],[79,82],[94,97],...
            [109,113],[129,134],[154,162],[185,193],[225,235]};
        
        % cell cycles corresponding to stages
        windowCellCycles = [3,4,5,6,7,8,9,10,11,12,13,14];
        
    end

    timeWindow_cell{kk} = timeWindows; 
    cellCycles_cell{kk} = windowCellCycles;
end


%% --- extraction of volumes

for kk = 1:numTimeLapses

    kk
    
    inputFile = kk;
    
    if inputFile == 1
        % Load first timelapse results
        
        sourceDir = '/Users/hilbert/academia/VastenhouwLab/WholeEmbryoHistoneRatioAnalysis/LiveImagingData/';
        saveFile = 'TrackSet_TimeLapse_0001_G3_Subset.czi_AnalysisOutput.mat';
                
    elseif inputFile == 2
        
        sourceDir = '/Users/hilbert/academia/VastenhouwLab/WholeEmbryoHistoneRatioAnalysis/LiveImagingData/';
        saveFile = 'TrackSet_TimeLapse_0001_G2_Subset.czi_AnalysisOutput.mat';
        
    elseif inputFile == 3
        
        sourceDir = '/Users/hilbert/academia/VastenhouwLab/WholeEmbryoHistoneRatioAnalysis/LiveImagingData/';
        saveFile = 'TrackSet_TimeLapse_0001_G1_Subset.czi_AnalysisOutput.mat';
         
    end
    
    loadStruct = load(fullfile(sourceDir,saveFile));
    
    keepTrackFlags = loadStruct.inLimitsFlags & ~loadStruct.rejectFlags;
    
    includedTracks = loadStruct.tracks(keepTrackFlags);
    
    %% --- plot the time courses of number of nuclei
    
    tt_vector = loadStruct.tt_vector;
    numTimePoints = numel(tt_vector);
    
    numNuc = zeros(1,numTimePoints);
    
    parfor ll = 1:numTimePoints
        
        numNuc(ll) = sum( ...
            cellfun(@(elmt)ismember(ll,elmt.timeInd),includedTracks));
        
    end
    
    tt_vec_cell{kk} = tt_vector;
    numNuc_vec_cell{kk} = numNuc;
    
    % Extract sum of nuclear volumes for different stages
    
    % -- set parameters and context
    
    timeWindows = timeWindow_cell{kk};
    
    % --- execute analysis
    
    numWindows = numel(timeWindows);
    
    windowNumNuc = zeros(1,numWindows);
    windowTotalVol = zeros(1,numWindows);
    windowIndivVol = zeros(1,numWindows);
    windowCytoVol = zeros(1,numWindows);
    
    windowSingleNucleiVolumes = cell(1,numWindows);
        
    for ww = 1:numWindows
                
        minTime = timeWindows{ww}(1);
        maxTime = timeWindows{ww}(2);
        
        minInd = find(tt_vector./60>=minTime,1,'first');
        maxInd = find(tt_vector./60<=maxTime,1,'last');
        inclInd = minInd:1:maxInd;
        
        windowCytoVol(ww) = max(loadStruct.cyto_vol_vec(inclInd));
        
        pickTraceMask = false(size(includedTracks));
        peakVols = zeros(size(includedTracks));
        
        for tt = 1:numel(includedTracks)
            
            thisTrace = includedTracks{tt};
            
            traceTimes = thisTrace.time./60;
            
            pickTraceMask(tt) =  max(traceTimes)>=minTime ...
                && min(traceTimes)<=maxTime;
            
            if pickTraceMask(tt)
                
                thisTimeInds = thisTrace.timeInd;
                
                volStartInd = max([minInd,min(thisTimeInds)]);
                volEndInd = min([maxInd,max(thisTimeInds)]);
                
                mappedStartInd = find(thisTimeInds==volStartInd);
                mappedEndInd = find(thisTimeInds==volEndInd);
                
                volMinInd = min(mappedStartInd,mappedEndInd);
                volMaxInd = max(mappedStartInd,mappedEndInd);
                                
                traceVols = thisTrace.volume(volMinInd:volMaxInd);
                
                peakVols(tt) = max(traceVols);
                
            end
            
        end
        
        windowNumNuc(ww) = sum(pickTraceMask);
        windowTotalVol(ww) = sum(peakVols(pickTraceMask));
        windowIndivVol(ww) = mean(peakVols(pickTraceMask));
        
        windowSingleNucleiVolumes{ww} = peakVols(pickTraceMask);
        
    end
    
    % Store results for this data set
    numNuc_cell{kk} = windowNumNuc;
    totalVol_cell{kk} = windowTotalVol;
    indivVol_cell{kk} = windowIndivVol;
    cytoVol_cell{kk} = windowCytoVol;
    
    singleNucleiVolumes_cell{kk} = windowSingleNucleiVolumes;
    
end
    

%% --- Save volume data

save([sourceDir,'LiveImagingVolumes'],...
    'numTimeLapses','numNuc_cell','totalVol_cell','indivVol_cell',...
    'cytoVol_cell','timeWindow_cell','cellCycles_cell')

%% --- Plot volume data

plotStrings = {'k-o','r-s','b-^'};

% --- single frame nuclei counts

figure(2)

clf

for kk = 1:numTimeLapses

    subplot(numTimeLapses,1,kk)
        
    for ll = 1:numel(timeWindow_cell{kk})
    
        patch([timeWindow_cell{kk}{ll}(1),timeWindow_cell{kk}{ll}(2),...
            timeWindow_cell{kk}{ll}(2),timeWindow_cell{kk}{ll}(1)],...
            [1,1,1e4,1e4],...
            [0.7,0.7,0.7],'EdgeColor','none')
        
        hold on
        
    end
    
    plot(tt_vec_cell{kk}./60,...
        numNuc_vec_cell{kk},plotStrings{kk})
    
    hold off
    
    xlabel('Time [min]')
    ylabel('Nuclei')
    
    set(gca,'XLim',tt_vec_cell{kk}([1,end])./60,'YScale','log')
    
end


figure(1)

subplot(4,1,1)

for kk = 1:numTimeLapses
    
    plot(cellCycles_cell{kk},numNuc_cell{kk},...
        plotStrings{kk})
    
    hold on

end

hold off

set(gca,'XTick',[3,4,5,6,7,8,9,10,11,12,13,14],...
    'XTickLabel',{'8','16','32','64','128','256','512','1K','Hi',...
    'Obl','Sph','Dome'},...
    'XLim',[2.5,14.5])

xlabel('Cell cycle')
ylabel('Nuclei count')




subplot(4,1,2)

plotStrings = {'k-o','r-s','b-^'};

for kk = 1:numTimeLapses
    
    plot(cellCycles_cell{kk},totalVol_cell{kk},plotStrings{kk})
    
    hold on

end

hold off

set(gca,'XTick',[3,4,5,6,7,8,9,10,11,12,13,14],...
    'XTickLabel',{'8','16','32','64','128','256','512','1K','Hi',...
    'Obl','Sph','Dome'},...
    'XLim',[2.5,14.5])

xlabel('Cell cycle')
ylabel('V_{nuc} [\mum^3]')


subplot(4,1,3)

plotStrings = {'k-o','r-s','b-^'};

for kk = 1:numTimeLapses
    
    plot(cellCycles_cell{kk},cytoVol_cell{kk},plotStrings{kk})
    
    hold on

end

hold off

set(gca,'XTick',[3,4,5,6,7,8,9,10,11,12,13,14],...
    'XTickLabel',{'8','16','32','64','128','256','512','1K','Hi',...
    'Obl','Sph','Dome'},...
    'XLim',[2.5,14.5])

xlabel('Cell cycle')
ylabel('V_{cyto} [\mum^3]')






subplot(4,1,4)

plotStrings = {'k-o','r-s','b-^'};

for kk = 1:numTimeLapses
    
    refVol = mean(cytoVol_cell{kk}...
        (cellCycles_cell{kk}>=8 & cellCycles_cell{kk}<=10));
    
    plot(cellCycles_cell{kk},...
        totalVol_cell{kk}./refVol,...
        plotStrings{kk})
    
    hold on

end

hold off

set(gca,'XTick',[3,4,5,6,7,8,9,10,11,12,13,14],...
    'XTickLabel',{'8','16','32','64','128','256','512','1K','Hi',...
    'Obl','Sph','Dome'},...
    'XLim',[2.5,14.5])

xlabel('Cell cycle')
ylabel('%V_{nuc}')

%% --- plot single nuclei distributions

stageNames = {'8','16','32','64','128','256','512','1K','Hi',...
    'Obl','Sph','Dome'};

figure(3)

clf

for kk = 1:numTimeLapses
    
    for ww = 1:numel(singleNucleiVolumes_cell{kk})

        subplot(numel(singleNucleiVolumes_cell{kk}),1,ww)

        plot_color = [1,0.6,0.6].*(kk-1) ...
            ./(numTimeLapses-1);
        
        xx_support = linspace(0,1e4,500);
        
        [yy,xx] = ksdensity(singleNucleiVolumes_cell{kk}{ww},...
            xx_support,...
            'Bandwidth',1e3);
        
%         [yy,xx] = hist(singleNucleiVolumes_cell{kk}{ww});
        
        plot(xx,yy,'k-','Color',plot_color)
        
        hold on
        
        title(stageNames{ww})
        
    end
    
    
end

for ww = 1:numel(singleNucleiVolumes_cell{kk})
    
    subplot(numel(singleNucleiVolumes_cell{kk}),1,ww)
    
    hold off
    
    ylabel('Bin count')
    
end

xlabel('Nulceus volume [\mum^3]')

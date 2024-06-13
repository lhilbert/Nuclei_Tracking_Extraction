clear all

%% --- values from mass spectrometry experiment

MassSpecStages = ...
    {'1-cell','8-cell','64-cell','256-cell',...
    '1000-cell','High','Oblong','Sphere'};

MassSpecVals = [4963.278739,...
5752.397177,...
7179.586696,...
9760.507518,...
13629.36916,...
13962.80005,...
15836.81046,...
16000]; % get exact last value from Shai

MassSpecCellCycles = [0,3,6,8,10,11,12,13];


%% --- load experimental data from live embryo experiments

numFiles = 4;
fileTargets = [1,2,3,3];

numTimeLapses = numel(unique(fileTargets));

cellCycles_cell = cell(1,numTimeLapses);
numNuc_cell = cell(1,numTimeLapses);
totalVol_cell = cell(1,numTimeLapses);
indivVol_cell = cell(1,numTimeLapses);
cytoVol_cell = cell(1,numTimeLapses);

for kk = 1:numTimeLapses
    
    cellCycles_cell{kk} = [];
    numNuc_cell{kk} = [];
    totalVol_cell{kk} = [];
    indivVol_cell{kk} = {};
    cytoVol_cell{kk} = [];

end

for kk = 1:numFiles

    inputFile = kk;
    
    if inputFile == 1
        % Load first timelapse results
        
        sourceDir = '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/';
        saveFile = 'LiveImagingVolumes_10_G1.mat';
                
    elseif inputFile == 2
        
        sourceDir = '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/';
        saveFile = 'LiveImagingVolumes_10_G2.mat';
        
    elseif inputFile == 3
        
        sourceDir = ...
            '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/PCNASubsets_5Feb2016/';
        saveFile = 'LiveImagingVolumes_5Feb_Early.mat';
        
    elseif inputFile == 4
        
        sourceDir = ...
            '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/PCNASubsets_5Feb2016/';
        saveFile = 'LiveImagingVolumes_5Feb_Late.mat';
        
    end
    
    loadStruct = load(fullfile(sourceDir,saveFile));
    
    cellCycles_cell{fileTargets(kk)} = ...
        [cellCycles_cell{fileTargets(kk)}, ...
        loadStruct.windowCellCycles];
    numNuc_cell{fileTargets(kk)} = ...
        [numNuc_cell{fileTargets(kk)}, ...
        loadStruct.windowNumNuc];
    totalVol_cell{fileTargets(kk)} = ...
        [totalVol_cell{fileTargets(kk)},...
        loadStruct.windowTotalVol];
    indivVol_cell{fileTargets(kk)} = ...
        {indivVol_cell{fileTargets(kk)}, ...
        loadStruct.windowIndivVol};
    cytoVol_cell{fileTargets(kk)} = ...
        [cytoVol_cell{fileTargets(kk)},...
        loadStruct.windowCytoVol];
    
end
    

figure(1)

subplot(1,3,1)

plotStrings = {'k-o','r-s','b-^'};

for kk = 1:numTimeLapses
    
    plot(cellCycles_cell{kk},numNuc_cell{kk},plotStrings{kk})
    
    hold on

end

hold off

xlabel('Cell cycle')
ylabel('Nuclei count')




subplot(2,3,2)

plotStrings = {'k-o','r-s','b-^'};

for kk = 1:numTimeLapses
    
    plot(cellCycles_cell{kk},totalVol_cell{kk},plotStrings{kk})
    
    hold on

end

hold off

xlabel('Cell cycle')
ylabel('V_{nuc} [\mum^3]')


subplot(2,3,5)

plotStrings = {'k-o','r-s','b-^'};

for kk = 1:numTimeLapses
    
    plot(cellCycles_cell{kk},cytoVol_cell{kk},plotStrings{kk})
    
    hold on

end

hold off

xlabel('Cell cycle')
ylabel('V_{cyto} [\mum^3]')






subplot(1,3,3)

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

xlabel('Cell cycle')
ylabel('%V_{nuc}')





%% --- load experimental data from fixed embryo IF experiments
    





%% --- calculate nuclear concentration





% cNuc = zeros(1,numAssmt);
%
%
%
%
% for kk = 1:numAssmt
%
%     kk
%     beginTime = volAssmntTimes(kk)-spanBkwd
%     endTime = volAssmntTimes(kk)+spanFwd
%
%     inclMask = cellfun(@(elmt) ...
%         sum(arrayfun(@(tt)tt>=beginTime&&tt<=endTime,elmt.time-tt_vector(1)))>0,...
%         validTraces);
%
%     numNucAssmt(kk) = sum(inclMask);
%
%     nucVolAssmt = sum(cellfun(@(elmt)max(elmt.volume),validTraces(inclMask)));
%
%     volFractionRatioAssmt(kk) = ...
%         nucVolAssmt./(cytoVolPlot(timeInds(kk)));
%
%     intRatioAssmt(kk) = ...
%         median(...
%         cellfun(@(elmt)max((elmt.intensity-I_offset)...
%         ./(elmt.cytoIntensity-I_offset)),...
%         validTraces(inclMask)));
%
%     medianNucVolAssmt(kk) = median(cellfun(@(elmt)max(elmt.volume), ...
%         validTraces(inclMask)));
%
% end
%
% clf
%
% subplot(3,1,1)
%
% plot(stageSupport,numNucAssmt,'ko')
% xlabel('Stage')
% ylabel('Cell number')
%
% set(gca,'XTick',stageSupport,...
%     'XTickLabel',stageLabels)
%
% subplot(3,1,2)
%
% plot(stageSupport,100.*volFractionRatioAssmt,'ko')
% xlabel('Stage')
%
% ylabel('Vol% Nucleus')
%
% set(gca,'XTick',stageSupport,...
%     'XTickLabel',stageLabels)
%
%
% subplot(3,1,3)
%
% plot(stageSupport,intRatioAssmt,'ko')
% xlabel('Stage')
%
% hold on
%
% Htot = interp1(MassSpecCellCounts,MassSpecVals,...
%     numNucAssmt);
% Vtot = cytoVol(timeInds);
% RR = intRatioAssmt;
% ff = volFractionRatioAssmt;
% Hbound_single = 1.0;
% Vnuc_single = medianNucVolAssmt;
%
% c_nuc_free = Htot./Vtot.*RR./(1-ff+RR.*ff)-Hbound_single./Vnuc_single;
%
% % plot(stageSupport,c_nuc_free,'k-',...
% %     'Color',[0.6,0.6,0.6])
%
% hold off
%
% ylabel('I_{nuc}/I_{cyto}')
%
% set(gca,'XTick',stageSupport,...
%     'XTickLabel',stageLabels)
%

clear all

plotIndividualTracks = false;

[load_sourceFile,load_sourceDir] = uigetfile('*.*');
thisPath = fullfile(load_sourceDir,load_sourceFile);

data = load(thisPath);

rejectFlags = data.rejectFlags;
inLimitsFlags = data.inLimitsFlags;

useInds = inLimitsFlags & ~rejectFlags;

tracks = data.tracks(useInds);

numTracks = numel(tracks);

numChannels = numel(tracks{1}.nucIntensity);

%% Run through time course and get average N/C ratio for all channels

numFrames = data.numFrames;
frameTimes = data.tt_vector;

NC_ratio_storage_cell = cell(numChannels,numFrames);

for kk = 1:numFrames
    
    disp(numFrames-kk);
    
    for cc = 1:numChannels
        
        NC_ratio_storage_cell{cc,kk} = [];
        
    end
    
    for ll = 1:numTracks
       
        thisTrack = tracks{ll};
        
        [correctTimeFlag,correctTimeInd] = ...
            ismember(kk,thisTrack.timeInd);
        
        if correctTimeFlag
           
            for cc = 1:numChannels
            
                nucInt = thisTrack.nucIntensity{cc}(correctTimeInd);
                cytoInt = thisTrack.cytoIntensity{cc}(correctTimeInd);
            
                NC_ratio_storage_cell{cc,kk} = ...
                    [NC_ratio_storage_cell{cc,kk},nucInt./cytoInt];
                
            end
            
        end
        
    end
    
end

numNuc = cellfun(@(elmt)numel(elmt),NC_ratio_storage_cell(1,:));
NC_ratio_storage_mean = ...
    cellfun(@(elmt)mean(elmt,'omitnan'),NC_ratio_storage_cell);
NC_ratio_storage_stdev = ...
    cellfun(@(elmt)std(elmt,'omitnan'),NC_ratio_storage_cell);




%% plotting

tt_vector = data.tt_vector;

figure(1)

clf

for cc = 1:numChannels
    
    subplot(1+2.*numChannels,1,(1+2.*(cc-1)):(2+2.*(cc-1)))
    
    if plotIndividualTracks
        
        for tt = 1:numTracks
            
            
            lh = plot(tracks{tt}.time./60, ...
                tracks{tt}.nucIntensity{cc}./tracks{tt}.cytoIntensity{cc},'k-');
            set(lh,'Color',[0,0,0,0.075])
            
            hold on
            
        end
        
    end
    
    plot(tt_vector./60,NC_ratio_storage_mean(cc,:),'k-')
    hold on
    plot(tt_vector./60,NC_ratio_storage_mean(cc,:)...
        +NC_ratio_storage_stdev(cc,:),'k-',...
        'Color',[0.6,0.6,0.6])
    plot(tt_vector./60,NC_ratio_storage_mean(cc,:)...
        -NC_ratio_storage_stdev(cc,:),'k-',...
        'Color',[0.6,0.6,0.6])
    hold off
    
    set(gca,'XLim',tt_vector([1,end])./60,'XTickLabel',[])
    xlabel('')
    ylabel(sprintf('N/C Channel %d',cc))
    
end

subplot(1+2.*numChannels,1,1+2.*numChannels)
plot(tt_vector./60,numNuc,'k-')
xlabel('Time [min]')
ylabel('Nuclei')

%% Write to .csv file

saveFile = fullfile(...
    load_sourceDir,[load_sourceFile(1:end-4),'_export.csv']);

columnLabels = cell(2.*(cc+1),1);

columnLabels{1} = 'Time(min)';
columnLabels{2} = 'numNuc';

writeArray = zeros(numFrames,2.*(cc+1));
writeArray(:,1) = tt_vector./60;
writeArray(:,2) = numNuc;

for cc = 1:numChannels
   
    columnLabels{1+2.*cc} = sprintf('Mean(Ch%d)',cc);
    columnLabels{2+2.*cc} = sprintf('STdev(Ch%d)',cc);
    
    writeArray(:,1+2.*cc) = NC_ratio_storage_mean(cc,:);
    writeArray(:,2+2.*cc) = NC_ratio_storage_stdev(cc,:);

end


fid = fopen(saveFile,'w');
for index = 1:(2.*(cc+1)-1)    
    fprintf(fid, '%s,', columnLabels{index});
end 
fprintf(fid, '%s\n', columnLabels{end});
fclose(fid);

dlmwrite(saveFile,writeArray,'-append', 'delimiter', ',');


% csvwrite(saveFile,writeArray,1,5);


%% -- write individual values to comma separated file


for cc = 1:numChannels

    saveFile = fullfile(...
        load_sourceDir,[load_sourceFile(1:end-4),...
        sprintf('_export_indNuclei_Ch%d.csv',cc)]);

    fid = fopen(saveFile,'w');
    
    for kk = 1:numFrames
    
        thisFrameRatios = sort(NC_ratio_storage_cell{cc,kk},'descend');
        
        thisFrameRatios = thisFrameRatios(~isnan(thisFrameRatios));
        
        fprintf(fid,'%6.6f,',thisFrameRatios);
        fprintf(fid,'\n');
        
    end

    fclose(fid);
    
end

% write out time stamps

saveFile = fullfile(...
    load_sourceDir,[load_sourceFile(1:end-4),...
    sprintf('_export_indNuclei_timestamps.csv')]);

fid = fopen(saveFile,'w');

for kk = 1:numFrames
        
    fprintf(fid,'%6.6f\n',tt_vector(kk)./60);
    
end

fclose(fid);
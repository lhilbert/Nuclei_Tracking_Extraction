clear all

%% --- rejection by intensity ratio

minRatio = 1.025; % minimum intensity increase in nucleus
% vs. surrounding cytoplasm

intChannel = 3; % Channel to get intensity ratios from

%% --- Load raw analysis results

[sourceFile,sourceDir] = uigetfile('*.*');
    
thisPath = fullfile(sourceDir,sourceFile);

load(thisPath);


%% --- Plot max projections

for ff = 1:numFrames
   
    disp(nucInt_cell{ff}{intChannel}./cytoInt_cell{ff}{intChannel})
    
    % Reject segmented objects according to criteria
    
    inLimFlags = ...
        (nucInt_cell{ff}{intChannel}...
        ./cytoInt_cell{ff}{intChannel})>=minRatio;
    
    
    figure(2)
    
    clf
      
   xxCoords = cellfun(@(elmt)elmt(1),centroid_cell{ff}(inLimFlags));
   yyCoords = cellfun(@(elmt)elmt(2),centroid_cell{ff}(inLimFlags));
   zzCoords = cellfun(@(elmt)elmt(3),centroid_cell{ff}(inLimFlags));
   
   subplot(2,3,1)
    
   imagesc([0,floor(rawStackSizeX./binning).*voxelSizeX],...
       [0,floor(rawStackSizeY./binning).*voxelSizeY],...
       zmaxProj_cell{ff})
   
   hold on
   
   plot(xxCoords,yyCoords,'r+')
   
   xlabel('x [\mum]')
   ylabel('y [\mum]')
   
   colormap(flipud(gray))
   
   axis equal
   axis tight
   

   
   subplot(2,3,2)
    
   imagesc([0,floor(rawStackSizeY./binning).*voxelSizeY],...
       [0,floor(rawStackSizeZ./binning).*voxelSizeZ],...
       xmaxProj_cell{ff}.')

   
   hold on
   
   plot(yyCoords,zzCoords,'r+')
   
   xlabel('y [\mum]')
   ylabel('z [\mum]')
   
   colormap(flipud(gray))
   
   axis equal
   axis tight

   set(gca,'YDir','normal')
   
   
   subplot(2,3,3)
    
   imagesc([0,floor(rawStackSizeX./binning).*voxelSizeX],...
       [0,floor(rawStackSizeZ./binning).*voxelSizeZ],...
       ymaxProj_cell{ff}.')

   
   hold on
   
   plot(xxCoords,zzCoords,'r+')
   
   xlabel('y [\mum]')
   ylabel('z [\mum]')
   
   colormap(flipud(gray))
   
   axis equal
   axis tight

   set(gca,'YDir','normal')
   
   
   
   
   subplot(2,3,4)
    
   imagesc([0,floor(rawStackSizeX./binning).*voxelSizeX],...
       [0,floor(rawStackSizeY./binning).*voxelSizeY],...
       zmaxProj_cell{ff})
   
   xlabel('x [\mum]')
   ylabel('y [\mum]')
   
   colormap(flipud(gray))
   
   axis equal
   axis tight
   

   
   subplot(2,3,5)
    
   imagesc([0,floor(rawStackSizeY./binning).*voxelSizeY],...
       [0,floor(rawStackSizeZ./binning).*voxelSizeZ],...
       (xmaxProj_cell{ff}).')

   xlabel('y [\mum]')
   ylabel('z [\mum]')
   
   colormap(flipud(gray))
   
   axis equal
   axis tight

   set(gca,'YDir','normal')
   
   
   subplot(2,3,6)
    
   imagesc([0,floor(rawStackSizeX./binning).*voxelSizeX],...
       [0,floor(rawStackSizeZ./binning).*voxelSizeZ],...
       (ymaxProj_cell{ff}).')

   xlabel('y [\mum]')
   ylabel('z [\mum]')
   
   colormap(flipud(gray))
   
   axis equal
   axis tight

   set(gca,'YDir','normal')
   
   waitforbuttonpress;
   
end
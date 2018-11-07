%% Create tables with n value data 
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'HbEmpty';'Kr2xProxEmpty';'KrDist17C';'Kr2xDist32C';'KrBoth17C'};    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

% Save number of nuclei in each AP bin for things calculated based on
% individual spots 
Filename=[DropboxFolder filesep 'Constructs' filesep 'BurstSize.mat'];
    load(Filename);
    SpotTable=[];
for cc=1:length(AvgAmpAllAP) 
    SpotTable=[SpotTable;(sum(~isnan(AvgAmpAllAP(cc).AllAmps)))];
end
SpotTable=num2cell(SpotTable);
for cc=1:length(ConstructList)
    SpotTable{cc,42}=ConstructList{cc};
end


FileLocation=[DropboxFolder filesep 'Constructs' filesep 'SingleSpotsTable.mat'];
save(FileLocation, 'SpotTable');
xlswrite([DropboxFolder filesep 'Constructs' filesep 'SingleSpotsTable'],SpotTable);
    
% Save # of embryos associated with each AP position for single spot data
% for each construct 
AllEmbTable=[];
for cc=1:length(AvgAmpAllAP)
    EmbTable=[];
    for ee=1:length(AvgAmpAllAP(cc).EmbryoAmp)
        EmbTable=[EmbTable; (~isnan(AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp(:)))'];
    end
    AllEmbTable=[AllEmbTable; sum(EmbTable,1)];
end

AllEmbTable=num2cell(AllEmbTable);
for cc=1:length(ConstructList)
    AllEmbTable{cc,42}=ConstructList{cc};
end
FileLocation2=[DropboxFolder filesep 'Constructs' filesep 'SingleSpotsEmbTable.mat'];
save(FileLocation2, 'AllEmbTable');
xlswrite([DropboxFolder filesep 'Constructs' filesep 'SingleSpotsEmbTable'],AllEmbTable);


%% Save information per AP bin on # of nuclei for 2 spot data
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C'}

Filename2=[DropboxFolder filesep 'Constructs' filesep 'AllInterNoise.mat'];
    load(Filename2);
    
NucleiTable=[];
for cc=1:length(WholeNoise)
    NucleiTable=[NucleiTable; (sum(~isnan(WholeNoise(cc).AllCoVarNoise)))];
end
NucleiTable=num2cell(NucleiTable);
for cc=1:length(ConstructList)
    NucleiTable{cc,42}=ConstructList{cc};
end
FileLocation=[DropboxFolder filesep 'Constructs' filesep 'BothSpotNucleiTable.mat'];
save(FileLocation, 'NucleiTable');
xlswrite([DropboxFolder filesep 'Constructs' filesep 'BothSpotNucleiTable'],NucleiTable);

% By Embryo

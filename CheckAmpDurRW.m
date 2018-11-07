%% check if # duration == # amplitude 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'Kr2xProxEmpty';'KrBoth';'KrBothEmpty'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
for cc=1:length(ConstructList)
    Data=LoadMS2SetsCS(ConstructList{cc});
    NEmbryos = length(Data);
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        BurstProps=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        load(BurstProps);
        for bb=1:length(BurstProperties)
            if length(BurstProperties(bb).BurstAmplitude) ~= length(BurstProperties(bb).Duration)
                disp(['Mismatch at',PrefixName,'spot', num2str(bb)]);
                return
            end
            if length(BurstProperties(bb).BurstAmplitude) ~=BurstProperties(bb).NBursts
                disp(['Mismatch # bursts at', PrefixName, 'spot', num2str(bb)]);
                return
            end
        end
    end
end

%% Run SlopeBurstCalling for all data sets
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'Kr2xProxEmpty';'KrBoth';'KrBothEmpty';'KrDist32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
for cc=1:length(ConstructList)
    Data=LoadMS2SetsCS(ConstructList{cc});
    NEmbryos = length(Data);
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        CompPars=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat']
        load(CompPars);
        Spotz=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat']
        load(Spotz);
        CheckSlopeBurstRW
        %FileName=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        %save(FileName, 'BurstProperties');
        clear BurstProperties PrefixName CompiledParticles FileName CompPars
    end
end
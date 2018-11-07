%% Correlation of spatial and temporal noise to total mRNA production 
% Load the saved temporal and spatial noise tables and the total mRNA
% structure
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

Filename1=[DropboxFolder filesep 'Constructs' filesep 'AllTemporalNoise'];
load(Filename1);
Filename2=[DropboxFolder filesep 'Constructs' filesep 'AllSpatialNoise'];
load(Filename2);
Filename3=[DropboxFolder filesep 'Constructs' filesep 'MeanSpatialNoise'];
load(Filename3);
Filename4=[DropboxFolder filesep 'Constructs' filesep 'AllTotalmRNAProd'];
load(Filename4);

%Build matrices to do correlation coefficient. 
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
for cc=1:length(ConstructList)
    for ee=1:length(WholeTimeTable(cc).Embryo);
        EmbryoTimeComp=[];
        EmbCorr=[];
        for bb=1:41
            EmbryoTimeComp(bb,1)=WholeTimeTable(cc).Embryo(ee).EmbryoAvg(bb);
            EmbryoTimeComp(bb,2)=AvgProdAllAP(cc).EmbryosProd(ee).MeanProd(bb);
        end
        EmbCorr=(corrcoef(EmbryoTimeComp,'Rows','complete'));
        CorrValues(ee,cc)=EmbCorr(1,2);
    end
end



 



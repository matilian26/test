%% calculate average mRNA produced per nucleus for each construct 
%load constructs
ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Data=Data(1:4);    %limit to 4 embryos for now since only have that many with 2xDist
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    MinParticles=5; MinEmbryos=1;
    %run integratemRNA
    [TotalProd,TotalProdError,TotalProdN,MeanTotalProd,SDTotalProd,SETotalProd] =IntegratemRNA(Data,MinParticles,MinEmbryos);
    TotalProdStruct(cc).TotalProd=TotalProd;
    TotalProdStruct(cc).TotalProdError=TotalProdError;
    TotalProdStruct(cc).MeanTotalProd=MeanTotalProd;
    TotalProdStruct(cc).SETotalProd=SETotalProd;
end

APtoUse=input('Which AP bin to look at?');
EggLength=APbinID(APtoUse)*100;
figure
for cc=1:length(TotalProdStruct)
    errorbar(cc,TotalProdStruct(cc).MeanTotalProd(APtoUse,14),TotalProdStruct(cc).SETotalProd(APtoUse,14),'o')
hold on
end
xlabel('Construct')
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
ylabel('Total nascent mRNA');
title(['Mean mRNA production',' ',num2str(EggLength),'% egg length']); 

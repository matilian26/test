%% Compare all bursting variables 
ConstructList= {'KrDist','KrProx','KrBothSep','KrDistEmpty','KrProxEmpty','KrDistDuplicN','KrProxDuplic','KrBoth'} %{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
%Count for each construct
AvgmRNAProd=[];
AllVarsList=[];
AllProdList=[];
for cc=1:length(ConstructList)
    ConVarsList=[];
    ConProdList=[];
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
     NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    
    for ee=1:NEmbryos
        EmbVarsList=[];
        EmbProdList=[];
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end 
        load(filename);
        for aa=1:length(APbinID)
            InfoAP=find([BurstProperties.APBin]==APbinID(aa));
            VarsList=[];
            ProdList=[];
            if isempty(InfoAP)
                VarsList(1,[1:6])=nan;
                ProdList=nan;
            else
            for bb=1:length(InfoAP)
                VarsList(bb,1)=nanmean([BurstProperties(InfoAP(bb)).BurstAmplitude]);
                VarsList(bb,2)=nanmean([BurstProperties(InfoAP(bb)).Duration]);
                VarsList(bb,3)=nanmean([BurstProperties(InfoAP(bb)).Frequency]);
                VarsList(bb,4)=[BurstProperties(InfoAP(bb)).FirstTimeOn];
                VarsList(bb,5)=cc;
                VarsList(bb,6)=aa;
                ProdList(bb)=[BurstProperties(InfoAP(bb)).TotalmRNA];
            end
            end
            EmbVarsList=[EmbVarsList; VarsList];
            EmbProdList=[EmbProdList, ProdList];
        end
        ConVarsList=[ConVarsList;EmbVarsList];
        ConProdList=[ConProdList, EmbProdList];
    end
    AllVarsList=[AllVarsList; ConVarsList];
    AllProdList=[AllProdList, ConProdList];
end

    %% Calculate autocorrelation times   
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info   for right now just going to use nc14
ncUse=input('Want to only use nc14?','s');
NumberNuclei=nan(20,41,[length(ConstructList)]);
% go through each embryo of each construct
for cc=1:length(ConstructList)
     Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        else
            Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelation.mat'];
        end
        load(Filename);
        APstuff=[SpotDiff.APBin];
        APAutoVal=[];
        for aa=1:length(APbinID)
            APsubset=[];
            APSpots=[];
            APsubset=SpotDiff(APstuff==APbinID(aa));
            if isempty(APsubset)
                    AutoVal=[nan];
            else
                    AutoVal=nan(length(APsubset),15)
            for bb=1:length(APsubset)
                    for TimeLag=1:30
                TopMult=[];
                    Bottom=[];
                for ss=1:length([APsubset(bb).SpotOne])
                    if length(APsubset(bb).SpotOne) >= (TimeLag +ss)
                    TopMult(ss)=((APsubset(bb).SpotOne(ss)-nanmean([APsubset(bb).SpotOne]))*(APsubset(bb).SpotOne(ss+TimeLag)-nanmean([APsubset(bb).SpotOne])));
                    Bottom(ss)=(APsubset(bb).SpotOne(ss)-nanmean([APsubset(bb).SpotOne]))^2;
                    else
                        TopMult(ss)=nan;
                        Bottom(ss)=nan;
                    end
                end
                    Bottom=nansum(Bottom);
                    AutoVal(bb,TimeLag)=((nansum(TopMult))/Bottom);

                end
                end
            end
                AllAutoVals(cc).Embryo(ee).APbin(aa).AutoVal=AutoVal;    
                end
            end
            end
            
        

    %% Plotting autocorrelation 
    APtoUse=input('Which AP bin to use?');
    EggLength=APbinID.*100;
    for cc=1:length(ConstructList)
        figure
        EmbAutoVals=[];
        for ee=1:length(AllAutoVals(cc).Embryo)
            if length(AllAutoVals(cc).Embryo(ee).APbin(APtoUse).AutoVal) > 1
            EmbAutoVals=[EmbAutoVals;[AllAutoVals(cc).Embryo(ee).APbin(APtoUse).AutoVal]]
            end
        end
        plot(nanmean(EmbAutoVals));
        title(ConstructList{cc});
    end
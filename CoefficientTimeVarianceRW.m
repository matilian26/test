%% calculate coefficient of variance for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info   for right now just going to use nc14
ncUse=input('Want to only use nc14?','s');

% go through each embryo of each construct
for cc=1:length(ConstructList)
     Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    APTable=[];
    TimeTable=[];
    for ee=1:NEmbryos
        
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        Filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            Filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end
        load(Filename);
        %if ncUse=='y'
%         nc_number=[CompiledParticles.nc];
%         CompiledParticles_14=CompiledParticles(nc_number==14);
        for bb=1:length(BurstProperties)
            TimeLength=[];
            AvgFluo=[];
            VarPart=[];
            SumVarPart=[];
            TimeLength=length(BurstProperties(bb).FluoFrames);
            AvgFluo=nanmean([BurstProperties(bb).Fluo]);
            for tt=1:TimeLength
                VarPart(tt)=((BurstProperties(bb).Fluo(tt) - AvgFluo)^2);
            end
            SumVarPart=sum(VarPart);
            TimeTable(bb,ee)=(sqrt((1/TimeLength)*(SumVarPart))/AvgFluo);
        end
    end
    TimeTable(TimeTable==0)=nan;
    WholeTimeTable(cc).TimeTable=TimeTable;
end
        
        
%% Visualizing 
DistalColor=[8 180 238] ./ 255;
DoubDistColor=[1 17 181] ./ 255;
ProxColor=[251 230 60] ./ 255;
DoubProxColor=[251 190 80] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothColor=[12 195 82] ./ 255;

Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DoubDistColor;
Colors(5).Color=DoubProxColor;
Colors(6).Color=BothColor;

for cc=1:length(ConstructList)
    %for ee=1:size(WholeTimeTable(cc).TimeTable,2)
    ConTable=[];
    
    ConTable=[WholeTimeTable(cc).TimeTable];
    for bb=1:numel(ConTable)
        AllCons(bb,cc)=ConTable(bb);
    end
end
AllCons(AllCons==0)=nan;
[p,tbl,stats]=anova1(AllCons);
        

figure
AvgAllcons=[];
for cc=1:length(ConstructList)
    MeanTimeVar=[];
    MeanTimeVar=nanmean(WholeTimeTable(cc).TimeTable);
    for ee=1:length(MeanTimeVar)
    AvgAllcons(ee,cc)=MeanTimeVar(ee);
    end
    for ee=1:length(MeanTimeVar)
        plot(cc, MeanTimeVar(ee), 'o', 'Color', Colors(cc).Color)
        hold on 
    end
end
AvgAllcons(AvgAllcons==0)=nan;
boxplot(AvgAllcons,'Colors','k');

xticks([1:6]);
xticklabels({'Distal', 'Proximal', 'Both Sep', '2x Distal', '2x Proximal', 'Both'});
xlabel('Construct');
ylabel('Mean coefficient of variation');
title('Noise across time');

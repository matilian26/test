%% calculate average burst frequency for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');

%Count for each construct
AvgNBursts=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgNBurstsCon=[];
    firsttime=1;
    ConNBurstsAllAP=[];
    ConNBurstsSE=[];
    ConNBurstsSD=[];
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end 
        load(filename);
        NumberBursts(ee,cc)=length([BurstProperties.Duration]);
        %seperate out by AP bin
        NBurstsAllAP=[];
        for aa=1:length(APbinID)
            NBurstsAP=[];
            NumberBurstsAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(NumberBurstsAP)
                NBurstsAP=[NBurstsAP; nan];
            else
            for bb=1:length(NumberBurstsAP)
                NBurstsAP=[NBurstsAP;[BurstProperties(NumberBurstsAP(bb)).NBursts]'];  %put all durations at a given AP value in a column going down
            end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(NBurstsAP)
                NBurstsAllAP(bb,aa,ee)=NBurstsAP(bb);
            end
            NBurstsAllAP(NBurstsAllAP==0)=nan;
            
            NumberBurstsSD(ee,aa,cc)=nanstd(NBurstsAllAP(:,aa,ee));
        NumberBurstsSE(ee,aa,cc)=NumberBurstsSD(ee,aa,cc)/sqrt(length(NumberBurstsAP));             
            clear NumberBurstsAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(NBurstsAllAP,3)
            ConNBurstsAllAP=[ConNBurstsAllAP; NBurstsAllAP(:,:,bb)];
        end
        for bb=1:size(NumberBurstsSD,3)
            ConNBurstsSD=[ConNBurstsSD; NumberBurstsSD(:,:,bb)];
        end
        for bb=1:size(NumberBurstsSE,3)
            ConNBurstsSE=[ConNBurstsSE;NumberBurstsSE(:,:,bb)];
        end
        AvgNBurstsAllAP(cc).EmbryoBursts(ee).AvgNBursts=nanmean(NBurstsAllAP(:,:,ee));
        AvgNBurstsAllAP(cc).EmbryoBursts(ee).SE=NumberBurstsSE(ee,:,cc);
        AvgNBurstsAllAP(cc).EmbryoBursts(ee).SD=NumberBurstsSD(ee,:,cc);
    end
        AvgNBurstsAllAP(cc).AvgNBursts=nanmean(ConNBurstsAllAP);  %Avg duration of all embryos of a construct
        AvgNBurstsAllAP(cc).NBurstSD=nanmean(ConNBurstsSD);
        AvgNBurstsAllAP(cc).NBurstSE=nanmean(ConNBurstsSE);
        AvgNBurstsAllAP(cc).AllNBursts=[ConNBurstsAllAP];

end

%%Plot mean of each embryo at specific AP position for each construct
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for ee=1:length(AvgNBurstsAllAP(cc).EmbryoBursts)
        APBoxing(ee,cc)=AvgNBurstsAllAP(cc).EmbryoBursts(ee).AvgNBursts(APtoUse);
    end
end
APBoxing(APBoxing==0)=nan;
figure
DistalColor=[8 180 238] ./ 255;
DoubDistColor=[1 17 181] ./ 255;
ProxColor=[251 230 60] ./ 255;
DoubProxColor=[251 190 80] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothColor=[12 195 82] ./ 255;
Colors=[BothColor; DistalColor; ProxColor; BothSepColor; DoubDistColor; DoubProxColor];

    for ee=1:length(AvgNBurstsAllAP(1).EmbryoBursts)
    plot(1,AvgNBurstsAllAP(1).EmbryoBursts(ee).AvgNBursts(APtoUse),'o','LineWidth',1.5,'Color',DistalColor)
    hold on
    end
    for ee=1:length(AvgNBurstsAllAP(2).EmbryoBursts)
    plot(2,AvgNBurstsAllAP(2).EmbryoBursts(ee).AvgNBursts(APtoUse),'o','LineWidth',1.5,'Color',ProxColor)
    hold on
    end
    for ee=1:length(AvgNBurstsAllAP(3).EmbryoBursts)
    plot(3,AvgNBurstsAllAP(3).EmbryoBursts(ee).AvgNBursts(APtoUse),'o','LineWidth',1.5,'Color',BothSepColor)
    hold on
    end
    for ee=1:length(AvgNBurstsAllAP(4).EmbryoBursts)
    plot(4,AvgNBurstsAllAP(4).EmbryoBursts(ee).AvgNBursts(APtoUse),'o','LineWidth',1.5,'Color',DoubDistColor)
    hold on
    end
    for ee=1:length(AvgNBurstsAllAP(5).EmbryoBursts)
    plot(5,AvgNBurstsAllAP(5).EmbryoBursts(ee).AvgNBursts(APtoUse),'o','LineWidth',1.5,'Color',DoubProxColor)
    hold on
    end
    for ee=1:length(AvgNBurstsAllAP(6).EmbryoBursts)
    plot(6,AvgNBurstsAllAP(6).EmbryoBursts(ee).AvgNBursts(APtoUse),'o','LineWidth',1.5,'Color',BothColor)
    hold on
    end
    boxplot(APBoxing,'Colors','k');
xlim([0 7]);
xlabel('Construct');
xticks([1:6]);
xticklabels({'Dist', 'Prox', 'Both Sep', '2x Dist', '2x Prox', 'Both'});
ylabel('Number of bursts nc14');
title(['Mean burst frequency',' ' ,num2str(EgglengthUse),'% egg length'])


%ANOVA at specific AP position of constructs against one another
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for bb=1:length([AvgNBurstsAllAP(cc).AllNBursts])
    NBurstComp(bb,cc)=AvgNBurstsAllAP(cc).AllNBursts(bb,APtoUse);
    end
    NBurstComp(NBurstComp==0)=nan;
end
[p,tbl,stats]=anova1(NBurstComp);
xlabel('Construct')
xticks([1:6]);
xticklabels({'Dist','Prox','BothSep','2xDist','2xProx','Both'});
ylabel('Number of bursts');
title(['Mean burst frequency AP bin',' ',num2str(APtoUse)]);



%plot constructs against one another as bar graphs at a single AP position
figure
for ii=1:length(AvgNBurstsAllAP)
    AvgNBurstatAP(ii)=AvgNBurstsAllAP(ii).AvgNBursts(APtoUse);

    errorbar(ii,AvgNBurstatAP(ii), AvgNBurstsAllAP(ii).NBurstSE(APtoUse),'o');
    hold on
end
bar(AvgNBurstatAP,'EdgeColor','k','LineWidth',1.5)
xticks([1:6]);
xticklabels({'Dist','Prox','BothSep','2xDist','2xProx','Both'});
xlim([0 7])
xlabel('Construct')
ylabel('Number of bursts nc14');
ylim([0 5]);
title(['Mean burst frequency',' ',num2str(EgglengthUse),'% egg length']);
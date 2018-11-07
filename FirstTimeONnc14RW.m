ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');

%Count for each construct
AvgFirstOn=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgFirstOnCon=[];
    firsttime=1;
    Con1OnAllAP=[];
    Con1OnSE=[];
    Con1OnSD=[];
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
        OneOnAllAP=[];
        for aa=1:length(APbinID)
            OneOnAP=[];
            FirstOnAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(FirstOnAP)
                OneOnAP=[OneOnAP; nan];
            else
            for bb=1:length(FirstOnAP)
                if ~isempty(BurstProperties(FirstOnAP(bb)).BurstAmplitude)
                OneOnAP=[OneOnAP;[BurstProperties(FirstOnAP(bb)).FirstTimeOn]'];  %put all durations at a given AP value in a column going down
                else
                    OneOnAP=[OneOnAP;nan];
                end
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(OneOnAP)
                OneOnAllAP(bb,aa,ee)=OneOnAP(bb);
            end
            OneOnAllAP(OneOnAllAP==0)=nan;
            
            FirstOnSD(ee,aa,cc)=nanstd(OneOnAllAP(:,aa,ee));
        FirstOnSE(ee,aa,cc)=FirstOnSD(ee,aa,cc)/sqrt(length(FirstOnAP));             
            clear FirstOnAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(OneOnAllAP,3)
            Con1OnAllAP=[Con1OnAllAP; OneOnAllAP(:,:,bb)];
        end
        for bb=1:size(FirstOnSD,3)
            Con1OnSD=[Con1OnSD; FirstOnSD(:,:,bb)];
        end
        for bb=1:size(FirstOnSE,3)
            Con1OnSE=[Con1OnSE;FirstOnSE(:,:,bb)];
        end
        AvgFirstOnAllAP(cc).Embryo1On(ee).Mean1On=nanmean(OneOnAllAP(:,:,ee));
        AvgFirstOnAllAP(cc).Embryo1On(ee).SD=FirstOnSD(ee,:,cc);
        AvgFirstOnAllAP(cc).Embryo1On(ee).SE=FirstOnSE(ee,:,cc);
    end
        AvgFirstOnAllAP(cc).Avg1On=nanmean(Con1OnAllAP);  %Avg duration of all embryos of a construct
        AvgFirstOnAllAP(cc).FirstSD=nanmean(Con1OnSD);
        AvgFirstOnAllAP(cc).FirstSE=nanmean(Con1OnSE);
        AvgFirstOnAllAP(cc).All1Ons=[Con1OnAllAP];

end

%Plot mean of each embryo for each construct
DistalColor=[8 180 238] ./ 255;
DoubDistColor=[1 17 181] ./ 255;
ProxColor=[251 220 50] ./ 255;
DoubProxColor=[251 190 100] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothColor=[12 195 82] ./ 255;

APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for ee=1:length(AvgFirstOnAllAP(cc).Embryo1On)
        APBoxing(ee,cc)=AvgFirstOnAllAP(cc).Embryo1On(ee).Mean1On(APtoUse);
    end
end
APBoxing(APBoxing==0)=nan;
figure 

    for ee=1:length(AvgFirstOnAllAP(1).Embryo1On)
    plot(1,AvgFirstOnAllAP(1).Embryo1On(ee).Mean1On(APtoUse),'o','LineWidth',1.5,'Color',DistalColor)
    hold on
    end
    for ee=1:length(AvgFirstOnAllAP(2).Embryo1On)
    plot(2,AvgFirstOnAllAP(2).Embryo1On(ee).Mean1On(APtoUse),'o','LineWidth',1.5,'Color',ProxColor)
        
    hold on
    end
    
    for ee=1:length(AvgFirstOnAllAP(3).Embryo1On)
    plot(3,AvgFirstOnAllAP(3).Embryo1On(ee).Mean1On(APtoUse),'o','LineWidth',1.5,'Color',BothSepColor)
    hold on
    end
    for ee=1:length(AvgFirstOnAllAP(4).Embryo1On)
    plot(4,AvgFirstOnAllAP(4).Embryo1On(ee).Mean1On(APtoUse),'o','LineWidth',1.5,'Color',DoubDistColor)
    hold on
    end
    for ee=1:length(AvgFirstOnAllAP(5).Embryo1On)
    plot(5,AvgFirstOnAllAP(5).Embryo1On(ee).Mean1On(APtoUse),'o','LineWidth',1.5,'Color',DoubProxColor)
    hold on
    end
    for ee=1:length(AvgFirstOnAllAP(6).Embryo1On)
    plot(6,AvgFirstOnAllAP(6).Embryo1On(ee).Mean1On(APtoUse),'o','LineWidth',1.5,'Color',BothColor);
    %errorbar(6,AvgAmpAllAP(6).EmbryoAmp(ee).MeanProd(APtoUse),AvgAmpAllAP(6).EmbryoAmp(ee).SE(APtoUse),'o','LineWidth',1.5,'Color',BothColor)
    hold on
    end
    boxplot(APBoxing,'Colors','k')
xlim([0 7]);
xlabel('Construct');
xticks([1:6]);
xticklabels({'Dist', 'Prox', 'Both Sep', '2x Dist', '2x Prox','Both'});
ylabel('Time into nc14');
title(['Mean first spot',' ' ,num2str(EgglengthUse),'% egg length'])

%ANOVA at specific AP position of constructs against one another
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for bb=1:length([AvgFirstOnAllAP(cc).All1Ons])
    OneOnComp(bb,cc)=AvgFirstOnAllAP(cc).All1Ons(bb,APtoUse);
    end
    OneOnComp(OneOnComp==0)=nan;
end
[p,tbl,stats]=anova1(OneOnComp);
xlabel('Construct')
xticks([1:6]);
xticklabels({'Dist','Prox','BothSep','2xDist','2xProx','Both'});
ylabel('Time into nc14');
title(['Avg time first spot',' ',num2str(APtoUse)]);

%plot constructs against one another as bar graphs at a single AP position
figure
for ii=1:length(AvgFirstOnAllAP)
    Avg1OnatAP(ii)=AvgFirstOnAllAP(ii).Avg1On(APtoUse);

    errorbar(ii,Avg1OnatAP(ii), AvgFirstOnAllAP(ii).FirstSE(APtoUse),'o');
    hold on
end
bar(Avg1OnatAP,'EdgeColor','k','LineWidth',1.5)
xticks([1:6]);
xticklabels({'Dist','Prox','BothSep','2xDist','2xProx','Both'});
xlim([0 7])
xlabel('Construct')
ylabel('Time into nc14');
title(['Mean time of first spot',' ',num2str(EgglengthUse),'% egg length']);
%% calculate average burst dynamics for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');

%Count for each construct
AvgDuration=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgDurationCon=[];
    firsttime=1;
    ConDurAllAP=[];
    ConDurSE=[];
    ConDurSD=[];
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
        DurAllAP=[];
        for aa=1:length(APbinID)
            DurAP=[];
            DurationAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(DurationAP)
                DurAP=[DurAP; nan];
            else
            for bb=1:length(DurationAP)
                if ~isempty(BurstProperties(DurationAP(bb)).Duration)
                DurAP=[DurAP;[BurstProperties(DurationAP(bb)).Duration]'];  %put all durations at a given AP value in a column going down
                else
                    DurAP=[DurAP; 0];
                end
               
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(DurAP)
                DurAllAP(bb,aa,ee)=DurAP(bb);
            end
            DurAllAP(DurAllAP==0)=nan;
            
            DurationSD(ee,aa,cc)=nanstd(DurAllAP(:,aa,ee));
        DurationSE(ee,aa,cc)=DurationSD(ee,aa,cc)/sqrt(length(DurationAP));             
            clear DurationAP
        end
%         AvgDurAllAP(cc).EmbryosDur(ee).MeanDuration=nanmean(DurAllAP(:,:,ee));
%         AvgDurAllAP(cc).EmbryosDur(ee).SD=DurationSD(ee,:,cc);
%         AvgDurAllAP(cc).EmbryosDur(ee).SE=DurationSE(ee,:,cc);
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(DurAllAP,3)
            ConDurAllAP=[ConDurAllAP; DurAllAP(:,:,bb)];
        end
        for bb=1:size(DurationSD,3)
            ConDurSD=[ConDurSD; DurationSD(:,:,bb)];
        end
        for bb=1:size(DurationSE,3)
            ConDurSE=[ConDurSE;DurationSE(:,:,bb)];
        end
        
    end
        AvgDurAllAP(cc).AvgDur=nanmean(ConDurAllAP);  %Avg duration of all embryos of a construct
        AvgDurAllAP(cc).ConSD=nanmean(ConDurSD);
        AvgDurAllAP(cc).ConSE=nanmean(ConDurSE);
        AvgDurAllAP(cc).AllDurs=[ConDurAllAP];

end

%Plot mean of each embryo at specific AP position for each construct
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
% for cc=1:length(ConstructList)
%     for ee=1:length(AvgDurAllAP(cc).EmbryosDur)
%         APBoxing(ee,cc)=AvgDurAllAP(cc).EmbryosDur(ee).MeanDuration(APtoUse);
%     end
% end
% APBoxing(APBoxing==0)=nan;
% figure
% DistalColor=[8 180 238] ./ 255;
% DoubDistColor=[1 17 181] ./ 255;
% ProxColor=[251 230 60] ./ 255;
% DoubProxColor=[251 190 80] ./ 255;
% BothSepColor=[94 250 81] ./ 255;
% BothColor=[12 195 82] ./ 255;
% Colors=[BothColor; DistalColor; ProxColor; BothSepColor; DoubDistColor; DoubProxColor];
% 
%     for ee=1:length(AvgDurAllAP(1).EmbryosDur)
%     errorbar(1,AvgDurAllAP(1).EmbryosDur(ee).MeanDuration(APtoUse),AvgDurAllAP(1).EmbryosDur(ee).SE(APtoUse),'o','LineWidth',1.5,'Color',DistalColor)
%     hold on
%     end
%     for ee=1:length(AvgDurAllAP(2).EmbryosDur)
%     errorbar(2,AvgDurAllAP(2).EmbryosDur(ee).MeanDuration(APtoUse),AvgDurAllAP(2).EmbryosDur(ee).SE(APtoUse),'o','LineWidth',1.5,'Color',ProxColor)
%     hold on
%     end
%     for ee=1:length(AvgDurAllAP(3).EmbryosDur)
%     errorbar(3,AvgDurAllAP(3).EmbryosDur(ee).MeanDuration(APtoUse),AvgDurAllAP(3).EmbryosDur(ee).SE(APtoUse),'o','LineWidth',1.5,'Color',BothSepColor)
%     hold on
%     end
%     for ee=1:length(AvgDurAllAP(4).EmbryosDur)
%     errorbar(4,AvgDurAllAP(4).EmbryosDur(ee).MeanDuration(APtoUse),AvgDurAllAP(4).EmbryosDur(ee).SE(APtoUse),'o','LineWidth',1.5,'Color',DoubDistColor)
%     hold on
%     end
%     for ee=1:length(AvgDurAllAP(5).EmbryosDur)
%     errorbar(5,AvgDurAllAP(5).EmbryosDur(ee).MeanDuration(APtoUse),AvgDurAllAP(5).EmbryosDur(ee).SE(APtoUse),'o','LineWidth',1.5,'Color',DoubProxColor)
%     hold on
%     end
%     for ee=1:length(AvgDurAllAP(6).EmbryosDur)
%     errorbar(6,AvgDurAllAP(6).EmbryosDur(ee).MeanDuration(APtoUse),AvgDurAllAP(6).EmbryosDur(ee).SE(APtoUse),'o','LineWidth',1.5,'Color',BothColor)
%     hold on
%     end
%     boxplot(APBoxing,'Colors','byg');
% xlim([0 7]);
% xlabel('Construct');
% xticks([1:6]);
% xticklabels({'Dist', 'Prox', 'Both Sep', '2x Dist', '2x Prox', 'Both'});
% ylabel('Burst length (min)');
% title(['Mean burst duration',' ' ,num2str(EgglengthUse),'% egg length'])

%ANOVA of each construct looking at variance with AP position
for cc=1:length(ConstructList)
    anova1(AvgDurAllAP(cc).AllDurs);
    xlabel('AP bin')
    ylabel('Avg burst duration');
    title(ConstructList{cc});
end

%ANOVA at specific AP position of constructs against one another

for cc=1:length(ConstructList)
    for bb=1:length([AvgDurAllAP(cc).AllDurs])
    ConComp(bb,cc)=AvgDurAllAP(cc).AllDurs(bb,APtoUse);
    end
    ConComp(ConComp==0)=nan;
end
[p,tbl,stats]=anova1(ConComp);
xlabel('Construct')
xticks([1:6]);
xticklabels({'Dist','Prox','BothSep','2xDist','2xProx','Both'});
ylabel('Burst length (min)');
title(['Avg Burst Duration AP bin',' ',num2str(APtoUse)]);

%plot the average across all embryos at each AP position for each construct
%as own graph
for cc=1:length(ConstructList)
    figure
    plot(1:length(APbinID),AvgDurAllAP(cc).AvgDur);
    hold on 
    errorbar(1:length(APbinID),AvgDurAllAP(cc).AvgDur,AvgDurAllAP(cc).ConSE,'o');
    xlabel('AP bin');
    ylabel('Burst length (min)');
    title(['Avg burst duration', ConstructList{cc}]);
    ylim([0 14]);
    xlim([0 41]);
end

%plot constructs against one another as bar graphs at a single AP position
figure
for ii=1:length(AvgDurAllAP)
    AvgDuratAP(ii)=AvgDurAllAP(ii).AvgDur(APtoUse);

    errorbar(ii,AvgDuratAP(ii), AvgDurAllAP(ii).ConSE(APtoUse),'o');
    hold on
end
bar(AvgDuratAP,'EdgeColor','k','LineWidth',1.5)
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
xlim([0 7])
xlabel('Construct')
ylabel('Burst length (min)');
ylim([0 10]);
title(['Mean burst duration',' ',num2str(EgglengthUse),'% of egg length']);

%Save duration structure info 
save([DropboxFolder filesep 'Constructs' filesep 'BurstDuration'],'AvgDurAllAP')

        
    
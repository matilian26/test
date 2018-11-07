%% calculate average total spot life for each construct 
%load constructs
ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');

%Count for each construct
AvgSpotLife=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgSpotLifeCon=[];
    firsttime=1;
    ConSLifeAllAP=[];
    ConSLifeSE=[];
    ConSLifeSD=[];
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
        SLifeAllAP=[];
        for aa=1:length(APbinID)
            SLifeAP=[];
            SpotLifeAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(SpotLifeAP)
                SLifeAP=[SLifeAP; nan];
            else
            for bb=1:length(SpotLifeAP)
                if ~isempty(BurstProperties(SpotLifeAP(bb)).Duration)
                SLifeAP=[SLifeAP;[BurstProperties(SpotLifeAP(bb)).SpotLife]'];  %put all durations at a given AP value in a column going down
                else
                    SLifeAP=[SLifeAP; 0];
                end
               
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(SLifeAP)
                SLifeAllAP(bb,aa,ee)=SLifeAP(bb);
            end
            SLifeAllAP(SLifeAllAP==0)=nan;
            
            SpotLifeSD(ee,aa,cc)=nanstd(SLifeAllAP(:,aa,ee));
        SpotLifeSE(ee,aa,cc)=SpotLifeSD(ee,aa,cc)/sqrt(length(SpotLifeAP));             
            clear SpotLifeAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(SLifeAllAP,3)
            ConSLifeAllAP=[ConSLifeAllAP; SLifeAllAP(:,:,bb)];
        end
        for bb=1:size(SpotLifeSD,3)
            ConSLifeSD=[ConSLifeSD; SpotLifeSD(:,:,bb)];
        end
        for bb=1:size(SpotLifeSE,3)
            ConSLifeSE=[ConSLifeSE;SpotLifeSE(:,:,bb)];
        end
        
    end
        AvgSLifeAllAP(cc).AvgSpotLife=nanmean(ConSLifeAllAP);  %Avg duration of all embryos of a construct
        AvgSLifeAllAP(cc).ConSD=nanmean(ConSLifeSD);
        AvgSLifeAllAP(cc).ConSE=nanmean(ConSLifeSE);
        AvgSLifeAllAP(cc).AllSpotLives=[ConSLifeAllAP];

end

%ANOVA of each construct looking at variance with AP position
for cc=1:length(ConstructList)
    anova1(AvgSLifeAllAP(cc).AllSpotLives);
    xlabel('AP bin')
    ylabel('Avg burst duration');
    title(ConstructList{cc});
end

%ANOVA at specific AP position of constructs against one another
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for bb=1:length([AvgSLifeAllAP(cc).AllSpotLives])
    ConComp(bb,cc)=AvgSLifeAllAP(cc).AllSpotLives(bb,APtoUse);
    end
    ConComp(ConComp==0)=nan;
end
[p,tbl,stats]=anova1(ConComp);
xlabel('Construct')
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
ylabel('Spot life (min)');
title(['Avg Spot Life AP bin',' ',num2str(APtoUse)]);

%plot the average across all embryos at each AP position for each construct
%as own graph
for cc=1:length(ConstructList)
    figure
    plot(1:length(APbinID),AvgSLifeAllAP(cc).AvgSpotLife);
    hold on 
    errorbar(1:length(APbinID),AvgSLifeAllAP(cc).AvgSpotLife,AvgSLifeAllAP(cc).ConSE,'o');
    xlabel('AP bin');
    ylabel('Life (min)');
    title(['Mean spot life', ConstructList{cc}]);
    ylim([0 30]);
    xlim([0 41]);
end

%plot constructs against one another as bar graphs at a single AP position
figure
for ii=1:length(AvgSLifeAllAP)
    AvgSLifeatAP(ii)=AvgSLifeAllAP(ii).AvgSpotLife(APtoUse);

    errorbar(ii,AvgSLifeatAP(ii), AvgSLifeAllAP(ii).ConSE(APtoUse),'o');
    hold on
end
bar(AvgSLifeatAP,'EdgeColor','k','LineWidth',1.5)
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
xlim([0 7])
xlabel('Construct')
ylabel('Spot life (min)');
ylim([0 30]);
title(['Mean spot life',' ',num2str(EgglengthUse),'% of egg length']);

%Save duration structure info 
%save([DropboxFolder filesep 'Constructs' filesep 'BurstDuration'],'AvgDurAllAP')
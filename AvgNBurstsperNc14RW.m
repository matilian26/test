%% calculate average burst frequency for each construct 
%load constructs
ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
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
        
    end
        AvgNBurstsAllAP(cc).AvgNBursts=nanmean(ConNBurstsAllAP);  %Avg duration of all embryos of a construct
        AvgNBurstsAllAP(cc).NBurstSD=nanmean(ConNBurstsSD);
        AvgNBurstsAllAP(cc).NBurstSE=nanmean(ConNBurstsSE);
        AvgNBurstsAllAP(cc).AllNBursts=[ConNBurstsAllAP];

end
%ANOVA of each construct looking at variance with AP position
for cc=1:length(ConstructList)
    anova1(AvgNBurstsAllAP(cc).AllNBursts);
    xlabel('AP bin')
    ylabel('Mean number of bursts');
    title(ConstructList{cc});
end

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
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
ylabel('Number of bursts');
title(['Mean burst frequency AP bin',' ',num2str(APtoUse)]);

%plot the average across all embryos at each AP position for each construct
%as own graph
for cc=1:length(ConstructList)
    figure
    plot(1:length(APbinID),AvgNBurstsAllAP(cc).AvgNBursts);
    hold on 
    errorbar(1:length(APbinID),AvgNBurstsAllAP(cc).AvgNBursts,AvgNBurstsAllAP(cc).NBurstSE,'o');
    xlabel('AP bin');
    ylabel('Number of bursts nc14');
    title(['Mean burst frequency', ConstructList{cc}]);
    ylim([0 5]);
    xlim([0 41]);
end

%plot constructs against one another as bar graphs at a single AP position
figure
for ii=1:length(AvgNBurstsAllAP)
    AvgNBurstatAP(ii)=AvgNBurstsAllAP(ii).AvgNBursts(APtoUse);

    errorbar(ii,AvgNBurstatAP(ii), AvgNBurstsAllAP(ii).NBurstSE(APtoUse),'o');
    hold on
end
bar(AvgNBurstatAP,'EdgeColor','k','LineWidth',1.5)
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
xlim([0 7])
xlabel('Construct')
ylabel('Number of bursts nc14');
ylim([0 5]);
title(['Mean burst frequency',' ',num2str(EgglengthUse),'% egg length']);

%save info 
save([DropboxFolder filesep 'Constructs' filesep 'BurstsperNC14'],'AvgNBurstsAllAP')

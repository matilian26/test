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
AvgAmplitude=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgAmplitudeCon=[];
    firsttime=1;
    ConAmpAllAP=[];
    ConAmpSE=[];
    ConAmpSD=[];
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
        AmpAllAP=[];
        for aa=1:length(APbinID)
            AmpAP=[];
            AmplitudeAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(AmplitudeAP)
                AmpAP=[AmpAP; nan];
            else
            for bb=1:length(AmplitudeAP)
                if ~isempty(BurstProperties(AmplitudeAP(bb)).BurstAmplitude)
                AmpAP=[AmpAP;[BurstProperties(AmplitudeAP(bb)).BurstAmplitude]'];  %put all durations at a given AP value in a column going down
                else
                    AmpAP=[AmpAP;nan];
                end
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(AmpAP)
                AmpAllAP(bb,aa,ee)=AmpAP(bb);
            end
            AmpAllAP(AmpAllAP==0)=nan;
            
            AmplitudeSD(ee,aa,cc)=nanstd(AmpAllAP(:,aa,ee));
        AmplitudeSE(ee,aa,cc)=AmplitudeSD(ee,aa,cc)/sqrt(length(AmplitudeAP));             
            clear AmplitudeAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(AmpAllAP,3)
            ConAmpAllAP=[ConAmpAllAP; AmpAllAP(:,:,bb)];
        end
        for bb=1:size(AmplitudeSD,3)
            ConAmpSD=[ConAmpSD; AmplitudeSD(:,:,bb)];
        end
        for bb=1:size(AmplitudeSE,3)
            ConAmpSE=[ConAmpSE;AmplitudeSE(:,:,bb)];
        end
        AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp=AmpAllAP(:,:,ee);
        AvgAmpAllAP(cc).EmbryoAmp(ee).SD=AmplitudeSD(ee,:,cc);
        AvgAmpAllAP(cc).EmbryoAmp(ee).SE=AmplitudeSE(ee,:,cc);
    end
        AvgAmpAllAP(cc).AvgAmp=nanmean(ConAmpAllAP);  %Avg duration of all embryos of a construct
        AvgAmpAllAP(cc).AmpSD=nanmean(ConAmpSD);
        AvgAmpAllAP(cc).AmpSE=nanmean(ConAmpSE);
        AvgAmpAllAP(cc).AllAmps=[ConAmpAllAP];

end
%ANOVA of each construct looking at variance with AP position
for cc=1:length(ConstructList)
    anova1(AvgAmpAllAP(cc).AllAmps);
    xlabel('AP bin')
    ylabel('Avg burst amplitude');
    title(ConstructList{cc});
end

%ANOVA at specific AP position of constructs against one another
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for bb=1:length([AvgAmpAllAP(cc).AllAmps])
    AmpComp(bb,cc)=AvgAmpAllAP(cc).AllAmps(bb,APtoUse);
    end
    AmpComp(AmpComp==0)=nan;
end
[p,tbl,stats]=anova1(AmpComp);
xlabel('Construct')
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
ylabel('Fluorescent intensity');
title(['Avg burst amplitude AP bin',' ',num2str(APtoUse)]);

%plot the average across all embryos at each AP position for each construct
%as own graph
for cc=1:length(ConstructList)
    figure
    plot(1:length(APbinID),AvgAmpAllAP(cc).AvgAmp);
    hold on 
    errorbar(1:length(APbinID),AvgAmpAllAP(cc).AvgAmp,AvgAmpAllAP(cc).AmpSE,'o');
    xlabel('AP bin');
    ylabel('Fluorescent intensity');
    title(['Avg burst amplitude', ConstructList{cc}]);
    ylim([0 60000]);
    xlim([0 41]);
end

%plot constructs against one another as bar graphs at a single AP position
figure
for ii=1:length(AvgAmpAllAP)
    AvgAmpatAP(ii)=AvgAmpAllAP(ii).AvgAmp(APtoUse);

    errorbar(ii,AvgAmpatAP(ii), AvgAmpAllAP(ii).AmpSE(APtoUse),'o');
    hold on
end
bar(AvgAmpatAP,'EdgeColor','k','LineWidth',1.5)
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
xlim([0 7])
xlabel('Construct')
ylabel('Fluorescent intensity (AU)');
ylim([0 50000]);
title(['Mean burst amplitude',' ',num2str(EgglengthUse),'% egg length']);


%save Amp info
save([DropboxFolder filesep 'Constructs' filesep 'BurstAmplitude'],'AvgAmpAllAP')

 

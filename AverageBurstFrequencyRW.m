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
AvgFrequency=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgFrequencyCon=[];
    firsttime=1;
    ConFreqAllAP=[];
    ConFreqSE=[];
    ConFreqSD=[];
    for ee=1:3%4%NEmbryos
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end 
        load(filename);
        NumberBursts(ee,cc)=length([BurstProperties.Duration]);
        %seperate out by AP bin
        FreqAllAP=[];
        for aa=1:length(APbinID)
            FreqAP=[];
            FrequencyAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(FrequencyAP)
                FreqAP=[FreqAP; nan];
            else
            for bb=1:length(FrequencyAP)
                FreqAP=[FreqAP;[BurstProperties(FrequencyAP(bb)).Frequency]'];  %put all durations at a given AP value in a column going down
            end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(FreqAP)
                FreqAllAP(bb,aa,ee)=FreqAP(bb);
            end
            FreqAllAP(FreqAllAP==0)=nan;
            
            FrequencySD(ee,aa,cc)=nanstd(FreqAllAP(:,aa,ee));
        FrequencySE(ee,aa,cc)=FrequencySD(ee,aa,cc)/sqrt(length(FrequencyAP));             
            clear FrequencyAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(FreqAllAP,3)
            ConFreqAllAP=[ConFreqAllAP; FreqAllAP(:,:,bb)];
        end
        for bb=1:size(FrequencySD,3)
            ConFreqSD=[ConFreqSD; FrequencySD(:,:,bb)];
        end
        for bb=1:size(FrequencySE,3)
            ConFreqSE=[ConFreqSE;FrequencySE(:,:,bb)];
        end
        
    end
        AvgFreqAllAP(cc).AvgFreq=nanmean(ConFreqAllAP);  %Avg duration of all embryos of a construct
        AvgFreqAllAP(cc).FreqSD=nanmean(ConFreqSD);
        AvgFreqAllAP(cc).FreqSE=nanmean(ConFreqSE);
        AvgFreqAllAP(cc).AllFreqs=[ConFreqAllAP];

end
%ANOVA of each construct looking at variance with AP position
for cc=1:length(ConstructList)
    anova1(AvgFreqAllAP(cc).AllFreqs);
    xlabel('AP bin')
    ylabel('Avg burst frequency');
    title(ConstructList{cc});
end

%ANOVA at specific AP position of constructs against one another
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for bb=1:length([AvgFreqAllAP(cc).AllFreqs])
    FreqComp(bb,cc)=AvgFreqAllAP(cc).AllFreqs(bb,APtoUse);
    end
    FreqComp(FreqComp==0)=nan;
end
[p,tbl,stats]=anova1(FreqComp);
xlabel('Construct')
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
ylabel('Bursts per min');
title(['Avg burst frequency AP bin',' ',num2str(APtoUse)]);

%plot the average across all embryos at each AP position for each construct
%as own graph
for cc=1:length(ConstructList)
    figure
    plot(1:length(APbinID),AvgFreqAllAP(cc).AvgFreq);
    hold on 
    errorbar(1:length(APbinID),AvgFreqAllAP(cc).AvgFreq,AvgFreqAllAP(cc).FreqSE,'o');
    xlabel('AP bin');
    ylabel('Burst per min');
    title(['Avg burst frequency', ConstructList{cc}]);
    ylim([0 1]);
    xlim([0 41]);
end

%plot constructs against one another as bar graphs at a single AP position
figure
for ii=1:length(AvgFreqAllAP)
    AvgFreqatAP(ii)=AvgFreqAllAP(ii).AvgFreq(APtoUse);

    errorbar(ii,AvgFreqatAP(ii), AvgFreqAllAP(ii).FreqSE(APtoUse),'o');
    hold on
end
bar(AvgFreqatAP)
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
xlim([0 7])
xlabel('Construct')
ylabel('Bursts per ');
ylim([0 0.5]);
title(['Avg burst frequency',' ',num2str(EgglengthUse),'% egg length']);


% %% Calculate average burst frequency for each construct 
% %load constructs
% ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
%     %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
%     
% [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
%  Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
%  ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
% %Ask if want whole movie or just nc14 info
% 
% ncUse=input('Want to just look at nc14?','s');
% 
% %Count for each construct
% AvgFrequency=[];
% for cc=1:length(ConstructList)
%     Data= LoadMS2SetsCS(ConstructList{cc});
%     Datalength(cc)=length(Data);
%     NEmbryos = length(Data);
%     Label = ConstructList(cc);
%     AvgFrequencyCon=[];
%     firsttime=1;
%     for ee=1:4%NEmbryos
%         PrefixName=Data(ee).Prefix;
%         if ncUse=='y'
%         filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
%         else
%             filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
%         end
%         load(filename);
%         for bb=1:length(BurstProperties)
%             if isempty(BurstProperties(bb).Duration)
%                 BurstProperties(bb).Frequency=[];
%             end
%         end
%         NumberBursts(ee,cc)=length([BurstProperties.Duration]);
%         
%         AvgFrequencyCon=[AvgFrequencyCon; [BurstProperties.Frequency]']; %need to change BurstProperties.Duration to a row vector of one column to vertcat all together into one long column
%         if firsttime==1
%             MeanFrequency(ee,cc)=nanmean(AvgFrequencyCon);
%         end
%         if ~isempty(AvgFrequencyCon)
%             if ee > 1
%                 MeanFrequency(ee,cc)=nanmean([BurstProperties.Frequency]);
%                 firsttime=0;
%             end
%         end
%         if ee==1
%         FrequencySD(ee,cc)=std(AvgFrequencyCon);
%         FrequencySE(ee,cc)=FrequencySD(ee,cc)/sqrt(length(AvgFrequencyCon));
%         else
%             FrequencySD(ee,cc)=std([BurstProperties.Frequency]);
%             FrequencySE(ee,cc)=FrequencySD(ee,cc)/sqrt(length([BurstProperties.Frequency]));
%         end
%     end
%     if isempty(AvgFrequency)
%     AvgFrequency=[AvgFrequency AvgFrequencyCon];  % making matrix where each column is a construct and rows going down are individual values of burst duration of all embryos for a construct
%     else
%        if length(AvgFrequency) > length(AvgFrequencyCon)
%            AvgFrequencyCon(end:length(AvgFrequency))=nan;   %making so all same length in number of observations for duration- this isn't the best way to do this
%            AvgFrequency=[AvgFrequency AvgFrequencyCon];
%        elseif length(AvgFrequency) < length(AvgFrequencyCon)
%         AvgFrequency=[AvgFrequency AvgFrequencyCon(1:length(AvgFrequency))];
%        end
%     end
% end
% 
% [p,tbl, stats]= anova1(AvgFrequency);
% xlabel('Construct')
% xticks([1:6])
% xticklabels({'Both', 'Dist', 'Prox', 'Sep', '2xDist', '2xProx'}) 
% ylabel('Avg burst frequency')
% 
% figure
% for cc=1:length(ConstructList)
%     for ee=1:4
%         %plot(cc,MeanDuration(ee,cc),'o')
%         errorbar(cc, MeanFrequency(ee,cc), FrequencySE(ee,cc), 'o')
%         hold on 
%     end
%     xlabel('Construct')
%     xticks([1:6]);
%     xticklabels({'Both', 'Dist', 'Prox', 'Sep', '2xDist', '2xProx'}) 
%     ylabel('Avg burst frequency')
%     xlim([0 7]);
%     
% end
% 
% figure
% for cc=1:length(ConstructList)
%     TotalAvg(cc)=mean(MeanFrequency(:,cc));
%     TotalSE(cc)=mean(FrequencySE(:,cc));
%     bar(TotalAvg);
%     hold on 
%     errorbar(cc,TotalAvg(cc),TotalSE(cc),'o','LineWidth',2);
%     title('Avg Frequency');
%     xlabel('Construct');
%     xticks([1:6]);
%     xticklabels({'Both', 'Dist', 'Prox', 'Sep', '2xDist', '2xProx'});
%     ylabel('bursts per minute');
% end
%% calculate average burst frequency for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'Kr2xProxEmpty';'KrDist17C';'KrBoth17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
SlopeUse=input('Want to use slope calculations?','s');
%Count for each construct
AvgFreq=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgFreqCon=[];
    firsttime=1;
    ConFreqAllAP=[];
    ConFreqSE=[];
    ConFreqSD=[];
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        if SlopeUse=='y'
            filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        elseif ncUse=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end 
        load(filename);
%         if SlopeUse=='y'
%         Frequency(ee,cc)=length([BurstProperties.FrequencyAct]);
%         else
          Frequency(ee,cc)=length([BurstProperties.Frequency]);  
        %end
        %seperate out by AP bin
        FreqAllAP=[];
        for aa=1:length(APbinID)
            FreqAP=[];
            FrequencyAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(FrequencyAP)
                FreqAP=[FreqAP; nan];
            else
            for bb=1:length(FrequencyAP)
%                 if SlopeUse=='y'
%                 FreqAP=[FreqAP;[BurstProperties(FrequencyAP(bb)).FrequencyAct]'];  %put all durations at a given AP value in a column going down
%                 else
                   FreqAP=[FreqAP;[BurstProperties(FrequencyAP(bb)).Frequency]']; 
                %end
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(FreqAP)
                FreqAllAP(bb,aa,ee)=FreqAP(bb);
            end
            FreqAllAP(FreqAllAP==0)=nan;
            
            FrequencySD(ee,aa,cc)=nanstd(FreqAllAP(:,aa,ee));
        FrequencySE(ee,aa,cc)=FrequencySD(ee,aa,cc)/sqrt(sum(~isnan(FreqAP)));             
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
        AvgFreqAllAP(cc).EmbryoFreq(ee).AvgFreq=nanmean(FreqAllAP(:,:,ee));
        AvgFreqAllAP(cc).EmbryoFreq(ee).SE=FrequencySE(ee,:,cc);
        AvgFreqAllAP(cc).EmbryoFreq(ee).SD=FrequencySD(ee,:,cc);
    end
        AvgFreqAllAP(cc).AvgFreq=nanmean(ConFreqAllAP);  %Avg duration of all embryos of a construct
        AvgFreqAllAP(cc).FreqSD=nanmean(ConFreqSD);
        AvgFreqAllAP(cc).FreqSE=nanmean(ConFreqSE);
        AvgFreqAllAP(cc).AllFreq=[ConFreqAllAP];
        AvgFreqAllAP(cc).AllFreqSD=nanstd(AvgFreqAllAP(cc).AllFreq);
        %AvgFreqAllAP(cc).AllFreqSE=(([AvgFreqAllAP(cc).AllFreqSD])/(sqrt(length(AvgFreqAllAP(cc).AllFreq))));
        for aa=1:length(APbinID)
        AvgFreqAllAP(cc).AllFreqSE(aa)=([AvgFreqAllAP(cc).AllFreqSD(aa)])/(sqrt(sum(~isnan(AvgFreqAllAP(cc).AllFreq(:,aa)))));
        end
        AvgFreqAllAP(cc).AllFreq95Conf=[AvgFreqAllAP(cc).AllFreqSE] .* 1.95;

end
% Get rid of single data points for an AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(AvgFreqAllAP(cc).AllFreq(:,aa)))==1
            AvgFreqAllAP(cc).AllFreqSD(aa)=nan;
            AvgFreqAllAP(cc).AvgFreq(aa)=nan;
            AvgProdAllAP(cc).AllFreqSE(aa)=nan;
            AvgProdAllAP(cc).AllFreq95Conf(aa)=nan;
        end
    end
end

%% Plot mean of each embryo at specific AP position for each construct
APtoUse=input('Which AP bin to compare constructs?');
Egglength=APbinID.*100;
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for ee=1:length(AvgFreqAllAP(cc).EmbryoFreq)
        APBoxing(ee,cc)=AvgFreqAllAP(cc).EmbryoFreq(ee).AvgFreq(APtoUse);
    end
end
APBoxing(APBoxing==0)=nan;

DistalColor=[1 64 172]./255;
DistalEmptyColor=[8 210 238] ./ 255;
Distal32CColor=[118 180 238] ./ 255;
DoubDistColor=[73 184 253] ./ 255;
ProxColor=[238 123 23]./255;
ProxEmptyColor=[251 250 50] ./255;
Proximal32CColor=[251 150 10] ./ 255;
DoubProxColor=[215 183 58] ./ 255;
DoubProxEmptyColor=[251 220 50] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[52 119 71]./255;
BothEmptyColor=[12 250 100] ./ 255;
Both32CColor=[120 195 82] ./ 255;
DoubProx32CColor=[200 150 100] ./ 255;


Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DistalEmptyColor;
Colors(5).Color=ProxEmptyColor;
Colors(6).Color=DoubDistColor;
Colors(7).Color=DoubProxColor;
Colors(8).Color=BothColor;
Colors(9).Color=BothEmptyColor;
Colors(10).Color=DistalColor
Colors(11).Color=ProxColor;
Colors(12).Color=BothSepColor;
Colors(13).Color=BothColor;
Colors(14).Color=DoubProxColor;
Colors(15).Color=DoubProxColor;
Colors(16).Color=DistalColor;
Colors(17).Color=BothColor;

fontsize=18;
fontname='Helvetica';
FigDirect=[DropboxFolder filesep 'Figures' filesep 'Transcriptional dynamics' filesep 'Frequency'];
%%
    for ee=1:length(AvgFreqAllAP(1).EmbryoFreq)
    plot(1,AvgFreqAllAP(1).EmbryoFreq(ee).AvgFreq(APtoUse),'o','LineWidth',1.5,'Color',DistalColor)
    hold on
    end
    for ee=1:length(AvgFreqAllAP(2).EmbryoFreq)
    plot(2,AvgFreqAllAP(2).EmbryoFreq(ee).AvgFreq(APtoUse),'o','LineWidth',1.5,'Color',ProxColor)
    hold on
    end
    for ee=1:length(AvgFreqAllAP(3).EmbryoFreq)
    plot(3,AvgFreqAllAP(3).EmbryoFreq(ee).AvgFreq(APtoUse),'o','LineWidth',1.5,'Color',BothSepColor)
    hold on
    end
    for ee=1:length(AvgFreqAllAP(4).EmbryoFreq)
    plot(4,AvgFreqAllAP(4).EmbryoFreq(ee).AvgFreq(APtoUse),'o','LineWidth',1.5,'Color',Colors(4).Color)
    hold on
    end
    for ee=1:length(AvgFreqAllAP(5).EmbryoFreq)
    plot(5,AvgFreqAllAP(5).EmbryoFreq(ee).AvgFreq(APtoUse),'o','LineWidth',1.5,'Color',Colors(5).Color)
    hold on
    end
    for ee=1:length(AvgFreqAllAP(6).EmbryoFreq)
    plot(6,AvgFreqAllAP(6).EmbryoFreq(ee).AvgFreq(APtoUse),'o','LineWidth',1.5,'Color',Colors(6).Color)
    hold on
    end
    for ee=1:length(AvgFreqAllAP(7).EmbryoFreq)
    plot(7,AvgFreqAllAP(7).EmbryoFreq(ee).AvgFreq(APtoUse),'o','LineWidth',1.5,'Color',Colors(7).Color)
    hold on
    end
     for ee=1:length(AvgFreqAllAP(8).EmbryoFreq)
    plot(8,AvgFreqAllAP(8).EmbryoFreq(ee).AvgFreq(APtoUse),'o','LineWidth',1.5,'Color',Colors(8).Color)
    hold on
     end
     for ee=1:length(AvgFreqAllAP(9).EmbryoFreq)
    plot(9,AvgFreqAllAP(9).EmbryoFreq(ee).AvgFreq(APtoUse),'o','LineWidth',1.5,'Color',Colors(9).Color)
    hold on
    end
    boxplot(APBoxing,'Colors','k');
xlim([0 12]);
xlabel('Construct');
xticks([1:11]);
xticklabels({'Dist', 'Prox', 'Both Sep','1x Dist','1x Prox', '2x Dist', '2x Prox', 'Both','1x Both','Dist32','BothSep32'});
ylabel('Bursts per minute');
title(['Mean burst frequency',' ' ,num2str(EgglengthUse),'% egg length'])


%ANOVA at specific AP position of constructs against one another
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for bb=1:length([AvgFreqAllAP(cc).AllFreq])
    FreqComp(bb,cc)=AvgFreqAllAP(cc).AllFreq(bb,APtoUse);
    end
    FreqComp(FreqComp==0)=nan;
end
[p,tbl,stats]=anova1(FreqComp);
xlabel('Construct')
xticks([1:9]);
xticklabels({'Dist','Prox','BothSep','1x Dist','1x Prox','2xDist','2xProx','Both','1x Both'});
ylabel('Frequency');
title(['Mean burst frequency AP bin',' ',num2str(APtoUse)]);

%Plot Construct means vs AP position
figure 
for cc=1:length(ConstructList);
    plot(Egglength,AvgFreqAllAP(cc).AvgFreq,'Color',Colors(cc).Color,'LineWidth',1.5);
    hold on
end
legend('Dist','Prox','Both Sep','1x Dist','1x Prox', '2x Dist','2x Prox', 'Both','1x Both');
xlabel('% Egg length');
xlim([0 100]);
ylabel('Bursts per min');
title('Mean burst frequency');
%%
%Singles vs Both Sep
figure 
errorbar(Egglength, AvgFreqAllAP(1).AvgFreq, AvgFreqAllAP(1).AllFreq95Conf,'Color', Colors(1).Color,'LineWidth',1.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(2).AvgFreq, AvgFreqAllAP(2).AllFreq95Conf,'Color', Colors(2).Color,'LineWidth',1.5)
errorbar(Egglength, AvgFreqAllAP(3).AvgFreq, AvgFreqAllAP(3).AllFreq95Conf,'Color', Colors(3).Color,'LineWidth',1.5)
legend('Dist', 'Prox', 'Both Sep')
xlabel('% Egg length');
xlim([0 100]);
ylabel('Bursts per minute');
title('Mean burst frequency');

figure 
errorbar(Egglength, AvgFreqAllAP(1).AvgFreq, AvgFreqAllAP(1).AllFreq95Conf,'Color', Colors(1).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(2).AvgFreq, AvgFreqAllAP(2).AllFreq95Conf,'Color', Colors(2).Color,'LineWidth',2.5)
%legend('Dist', 'Prox')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
%title('Mean burst frequency');

%% Hemizygotes
figure 
errorbar(Egglength, AvgFreqAllAP(1).AvgFreq, AvgFreqAllAP(1).AllFreq95Conf,'Color', Colors(1).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(4).AvgFreq, AvgFreqAllAP(4).AllFreq95Conf,'Color', Colors(1).Color,'LineWidth',2.5,'LineStyle',':')
%legend('Distal', '1x Distal')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print( [FigDirect filesep 'HemiDistFreq'],'-dsvg');

figure 
errorbar(Egglength, AvgFreqAllAP(2).AvgFreq, AvgFreqAllAP(2).AllFreq95Conf,'Color', Colors(2).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(5).AvgFreq, AvgFreqAllAP(5).AllFreq95Conf,'Color', Colors(2).Color,'LineWidth',2.5,'LineStyle',':')
%legend('Proximal', '1x Proximal')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print( [FigDirect filesep 'HemiProxFreq'],'-dsvg');

figure 
errorbar(Egglength, AvgFreqAllAP(8).AvgFreq, AvgFreqAllAP(8).AllFreq95Conf,'Color', Colors(8).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(9).AvgFreq, AvgFreqAllAP(9).AllFreq95Conf,'Color', Colors(8).Color,'LineWidth',2.5,'LineStyle',':')
%legend('Both', '1x Both')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print( [FigDirect filesep 'HemiBothFreq'],'-dsvg');

figure 
errorbar(Egglength, AvgFreqAllAP(7).AvgFreq, AvgFreqAllAP(7).AllFreq95Conf,'Color', Colors(7).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(15).AvgFreq, AvgFreqAllAP(15).AllFreq95Conf,'Color', Colors(15).Color,'LineWidth',2.5,'LineStyle',':')
%legend('Both', '1x Both')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print( [FigDirect filesep 'Hemi2xProxFreq'],'-dsvg');
%%
figure 
errorbar(Egglength, AvgFreqAllAP(1).AvgFreq, AvgFreqAllAP(1).AllFreq95Conf,'Color', Colors(1).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(2).AvgFreq, AvgFreqAllAP(2).AllFreq95Conf,'Color', Colors(2).Color,'LineWidth',2.5)
errorbar(Egglength, AvgFreqAllAP(8).AvgFreq, AvgFreqAllAP(8).AllFreq95Conf,'Color', Colors(8).Color,'LineWidth',2.5)
%legend('Distal', 'Proximal','Both')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print( [FigDirect filesep 'SinglesvBothFreq'],'-dsvg');
%%
%Duplicates vs SE
figure 
errorbar(Egglength,AvgFreqAllAP(6).AvgFreq,AvgFreqAllAP(6).AllFreqSE,'Color',Colors(6).Color,'LineWidth',2.5);
hold on
errorbar(Egglength,AvgFreqAllAP(7).AvgFreq,AvgFreqAllAP(7).AllFreqSE,'Color',Colors(7).Color,'LineWidth',2.5);
errorbar(Egglength,AvgFreqAllAP(8).AvgFreq,AvgFreqAllAP(8).AllFreqSE,'Color',Colors(8).Color,'LineWidth',2.5);
%legend('2x Prox', 'Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
xlim([0 100]);
ylabel('Bursts per min');
%title('Mean burst frequency');

%Doubles vs singles
figure 
errorbar(Egglength,AvgFreqAllAP(1).AvgFreq,AvgFreqAllAP(1).AllFreqSE,'Color',Colors(1).Color,'LineWidth',1.5);
hold on
errorbar(Egglength,AvgFreqAllAP(6).AvgFreq,AvgFreqAllAP(6).AllFreqSE,'Color',Colors(6).Color,'LineWidth',1.5);
legend('Dist','2x Dist');
xlabel('% Egg length');
xlim([0 100]);
ylabel('Bursts per min');
title('Mean burst frequency');

figure 
errorbar(Egglength,AvgFreqAllAP(2).AvgFreq,AvgFreqAllAP(2).AllFreqSE,'Color',Colors(2).Color,'LineWidth',1.5);
hold on
errorbar(Egglength,AvgFreqAllAP(7).AvgFreq,AvgFreqAllAP(7).AllFreqSE,'Color',Colors(7).Color,'LineWidth',1.5);
legend('Prox','2x Prox');
xlabel('% Egg length');
xlim([0 100]);
ylabel('Bursts per min');
title('Mean burst frequency');

figure 
errorbar(Egglength,AvgFreqAllAP(3).AvgFreq,AvgFreqAllAP(3).AllFreqSE,'Color',Colors(3).Color,'LineWidth',1.5);
hold on
errorbar(Egglength,AvgFreqAllAP(8).AvgFreq,AvgFreqAllAP(8).AllFreqSE,'Color',Colors(8).Color,'LineWidth',1.5);
legend('Both Sep','Both');
xlabel('% Egg length');
xlim([0 100]);
ylabel('Bursts per min');
title('Mean burst frequency');
%%
figure 
errorbar(Egglength, AvgFreqAllAP(4).AvgFreq, AvgFreqAllAP(4).AllFreq95Conf,'Color', Colors(4).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(8).AvgFreq, AvgFreqAllAP(8).AllFreq95Conf,'Color', Colors(8).Color,'LineWidth',2.5)
legend('1x Distal', 'Both')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
xlim([0 100]);
ylabel('Bursts per minute');
%title('Mean burst frequency');

figure 
errorbar(Egglength, AvgFreqAllAP(5).AvgFreq, AvgFreqAllAP(5).AllFreq95Conf,'Color', Colors(5).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(8).AvgFreq, AvgFreqAllAP(8).AllFreq95Conf,'Color', Colors(8).Color,'LineWidth',2.5)
legend('1x Proximal', 'Both')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
xlim([0 100]);
ylabel('Bursts per minute');
%% Temperature comparisions 
figure 
errorbar(Egglength, AvgFreqAllAP(1).AvgFreq, AvgFreqAllAP(1).AllFreq95Conf,'Color', Colors(1).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(10).AvgFreq, AvgFreqAllAP(10).AllFreq95Conf,'Color', Colors(10).Color,'LineWidth',2.5,'LineStyle','-.')
%legend('Distal', 'Distal 32C')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print( [FigDirect filesep 'DistTCompFreq'],'-dsvg');

figure 
errorbar(Egglength, AvgFreqAllAP(2).AvgFreq, AvgFreqAllAP(2).AllFreq95Conf,'Color', Colors(2).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(11).AvgFreq, AvgFreqAllAP(11).AllFreq95Conf,'Color', Colors(11).Color,'LineWidth',2.5,'LineStyle','-.')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print( [FigDirect filesep 'ProxTCompFreq'],'-dsvg');

figure 
errorbar(Egglength, AvgFreqAllAP(3).AvgFreq, AvgFreqAllAP(3).AllFreq95Conf,'Color', Colors(3).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(12).AvgFreq, AvgFreqAllAP(12).AllFreq95Conf,'Color', Colors(12).Color,'LineWidth',2.5,'LineStyle','-.')
%legend('Both Sep', 'Both Sep 32C')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');

figure 
errorbar(Egglength, AvgFreqAllAP(8).AvgFreq, AvgFreqAllAP(8).AllFreq95Conf,'Color', Colors(8).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(13).AvgFreq, AvgFreqAllAP(13).AllFreq95Conf,'Color', Colors(13).Color,'LineWidth',2.5,'LineStyle','-.')
%legend('Both', 'Both 32C')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print( [FigDirect filesep 'BothTCompFreq'],'-dsvg');

figure 
errorbar(Egglength, AvgFreqAllAP(7).AvgFreq, AvgFreqAllAP(7).AllFreq95Conf,'Color', Colors(7).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(14).AvgFreq, AvgFreqAllAP(14).AllFreq95Conf,'Color', Colors(14).Color,'LineWidth',2.5,'LineStyle','-.')
%legend('2x Proximal', '2x Proximal 32C','Location','best')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print( [FigDirect filesep '2xProxTCompFreq'],'-dsvg');
%% 17C Comparisions
figure 
errorbar(Egglength, AvgFreqAllAP(1).AvgFreq, AvgFreqAllAP(1).AllFreq95Conf,'Color', Colors(1).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(16).AvgFreq, AvgFreqAllAP(16).AllFreq95Conf,'Color', Colors(16).Color,'LineWidth',2.5,'LineStyle','-.')
%legend('Distal', 'Distal 32C')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
print( [FigDirect filesep 'Dist17CCompFreq'],'-dsvg');

figure 
errorbar(Egglength, AvgFreqAllAP(8).AvgFreq, AvgFreqAllAP(8).AllFreq95Conf,'Color', Colors(8).Color,'LineWidth',2.5)
hold on 
errorbar(Egglength, AvgFreqAllAP(17).AvgFreq, AvgFreqAllAP(17).AllFreq95Conf,'Color', Colors(17).Color,'LineWidth',2.5,'LineStyle','-.')
%legend('Distal', 'Distal 32C')
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('bursts per minute');
%%

%save frequency info
%save([DropboxFolder filesep 'Constructs' filesep 'BurstFrequency'],'AvgFreqAllAP');
save([DropboxFolder filesep 'Constructs' filesep 'BurstFrequencySlope'],'AvgFreqAllAP');

%% 2way anova of frequency vs AP bin vs construct 
%% 2way anova of Duration vs AP position vs genotype 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
%ncUse=input('Want to only use nc14?','s');

%Count for each construct

FrequencyManovaVect=[];
APManovaVect=[];
ConManovaVect=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end 
        load(filename);
        NumberBursts(ee,cc)=length([BurstProperties.Duration]);
        for nn=1:length(BurstProperties)
            FrequencyManovaVect=[FrequencyManovaVect, nanmean([BurstProperties(nn).Frequency])];
            APManovaVect=[APManovaVect, BurstProperties(nn).APBin];
            ConManovaVect=[ConManovaVect, cc];
        end
    end
end

[p,tbl,stats]=anovan(FrequencyManovaVect,{APManovaVect, ConManovaVect},'sstype',1,'varnames',{'AP bin','Construct'})
%multcompare(stats,'Dimension',2)
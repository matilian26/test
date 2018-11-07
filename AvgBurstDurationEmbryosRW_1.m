%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'Kr2xProxEmpty';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32c';'KrDist17C';'KrBoth17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
SlopeUse=input('Want to use slope calculations?','s');

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
    for ee=1:NEmbryos %round(NEmbryos/2)
        PrefixName=Data(ee).Prefix;
        if SlopeUse=='y'
            filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        elseif ncUse=='y'
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
                    DurAP=[DurAP; nan];
                end
               
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(DurAP)
                DurAllAP(bb,aa,ee)=DurAP(bb);
            end
            DurAllAP(DurAllAP==0)=nan;
            
            DurationSD(ee,aa,cc)=nanstd(DurAllAP(:,aa,ee));
        %DurationSE(ee,aa,cc)=DurationSD(ee,aa,cc)/sqrt(length(DurationAP));
        %7.16.18
        DurationSE(ee,aa,cc)=DurationSD(ee,aa,cc)/sqrt(sum(~isnan(DurAP)));
            clear DurationAP
        end
        AvgDurAllAP(cc).EmbryosDur(ee).MeanDuration=nanmean(DurAllAP(:,:,ee));
        AvgDurAllAP(cc).EmbryosDur(ee).SD=DurationSD(ee,:,cc);
        AvgDurAllAP(cc).EmbryosDur(ee).SE=DurationSE(ee,:,cc);
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
        AvgDurAllAP(cc).AllSD=nanstd(ConDurAllAP);
        %AvgDurAllAP(cc).AllSE=([AvgDurAllAP(cc).AllSD])/(sqrt(length(AvgDurAllAP(cc).AllDurs)));
        for aa=1:length(APbinID)
            AvgDurAllAP(cc).AllSE(aa)=([AvgDurAllAP(cc).AllSD(aa)])/(sqrt(sum(~isnan(AvgDurAllAP(cc).AllDurs(:,aa)))));
        end
        AvgDurAllAP(cc).All95Conf=([AvgDurAllAP(cc).AllSE]).*1.95;

end
%Get rid of single data points of a whole AP bin
for cc=1:length(ConstructList)
    AvgDurAllAP(cc).Construct=ConstructList{cc}
    for aa=1:length(APbinID)
        if sum(~isnan(AvgDurAllAP(cc).AllDurs(:,aa)))==1
            AvgDurAllAP(cc).AllSD(aa)=nan;
            AvgDurAllAP(cc).AvgDur(aa)=nan;
            AvgDurAllAP(cc).AllSE(aa)=nan;
            AvgDurAllAP(cc).All95Conf(aa)=nan;
        end
    end
end

%% Plot mean of each embryo at specific AP position for each construct
APtoUse=input('Which AP bin to compare constructs?');
Egglength=APbinID.*100;
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for ee=1:length(AvgDurAllAP(cc).EmbryosDur)
        APBoxing(ee,cc)=AvgDurAllAP(cc).EmbryosDur(ee).MeanDuration(APtoUse);
    end
end
APBoxing(APBoxing==0)=nan;
figure
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
Colors(8).Color=DoubProxEmptyColor;
Colors(9).Color=BothColor;
Colors(10).Color=BothEmptyColor;
Colors(11).Color=DistalColor;
Colors(12).Color=ProxColor;
Colors(13).Color=BothSepColor;
Colors(14).Color=BothColor;
Colors(15).Color=DoubProxColor;
Colors(16).Color=DistalColor;
Colors(17).Color=BothColor;

fontsize=18;
fontname='Helvetica';
FigDirect=[DropboxFolder filesep 'Figures'];
%%
    for ee=1:length(AvgDurAllAP(1).EmbryosDur)
    plot(1,AvgDurAllAP(1).EmbryosDur(ee).MeanDuration(APtoUse),'o','LineWidth',1.5,'Color',DistalColor)
    hold on
    end
    for ee=1:length(AvgDurAllAP(2).EmbryosDur)
    plot(2,AvgDurAllAP(2).EmbryosDur(ee).MeanDuration(APtoUse),'o','LineWidth',1.5,'Color',ProxColor)
    hold on
    end
    for ee=1:length(AvgDurAllAP(3).EmbryosDur)
    plot(3,AvgDurAllAP(3).EmbryosDur(ee).MeanDuration(APtoUse),'o','LineWidth',1.5,'Color',BothSepColor)
    hold on
    end
    for ee=1:length(AvgDurAllAP(4).EmbryosDur)
    plot(4,AvgDurAllAP(4).EmbryosDur(ee).MeanDuration(APtoUse),'o','LineWidth',1.5,'Color',Colors(4).Color)
    hold on
    end
    for ee=1:length(AvgDurAllAP(5).EmbryosDur)
    plot(5,AvgDurAllAP(5).EmbryosDur(ee).MeanDuration(APtoUse),'o','LineWidth',1.5,'Color',Colors(5).Color)
    hold on
    end
    for ee=1:length(AvgDurAllAP(6).EmbryosDur)
    plot(6,AvgDurAllAP(6).EmbryosDur(ee).MeanDuration(APtoUse),'o','LineWidth',1.5,'Color',Colors(6).Color)
    hold on
    end
    for ee=1:length(AvgDurAllAP(7).EmbryosDur)
    plot(7,AvgDurAllAP(7).EmbryosDur(ee).MeanDuration(APtoUse),'o','LineWidth',1.5,'Color',BothColor)
    hold on
    end
     for ee=1:length(AvgDurAllAP(8).EmbryosDur)
    plot(8,AvgDurAllAP(8).EmbryosDur(ee).MeanDuration(APtoUse),'o','LineWidth',1.5,'Color',Colors(8).Color)
    hold on
     end
    for ee=1:length(AvgDurAllAP(9).EmbryosDur)
    plot(9,AvgDurAllAP(9).EmbryosDur(ee).MeanDuration(APtoUse),'o','LineWidth',1.5,'Color',Colors(9).Color)
    hold on
    end
    boxplot(APBoxing,'Colors','k');
xlim([0 10]);
xlabel('Construct');
xticks([1:9]);
xticklabels({'Dist', 'Prox', 'Both Sep', '1x Dist','1x Prox', '2x Dist', '2x Prox', 'Both','1x Both'});
ylabel('Burst length (min)');
title(['Mean burst duration',' ' ,num2str(EgglengthUse),'% egg length'])

%ANOVA at specific AP position of constructs against one another

for cc=1:length(ConstructList)
    for bb=1:length([AvgDurAllAP(cc).AllDurs])
    ConComp(bb,cc)=AvgDurAllAP(cc).AllDurs(bb,APtoUse);
    end
    ConComp(ConComp==0)=nan;
end
[p,tbl,stats]=anova1(ConComp);
xlabel('Construct')
xticks([1:11]);
xticklabels({'Dist','Prox','BothSep','1xDist','1xProx','2xDist','2xProx','Both','1xBoth','Dist32','BothSep32'});
ylabel('Burst length (min)');
title(['Avg Burst Duration AP bin',' ',num2str(APtoUse)]);

% plot construct means vs AP
figure 
for cc=1:length(ConstructList)
    plot(Egglength,AvgDurAllAP(cc).AvgDur,'Color',Colors(cc).Color,'LineWidth',1.5);
    hold on
end
legend('Dist','Prox','Both Sep','1x Dist','1x Prox', '2x Dist', '2x Prox', 'Both','1x Both');
title('Mean burst duration');
xlabel('% Egg length');
ylabel('Time (min');
xlim([0 100]);
%%
%doubles vs SE
figure
%errorbar(Egglength,AvgDurAllAP(6).AvgDur,AvgDurAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength,AvgDurAllAP(7).AvgDur,AvgDurAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
errorbar(Egglength,AvgDurAllAP(9).AvgDur,AvgDurAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
legend('2x Prox', 'Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
title('Mean burst duration');
xlabel('% Egg length');
ylabel('Time (min');
xlim([0 100]);

%singles vs each other or Both sep
figure
errorbar(Egglength,AvgDurAllAP(1).AvgDur,AvgDurAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(2).AvgDur,AvgDurAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
%legend('Distal','Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
%title('Mean burst duration');
xlabel('% egg length');
ylabel('burst duration (min)');
%ylabel('Time (min)');
xlim([0 100]);
print( [FigDirect filesep 'SinglesDuration'],'-dsvg');

%%
figure
errorbar(Egglength,AvgDurAllAP(1).AvgDur,AvgDurAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(4).AvgDur,AvgDurAllAP(4).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5,'LineStyle',':');
legend('Distal','1x Distal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
%title('Mean burst duration');
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
print( [FigDirect filesep 'DistHemiCompDuration'],'-dsvg');

figure
errorbar(Egglength,AvgDurAllAP(2).AvgDur,AvgDurAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(5).AvgDur,AvgDurAllAP(5).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5,'LineStyle',':');
%legend('Proximal','1x Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
%title('Mean burst duration');
xlabel('% egg length');
ylabel('burst duration (min)');
%ylabel('Time (min)');
xlim([0 100]);
print( [FigDirect filesep 'ProxHemiCompDuration'],'-dsvg');


figure
errorbar(Egglength,AvgDurAllAP(9).AvgDur,AvgDurAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(10).AvgDur,AvgDurAllAP(10).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle',':');
%legend('Both','1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
%title('Mean burst duration');
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
print( [FigDirect filesep 'BothHemiCompDuration'],'-dsvg');


figure
errorbar(Egglength,AvgDurAllAP(7).AvgDur,AvgDurAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(8).AvgDur,AvgDurAllAP(8).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5,'LineStyle',':');
%legend('Both','1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
%title('Mean burst duration');
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
print( [FigDirect filesep '2xProxHemiCompDuration'],'-dsvg');

%%
figure
errorbar(Egglength,AvgDurAllAP(1).AvgDur,AvgDurAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength,AvgDurAllAP(2).AvgDur,AvgDurAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',1.5);
errorbar(Egglength,AvgDurAllAP(3).AvgDur,AvgDurAllAP(3).All95Conf,'Color',Colors(3).Color,'LineWidth',1.5);
legend('Distal','Proximal','Both Sep');
title('Mean burst duration');
xlabel('% Egg length');
ylabel('Time (min');
xlim([0 100]);

figure
errorbar(Egglength,AvgDurAllAP(1).AvgDur,AvgDurAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(2).AvgDur,AvgDurAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
errorbar(Egglength,AvgDurAllAP(9).AvgDur,AvgDurAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
%legend('Distal','Proximal','Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
%title('Mean burst duration');
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
print( [FigDirect filesep 'SinglesvBothDuration'],'-dsvg');

%%
%doubles vs singles
figure
errorbar(Egglength,AvgDurAllAP(1).AvgDur,AvgDurAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength,AvgDurAllAP(6).AvgDur,AvgDurAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',1.5);
legend('Dist','2x Dist');
title('Mean burst duration');
xlabel('% Egg length');
ylabel('Time (min');
xlim([0 100]);

figure
errorbar(Egglength,AvgDurAllAP(2).AvgDur,AvgDurAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength,AvgDurAllAP(7).AvgDur,AvgDurAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',1.5);
legend('Prox','2x Prox');
title('Mean burst duration');
xlabel('% Egg length');
ylabel('Time (min');
xlim([0 100]);

figure
errorbar(Egglength,AvgDurAllAP(3).AvgDur,AvgDurAllAP(3).All95Conf,'Color',Colors(3).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength,AvgDurAllAP(9).AvgDur,AvgDurAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',1.5);
legend('Both Sep','Both');
title('Mean burst duration');
xlabel('% Egg length');
ylabel('Time (min');
xlim([0 100]);
%%
figure
errorbar(Egglength,AvgDurAllAP(4).AvgDur,AvgDurAllAP(4).All95Conf,'Color',Colors(4).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength,AvgDurAllAP(9).AvgDur,AvgDurAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',1.5);
legend('1x Distal','Both');
%title('Mean burst duration');
xlabel('% Egg length');
ylabel('Burst duration (min');
xlim([0 100]);

figure
errorbar(Egglength,AvgDurAllAP(5).AvgDur,AvgDurAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength,AvgDurAllAP(9).AvgDur,AvgDurAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',1.5);
legend('1x Proximal','Both');
%title('Mean burst duration');
xlabel('% Egg length');
ylabel('Burst duration (min');
xlim([0 100]);
%%
% Temperature 
figure
errorbar(Egglength,AvgDurAllAP(1).AvgDur,AvgDurAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(11).AvgDur,AvgDurAllAP(11).All95Conf,'Color',Colors(11).Color,'LineWidth',2.5,'LineStyle','-.');
%legend('Distal','Distal 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
print( [FigDirect filesep 'DistalTCompDuration'],'-dsvg');

figure
errorbar(Egglength,AvgDurAllAP(2).AvgDur,AvgDurAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(12).AvgDur,AvgDurAllAP(12).All95Conf,'Color',Colors(12).Color,'LineWidth',2.5,'LineStyle','-.');
%legend('Distal','Distal 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
print( [FigDirect filesep 'ProxTCompDuration'],'-dsvg');

figure
errorbar(Egglength,AvgDurAllAP(3).AvgDur,AvgDurAllAP(3).All95Conf,'Color',Colors(3).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(13).AvgDur,AvgDurAllAP(13).All95Conf,'Color',Colors(13).Color,'LineWidth',2.5,'LineStyle','-.');
legend('BothSep','BothSep 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Burst duration (min)');
xlim([0 100]);

figure
errorbar(Egglength,AvgDurAllAP(9).AvgDur,AvgDurAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(14).AvgDur,AvgDurAllAP(14).All95Conf,'Color',Colors(14).Color,'LineWidth',2.5,'LineStyle','-.');
%legend('Both','Both 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
print( [FigDirect filesep 'BothTCompDuration'],'-dsvg');

figure
errorbar(Egglength,AvgDurAllAP(7).AvgDur,AvgDurAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(15).AvgDur,AvgDurAllAP(15).All95Conf,'Color',Colors(15).Color,'LineWidth',2.5,'LineStyle','-.');
%legend('2x Proximal','2x Proximal 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
print( [FigDirect filesep '2xProxTCompDuration'],'-dsvg');
%% 17C Comp
figure
errorbar(Egglength,AvgDurAllAP(1).AvgDur,AvgDurAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(11).AvgDur,AvgDurAllAP(11).All95Conf,'Color',Colors(11).Color,'LineWidth',2.5,'LineStyle','-.');
%legend('Distal','Distal 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
print( [FigDirect filesep 'Dist17CCompDuration'],'-dsvg');

figure
errorbar(Egglength,AvgDurAllAP(9).AvgDur,AvgDurAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength,AvgDurAllAP(17).AvgDur,AvgDurAllAP(17).All95Conf,'Color',Colors(17).Color,'LineWidth',2.5,'LineStyle','-.');
%legend('Distal','Distal 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('burst duration (min)');
xlim([0 100]);
%%


%save duration info
%save([DropboxFolder filesep 'Constructs' filesep 'BurstDuration'],'AvgDurAllAP');
%save([DropboxFolder filesep 'Constructs' filesep 'BurstDurationSlope'],'AvgDurAllAP');
%% 2way anova of Duration vs AP position vs genotype 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
%ncUse=input('Want to only use nc14?','s');

%Count for each construct

DurationManovaVect=[];
APManovaVect=[];
ConManovaVect=[];
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
        for nn=1:length(BurstProperties)
            DurationManovaVect=[DurationManovaVect, nanmean([BurstProperties(nn).Duration])];
            APManovaVect=[APManovaVect, BurstProperties(nn).APBin];
            ConManovaVect=[ConManovaVect, cc];
        end
    end
end

[p,tbl,stats]=anovan(DurationManovaVect,{APManovaVect, ConManovaVect},'sstype',1,'varnames',{'AP bin','Construct'})
%multcompare(stats,'Dimension',2)


% %plot constructs against one another as bar graphs at a single AP position
% figure
% for ii=1:length(AvgDurAllAP)
%     AvgDuratAP(ii)=AvgDurAllAP(ii).AvgDur(APtoUse);
% 
%     errorbar(ii,AvgDuratAP(ii), AvgDurAllAP(ii).ConSE(APtoUse),'o');
%     hold on
% end
% bar(AvgDuratAP,'EdgeColor','k','LineWidth',1.5)
% xticks([1:6]);
% xticklabels({'Dist','Prox','BothSep','2xDist','2xProx','Both'});
% xlim([0 7])
% xlabel('Construct')
% ylabel('Burst length (min)');
% ylim([0 10]);
% title(['Mean burst duration',' ',num2str(EgglengthUse),'% of egg length']);
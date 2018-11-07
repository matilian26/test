ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'HbEmpty';'Kr2xProxEmpty';'KrDist17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
SlopeUse=input('Want to use slope calculations?','s');
Halves=input('Do halves calculations?','s');
%if Halves ~='y'
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
        %AmplitudeSE(ee,aa,cc)=AmplitudeSD(ee,aa,cc)/sqrt(length(AmplitudeAP));   
        AmplitudeSE(ee,aa,cc)=AmplitudeSD(ee,aa,cc)/sqrt(sum(~isnan(AmpAP)));
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
        AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp=nanmean(AmpAllAP(:,:,ee));
        AvgAmpAllAP(cc).EmbryoAmp(ee).SD=AmplitudeSD(ee,:,cc);
        AvgAmpAllAP(cc).EmbryoAmp(ee).SE=AmplitudeSE(ee,:,cc);
        AvgAmpAllAP(cc).EmbryoAmp(ee).Name=PrefixName;
    end
        AvgAmpAllAP(cc).AvgAmp=nanmean(ConAmpAllAP);  %Avg duration of all embryos of a construct
        AvgAmpAllAP(cc).AmpSD=nanmean(ConAmpSD);
        AvgAmpAllAP(cc).AmpSE=nanmean(ConAmpSE);
        AvgAmpAllAP(cc).AllAmps=[ConAmpAllAP];
        AvgAmpAllAP(cc).AllSD=nanstd([AvgAmpAllAP(cc).AllAmps]);
        %AvgAmpAllAP(cc).AllSE=((AvgAmpAllAP(cc).AllSD)/(sqrt(length(AvgAmpAllAP(cc).AllAmps))));
        for aa=1:length(APbinID)
        AvgAmpAllAP(cc).AllSE(aa)=(AvgAmpAllAP(cc).AllSD(aa))/(sqrt(sum(~isnan(AvgAmpAllAP(cc).AllAmps(:,aa)))));
        end
        AvgAmpAllAP(cc).All95Conf=(AvgAmpAllAP(cc).AllSE) .*1.95;

end

% Get rid of points that are only one data point for the whole AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(AvgAmpAllAP(cc).AllAmps(:,aa))) ==1
            AvgAmpAllAP(cc).AllSD(aa)=nan;
            AvgAmpAllAP(cc).AvgAmp(aa)=nan;
            AvgAmpAllAP(cc).AllSE(aa)=nan;
            AvgAmpAllAP(cc).All95Conf(aa)=nan;
        end
    end
end


%% ploting
%if Halves ~='y'
%Plot mean of each embryo for each construct
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
HbColor='k';

Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DistalEmptyColor;
Colors(5).Color=ProxEmptyColor;
Colors(6).Color=DoubDistColor;
Colors(7).Color=DoubProxColor;
Colors(8).Color=BothColor;
Colors(9).Color=BothEmptyColor;
Colors(10).Color=DistalColor;
Colors(11).Color=ProxColor;
Colors(12).Color=BothColor;
Colors(13).Color=BothColor;
Colors(14).Color=DoubProxColor;
Colors(15).Color=HbColor;
Colors(16).Color=DoubProxColor;
Colors(17).Color=DistalColor;

fontsize=18;
fontname='Helvetica';
FigDirect=[DropboxFolder filesep 'Figures'];

Egglength=APbinID .*100;
MS2Convers=377;

%%
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for ee=1:length(AvgAmpAllAP(cc).EmbryoAmp)
        APBoxing(ee,cc)=AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp(APtoUse);
    end
end
APBoxing(APBoxing==0)=nan;
figure 
MaxFluoVal=0;
for cc=1:length(ConstructList)
    for ee=1:length(AvgAmpAllAP(cc).EmbryoAmp)
        yyaxis left
        plot(cc,AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp(APtoUse),'o','LineWidth',1.5,'Color',Colors(cc).Color)
        yyaxis right
        plot(cc, (AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp(APtoUse)/MS2Convers),'o','LineWidth',1.5,'Color',Colors(cc).Color);
        hold on
        if AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp(APtoUse) > MaxFluoVal
            MaxFluoVal=AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp(APtoUse);
        end
    end
end
yyaxis left    
    boxplot(APBoxing,'Colors','k')
xlim([0 16]);
xlabel('Construct');
xticks([1:15]);
xticklabels({'Dist', 'Prox', 'Both Sep','1x Dist', '1x Prox', '2x Dist', '2x Prox','Both','1x Both','Dist32C','BothSep32C','Both32C','2xProx32C'});
ylabel('Fluorescence intensity (AU)');
title(['Mean burst amplitude',' ' ,num2str(EgglengthUse),'% egg length'])
ylim([0 (MaxFluoVal+1000)]);
yyaxis right
ylabel('# RNAP mol');
ylim([0 ((MaxFluoVal+1000)/MS2Convers)]);
 


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
xticks([1:15]);
xticklabels({'Dist','Prox','BothSep','1xDist','1xProx','2xDist','2xProx','Both','1xBoth','Dist32','BothSep32'});
ylabel('Fluorescent intensity');
title(['Avg burst amplitude AP bin',' ',num2str(APtoUse)]);

%Plot construct means vs AP position
figure 
for cc=1:length(ConstructList)
    yyaxis left
    plot(Egglength,AvgAmpAllAP(cc).AvgAmp,'Color',Colors(cc).Color,'LineWidth',1.5,'LineStyle','-');
    hold on 
    yyaxis right
    plot(Egglength,(AvgAmpAllAP(cc).AvgAmp./MS2Convers),'Color',Colors(cc).Color,'LineWidth',1.5,'LineStyle','-');
end
yyaxis left
legend('Dist', 'Prox', 'Both Sep', '1x Dist','1x Prox', '2x Dist', '2x Prox', 'Both','1x Both');
xlabel('% Egg length');
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);
ylim([0 50000]);
yyaxis right 
ylabel('# active PolII molecules');
ylim([0 (50000/MS2Convers)]);

%%
%doubles vs SE
figure 
%errorbar(Egglength, AvgAmpAllAP(6).AvgAmp, AvgAmpAllAP(6).All95Conf, 'Color', Colors(6).Color, 'LineWidth', 1.5);
hold on
yyaxis left
errorbar(Egglength, AvgAmpAllAP(7).AvgAmp, AvgAmpAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth', 2.5);
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
%legend('2x Prox', 'Both');
xlabel('% Egg length');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Fluorescence intensity');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 45000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(7).AvgAmp./MS2Convers), (AvgAmpAllAP(7).All95Conf./MS2Convers), 'Color', Colors(7).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(8).AvgAmp./MS2Convers), (AvgAmpAllAP(8).All95Conf./MS2Convers), 'Color', Colors(8).Color, 'LineWidth', 2.5,'LineStyle','-');
ylabel('Number PolII molecules');
ylim([0 (45000/MS2Convers)]);
legend('2x Prox', 'Both');

figure 
errorbar(Egglength, AvgAmpAllAP(6).AvgAmp, AvgAmpAllAP(6).All95Conf, 'Color', Colors(6).Color, 'LineWidth', 2.5);
hold on
yyaxis left
errorbar(Egglength, AvgAmpAllAP(7).AvgAmp, AvgAmpAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth', 2.5);
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
%legend('2x Prox', 'Both');
xlabel('% Egg length');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Fluorescence intensity');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 25000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(6).AvgAmp./MS2Convers), (AvgAmpAllAP(6).All95Conf./MS2Convers),'Color',Colors(6).Color,'LineWidth',2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(7).AvgAmp./MS2Convers), (AvgAmpAllAP(7).All95Conf./MS2Convers), 'Color', Colors(7).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(8).AvgAmp./MS2Convers), (AvgAmpAllAP(8).All95Conf./MS2Convers), 'Color', Colors(8).Color, 'LineWidth', 2.5,'LineStyle','-');
ylabel('Number PolII molecules');
ylim([0 (25000/MS2Convers)]);



%%
%Singles vs each other or SE
figure 
yyaxis left
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
%legend('Distal', 'Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname,'YColor','k');
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 25000]);
yyaxis right
set(gca,'YColor','k');
errorbar(Egglength, (AvgAmpAllAP(1).AvgAmp./MS2Convers), (AvgAmpAllAP(1).All95Conf./MS2Convers), 'Color', Colors(1).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(2).AvgAmp./MS2Convers), (AvgAmpAllAP(2).All95Conf./MS2Convers), 'Color', Colors(2).Color, 'LineWidth', 2.5,'LineStyle','-');
ylabel('number polII molecules');
ylim([0 (25000/MS2Convers)]);
print( [FigDirect filesep 'SinglesAmpmRNA'],'-dsvg');

%%
figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(4).AvgAmp, AvgAmpAllAP(4).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5,'LineStyle',':');
%legend('Distal', '1x Distal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 30000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(1).AvgAmp./MS2Convers), (AvgAmpAllAP(1).All95Conf./MS2Convers), 'Color', Colors(1).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(4).AvgAmp./MS2Convers), (AvgAmpAllAP(4).All95Conf./MS2Convers), 'Color', Colors(1).Color, 'LineWidth', 2.5,'LineStyle',':');
ylabel('number polII molecules');
ylim([0 (30000/MS2Convers)]);
print( [FigDirect filesep 'HemiDistAmp'],'-dsvg');

figure 
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(5).AvgAmp, AvgAmpAllAP(5).All95Conf, 'Color', Colors(5).Color, 'LineWidth', 2.5,'LineStyle',':');
%legend('Proximal', '1x Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 30000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(2).AvgAmp./MS2Convers), (AvgAmpAllAP(2).All95Conf./MS2Convers), 'Color', Colors(2).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(5).AvgAmp./MS2Convers), (AvgAmpAllAP(5).All95Conf./MS2Convers), 'Color', Colors(2).Color, 'LineWidth', 2.5,'LineStyle',':');
ylabel('number polII molecules');
ylim([0 (30000/MS2Convers)]);
print( [FigDirect filesep 'HemiProxAmp'],'-dsvg');

figure 
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(9).AvgAmp, AvgAmpAllAP(9).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5,'LineStyle',':');
%legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 30000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(8).AvgAmp./MS2Convers), (AvgAmpAllAP(8).All95Conf./MS2Convers), 'Color', Colors(8).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(9).AvgAmp./MS2Convers), (AvgAmpAllAP(9).All95Conf./MS2Convers), 'Color', Colors(8).Color, 'LineWidth', 2.5,'LineStyle',':');
ylabel('number polII molecules');
ylim([0 (30000/MS2Convers)]);
print( [FigDirect filesep 'HemiBothAmp'],'-dsvg');

figure 
errorbar(Egglength, AvgAmpAllAP(7).AvgAmp, AvgAmpAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(16).AvgAmp, AvgAmpAllAP(16).All95Conf, 'Color', Colors(16).Color, 'LineWidth', 2.5,'LineStyle',':');
%legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 30000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(7).AvgAmp./MS2Convers), (AvgAmpAllAP(7).All95Conf./MS2Convers), 'Color', Colors(7).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(16).AvgAmp./MS2Convers), (AvgAmpAllAP(16).All95Conf./MS2Convers), 'Color', Colors(16).Color, 'LineWidth', 2.5,'LineStyle',':');
ylabel('number polII molecules');
ylim([0 (30000/MS2Convers)]);
print( [FigDirect filesep 'Hemi2xProxAmp'],'-dsvg');
%%
figure 
yyaxis left
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 1.5);
hold on
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 1.5);
errorbar(Egglength, AvgAmpAllAP(3).AvgAmp, AvgAmpAllAP(3).All95Conf, 'Color', Colors(3).Color, 'LineWidth', 1.5);
xlabel('% Egg length');
ylabel('Fluorescence intensity');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 45000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(1).AvgAmp./MS2Convers), (AvgAmpAllAP(1).All95Conf./MS2Convers), 'Color', Colors(1).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(2).AvgAmp./MS2Convers), (AvgAmpAllAP(2).All95Conf./MS2Convers), 'Color', Colors(2).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(3).AvgAmp./MS2Convers), (AvgAmpAllAP(3).All95Conf./MS2Convers), 'Color', Colors(3).Color, 'LineWidth', 2.5,'LineStyle','-');
ylabel('number polII molecules');
ylim([0 (45000/MS2Convers)]);
legend('Distal', 'Proximal','Both Sep');
set(gca, 'FontSize', fontsize, 'FontName', fontname);

figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 25000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(1).AvgAmp./MS2Convers), (AvgAmpAllAP(1).All95Conf./MS2Convers), 'Color', Colors(1).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(2).AvgAmp./MS2Convers), (AvgAmpAllAP(2).All95Conf./MS2Convers), 'Color', Colors(2).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(8).AvgAmp./MS2Convers), (AvgAmpAllAP(8).All95Conf./MS2Convers), 'Color', Colors(8).Color, 'LineWidth', 2.5,'LineStyle','-');
ylabel('number polII molecules','Color','k');
set(gca,'YColor','k');
ylim([0 (25000/MS2Convers)]);
%legend('Distal', 'Proximal','Both');
print( [FigDirect filesep 'SinglesvBothAmp'],'-dsvg');



%%
%singles vs doubles
figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).AmpSE, 'Color', Colors(1).Color, 'LineWidth', 1.5);
hold on
errorbar(Egglength, AvgAmpAllAP(6).AvgAmp, AvgAmpAllAP(6).AmpSE, 'Color', Colors(6).Color, 'LineWidth', 1.5);
legend('Dist', '2x Dist');
xlabel('% Egg length');
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).AmpSE, 'Color', Colors(2).Color, 'LineWidth', 1.5);
hold on
errorbar(Egglength, AvgAmpAllAP(7).AvgAmp, AvgAmpAllAP(7).AmpSE, 'Color', Colors(7).Color, 'LineWidth', 1.5);
%legend('Prox', '2x Prox');
xlabel('% egg length');
ylabel('fluorescence intensity');
%title('Mean burst amplitude');
xlim([0 100]);
print( [FigDirect filesep 'ProxvsDuplicAmp'],'-dsvg');


figure 
errorbar(Egglength, AvgAmpAllAP(3).AvgAmp, AvgAmpAllAP(3).AmpSE, 'Color', Colors(3).Color, 'LineWidth', 1.5);
hold on
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).AmpSE, 'Color', Colors(8).Color, 'LineWidth', 1.5);
legend('Both Sep', 'Both');
xlabel('% Egg length');
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);
%%
% single allele stuff 
figure 
errorbar(Egglength, AvgAmpAllAP(4).AvgAmp, AvgAmpAllAP(4).All95Conf, 'Color', Colors(4).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
legend('1x Distal', 'Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(5).AvgAmp, AvgAmpAllAP(5).All95Conf, 'Color', Colors(5).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
legend('1x Proximal', 'Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
%%
% 32C vs room temp
figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(10).AvgAmp, AvgAmpAllAP(10).All95Conf, 'Color', Colors(10).Color, 'LineWidth', 2.5,'LineStyle','-.');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 25000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(1).AvgAmp./MS2Convers), (AvgAmpAllAP(1).All95Conf./MS2Convers), 'Color', Colors(1).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(10).AvgAmp./MS2Convers), (AvgAmpAllAP(10).All95Conf./MS2Convers), 'Color', Colors(10).Color, 'LineWidth', 2.5,'LineStyle','-.');
ylabel('number polII molecules','Color','k');
ylim([0 (25000/MS2Convers)]);
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep 'DistalTempCompAmp'],'-dsvg');


figure 
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(11).AvgAmp, AvgAmpAllAP(11).All95Conf, 'Color', Colors(11).Color, 'LineWidth', 2.5,'LineStyle','-.');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 25000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(2).AvgAmp./MS2Convers), (AvgAmpAllAP(2).All95Conf./MS2Convers), 'Color', Colors(2).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(11).AvgAmp./MS2Convers), (AvgAmpAllAP(11).All95Conf./MS2Convers), 'Color', Colors(11).Color, 'LineWidth', 2.5,'LineStyle','-.');
ylabel('number polII molecules','Color','k');
ylim([0 (25000/MS2Convers)]);
%legend('Proximal', 'Proximal 32C','Location','best');
print( [FigDirect filesep 'ProxTCompAmp'],'-dsvg');


figure 
errorbar(Egglength, AvgAmpAllAP(3).AvgAmp, AvgAmpAllAP(3).All95Conf, 'Color', Colors(3).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(12).AvgAmp, AvgAmpAllAP(12).All95Conf, 'Color', Colors(12).Color, 'LineWidth', 2.5,'LineStyle','-.');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 25000]);
yyaxis right 
errorbar(Egglength, (AvgAmpAllAP(3).AvgAmp./MS2Convers), (AvgAmpAllAP(3).All95Conf./MS2Convers), 'Color', Colors(3).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(12).AvgAmp./MS2Convers), (AvgAmpAllAP(12).All95Conf./MS2Convers), 'Color', Colors(12).Color, 'LineWidth', 2.5,'LineStyle','-.');
ylabel('number polII molecules','Color','k');
ylim([0 (25000/MS2Convers)]);
%legend('Both Sep', 'Both Sep 32C','Location','Best');
print( [FigDirect filesep 'SepTCompAmp'],'-dsvg');


figure 
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(13).AvgAmp, AvgAmpAllAP(13).All95Conf, 'Color', Colors(13).Color, 'LineWidth', 2.5,'LineStyle','-.');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
xlim([0 100]);
ylim([0 25000]);
yyaxis right 
errorbar(Egglength, (AvgAmpAllAP(8).AvgAmp./MS2Convers), (AvgAmpAllAP(8).All95Conf./MS2Convers), 'Color', Colors(8).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(13).AvgAmp./MS2Convers), (AvgAmpAllAP(13).All95Conf./MS2Convers), 'Color', Colors(13).Color, 'LineWidth', 2.5,'LineStyle','-.');
ylabel('number polII molecules','Color','k');
ylim([0 (25000/MS2Convers)]);
%legend('both', 'both 32C','Location','Best');
print( [FigDirect filesep 'BothTCompAmp'],'-dsvg');


figure 
errorbar(Egglength, AvgAmpAllAP(7).AvgAmp, AvgAmpAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(14).AvgAmp, AvgAmpAllAP(14).All95Conf, 'Color', Colors(14).Color, 'LineWidth', 2.5,'LineStyle','-.');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
xlim([0 100]);
ylim([0 25000]);
yyaxis right 
errorbar(Egglength, (AvgAmpAllAP(7).AvgAmp./MS2Convers), (AvgAmpAllAP(7).All95Conf./MS2Convers), 'Color', Colors(7).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(14).AvgAmp./MS2Convers), (AvgAmpAllAP(14).All95Conf./MS2Convers), 'Color', Colors(14).Color, 'LineWidth', 2.5,'LineStyle','-.');
ylabel('number polII molecules','Color','k');
ylim([0 (25000/MS2Convers)]);
%legend('2x Proximal', '2xProximal 32C','Location','best');
print( [FigDirect filesep '2xProxTCompAmp'],'-dsvg');

%% 17C vs RT
figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(17).AvgAmp, AvgAmpAllAP(17).All95Conf, 'Color', Colors(17).Color, 'LineWidth', 2.5,'LineStyle','-.');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);
ylim([0 25000]);
yyaxis right
errorbar(Egglength, (AvgAmpAllAP(1).AvgAmp./MS2Convers), (AvgAmpAllAP(1).All95Conf./MS2Convers), 'Color', Colors(1).Color, 'LineWidth', 2.5,'LineStyle','-');
errorbar(Egglength, (AvgAmpAllAP(10).AvgAmp./MS2Convers), (AvgAmpAllAP(10).All95Conf./MS2Convers), 'Color', Colors(10).Color, 'LineWidth', 2.5,'LineStyle','-.');
ylabel('number polII molecules','Color','k');
ylim([0 (25000/MS2Convers)]);
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep 'Dist17CCompAmp'],'-dsvg');
%%
figure
yyaxis left
errorbar(Egglength, AvgAmpAllAP(15).AvgAmp, AvgAmpAllAP(15).All95Conf, 'Color',Colors(15).Color,'LineWidth',2.5);
ylabel('Fluorescence intensity (AU)');
MaxVal=max(AvgAmpAllAP(14).AvgAmp);
ylim([0 (MaxVal+5000)]);
yyaxis right 
errorbar(Egglength, (AvgAmpAllAP(15).AvgAmp./MS2Convers), (AvgAmpAllAP(15).All95Conf./MS2Convers),'Color',Colors(15).Color,'LineWidth',2.5);
ylabel('# active RNAP molecules');
ylim([0 ((MaxVal+5000)/MS2Convers)]);


%Save amp info 
save([DropboxFolder filesep 'Constructs' filesep 'BurstAmplitude'],'AvgAmpAllAP');

%% 2way anova of Duration vs AP position vs genotype 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
%ncUse=input('Want to only use nc14?','s');

%Count for each construct

AmpManovaVect=[];
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
            AmpManovaVect=[AmpManovaVect, nanmean([BurstProperties(nn).BurstAmplitude])];
            APManovaVect=[APManovaVect, BurstProperties(nn).APBin];
            ConManovaVect=[ConManovaVect, cc];
        end
    end
end

[p,tbl,stats]=anovan(AmpManovaVect,{APManovaVect, ConManovaVect},'sstype',1,'varnames',{'AP bin','Construct'})
%multcompare(stats,'Dimension',2)




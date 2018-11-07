ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty';'KrDist32C';'KrBothSep32';'KrBoth32C';'Kr2xProx32C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
SlopeUse=input('Want to use slope calculations?','s');
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
    for ee=1:2:NEmbryos
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
        AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp=nanmean(AmpAllAP(:,:,ee));
        AvgAmpAllAP(cc).EmbryoAmp(ee).SD=AmplitudeSD(ee,:,cc);
        AvgAmpAllAP(cc).EmbryoAmp(ee).SE=AmplitudeSE(ee,:,cc);
    end
        AvgAmpAllAP(cc).AvgAmp=nanmean(ConAmpAllAP);  %Avg duration of all embryos of a construct
        AvgAmpAllAP(cc).AmpSD=nanmean(ConAmpSD);
        AvgAmpAllAP(cc).AmpSE=nanmean(ConAmpSE);
        AvgAmpAllAP(cc).AllAmps=[ConAmpAllAP];
        AvgAmpAllAP(cc).AllSD=nanstd([AvgAmpAllAP(cc).AllAmps]);
        AvgAmpAllAP(cc).AllSE=((AvgAmpAllAP(cc).AllSD)/(sqrt(length(AvgAmpAllAP(cc).AllAmps))));
        AvgAmpAllAP(cc).All95Conf=(AvgAmpAllAP(cc).AllSE) .*1.95;

end
% Do for 2nd half of data 
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
    for ee=2:2:NEmbryos
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
        AvgAmpAllAP2(cc).EmbryoAmp(ee).MeanAmp=nanmean(AmpAllAP(:,:,ee));
        AvgAmpAllAP2(cc).EmbryoAmp(ee).SD=AmplitudeSD(ee,:,cc);
        AvgAmpAllAP2(cc).EmbryoAmp(ee).SE=AmplitudeSE(ee,:,cc);
    end
        AvgAmpAllAP2(cc).AvgAmp=nanmean(ConAmpAllAP);  %Avg duration of all embryos of a construct
        AvgAmpAllAP2(cc).AmpSD=nanmean(ConAmpSD);
        AvgAmpAllAP2(cc).AmpSE=nanmean(ConAmpSE);
        AvgAmpAllAP2(cc).AllAmps=[ConAmpAllAP];
        AvgAmpAllAP2(cc).AllSD=nanstd([AvgAmpAllAP2(cc).AllAmps]);
        AvgAmpAllAP2(cc).AllSE=((AvgAmpAllAP2(cc).AllSD)/(sqrt(length(AvgAmpAllAP2(cc).AllAmps))));
        AvgAmpAllAP2(cc).All95Conf=(AvgAmpAllAP2(cc).AllSE) .*1.95;

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

% Do for 2nd half as well
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(AvgAmpAllAP2(cc).AllAmps(:,aa))) ==1
            AvgAmpAllAP2(cc).AllSD(aa)=nan;
            AvgAmpAllAP2(cc).AvgAmp(aa)=nan;
            AvgAmpAllAP2(cc).AllSE(aa)=nan;
            AvgAmpAllAP2(cc).All95Conf(aa)=nan;
        end
    end
end


%% ploting
%Plot mean of each embryo for each construct
DistalColor=[8 180 238] ./ 255;
DistalEmptyColor=[8 210 238] ./ 255; 
Distal32CColor=[118 180 238] ./ 255;
DoubDistColor=[1 17 181] ./ 255;
ProxColor=[251 220 50] ./ 255;
ProxEmptyColor=[251 250 50] ./255;
DoubProxColor=[251 190 100] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[12 195 82] ./ 255;
BothEmptyColor=[12 250 150] ./ 255;
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
Colors(10).Color=Distal32CColor;
Colors(11).Color=BothSep32CColor;
Colors(12).Color=Both32CColor;
Colors(13).Color=DoubProx32CColor;

fontsize=18;
fontname='Helvetica';

Egglength=APbinID .*100;
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for ee=1:2:length(AvgAmpAllAP(cc).EmbryoAmp)
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
        plot(cc, (AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp(APtoUse)/350),'o','LineWidth',1.5,'Color',Colors(cc).Color);
        hold on
        if AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp(APtoUse) > MaxFluoVal
            MaxFluoVal=AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp(APtoUse);
        end
        
    end
end
yyaxis left    
    boxplot(APBoxing,'Colors','k')
xlim([0 14]);
xlabel('Construct');
xticks([1:13]);
xticklabels({'Dist', 'Prox', 'Both Sep','1x Dist', '1x Prox', '2x Dist', '2x Prox','Both','1x Both','Dist32C','BothSep32C','Both32C','2xProx32C'});
title(['Mean burst amplitude',' ' ,num2str(EgglengthUse),'% egg length'])
yyaxis left
ylabel('Fluorescence intensity (AU)');
ylim([0 (MaxFluoVal+1000)]);
yyaxis right
ylabel('# RNAP mol');
ylim([0 ((MaxFluoVal+1000)/350)]);

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
xticks([1:11]);
xticklabels({'Dist','Prox','BothSep','1xDist','1xProx','2xDist','2xProx','Both','1xBoth','Dist32','BothSep32'});
ylabel('Fluorescent intensity');
title(['Avg burst amplitude AP bin',' ',num2str(APtoUse)]);

%Plot construct means vs AP position
figure 
for cc=1:length(ConstructList)
    yyaxis left
    plot(Egglength,AvgAmpAllAP(cc).AvgAmp,'Color',Colors(cc).Color,'LineWidth',1.5);
    yyaxis right
    plot(Egglength,((AvgAmpAllAP(cc).AvgAmp)/350),'Color',Colors(cc).Color,'LineWidth',1.5);
    hold on 
end
yyaxis left
cla
legend('Dist', 'Prox', 'Both Sep', '1x Dist','1x Prox', '2x Dist', '2x Prox', 'Both','1x Both');
xlabel('% Egg length');
title('Mean burst amplitude');
xlim([0 100]);
yyaxis right
ylabel('# active RNAP mol');

%doubles vs SE
figure 
%errorbar(Egglength, AvgAmpAllAP(6).AvgAmp, AvgAmpAllAP(6).All95Conf, 'Color', Colors(6).Color, 'LineWidth', 1.5);
hold on
errorbar(Egglength, AvgAmpAllAP(7).AvgAmp, AvgAmpAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth', 2.5);
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
legend('2x Prox', 'Both');
xlabel('% Egg length');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(9).AvgAmp, AvgAmpAllAP(9).All95Conf, 'Color', Colors(9).Color, 'LineWidth', 2.5);
legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

%Singles vs each other or SE
figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
%legend('Distal', 'Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(4).AvgAmp, AvgAmpAllAP(4).All95Conf, 'Color', Colors(4).Color, 'LineWidth', 2.5);
legend('Distal', '1x Distal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(5).AvgAmp, AvgAmpAllAP(5).All95Conf, 'Color', Colors(5).Color, 'LineWidth', 2.5);
legend('Proximal', '1x Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 1.5);
hold on
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 1.5);
errorbar(Egglength, AvgAmpAllAP(3).AvgAmp, AvgAmpAllAP(3).All95Conf, 'Color', Colors(3).Color, 'LineWidth', 1.5);
legend('Distal', 'Proximal','Both Sep');
xlabel('% Egg length');
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
legend('Distal', 'Proximal','Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

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
legend('Prox', '2x Prox');
xlabel('% Egg length');
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(3).AvgAmp, AvgAmpAllAP(3).AmpSE, 'Color', Colors(3).Color, 'LineWidth', 1.5);
hold on
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).AmpSE, 'Color', Colors(8).Color, 'LineWidth', 1.5);
legend('Both Sep', 'Both');
xlabel('% Egg length');
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);

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

% 32C vs room temp
figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(10).AvgAmp, AvgAmpAllAP(10).All95Conf, 'Color', Colors(10).Color, 'LineWidth', 2.5);
legend('Distal', 'Distal 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(3).AvgAmp, AvgAmpAllAP(3).All95Conf, 'Color', Colors(3).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(11).AvgAmp, AvgAmpAllAP(11).All95Conf, 'Color', Colors(11).Color, 'LineWidth', 2.5);
legend('Both Sep', 'Both Sep 32C','Location','Best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(12).AvgAmp, AvgAmpAllAP(12).All95Conf, 'Color', Colors(12).Color, 'LineWidth', 2.5);
legend('Both', 'Both 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(7).AvgAmp, AvgAmpAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(13).AvgAmp, AvgAmpAllAP(13).All95Conf, 'Color', Colors(13).Color, 'LineWidth', 2.5);
legend('2x Proximal', '2xProximal 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
xlim([0 100]);

%Save amp info 
save([DropboxFolder filesep 'Constructs' filesep 'BurstAmplitude'],'AvgAmpAllAP');
%plot constructs against one another as bar graphs at a single AP position
% figure
% for ii=1:length(AvgAmpAllAP)
%     AvgAmpatAP(ii)=AvgAmpAllAP(ii).AvgAmp(APtoUse);
% 
%     errorbar(ii,AvgAmpatAP(ii), AvgAmpAllAP(ii).AmpSE(APtoUse),'o');
%     hold on
% end
% bar(AvgAmpatAP,'EdgeColor','k','LineWidth',1.5)
% xticks([1:6]);
% xticklabels({'Dist','Prox','BothSep','2xDist','2xProx','Both'});
% xlim([0 7])
% xlabel('Construct')
% ylabel('Fluorescent intensity (AU)');
% ylim([0 50000]);
% title(['Mean burst amplitude',' ',num2str(EgglengthUse),'% egg length']);

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


%% run 2nd half of embryos
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
    for ee=2:2:NEmbryos
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
        AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp=nanmean(AmpAllAP(:,:,ee));
        AvgAmpAllAP(cc).EmbryoAmp(ee).SD=AmplitudeSD(ee,:,cc);
        AvgAmpAllAP(cc).EmbryoAmp(ee).SE=AmplitudeSE(ee,:,cc);
    end
        AvgAmpAllAP(cc).AvgAmp=nanmean(ConAmpAllAP);  %Avg duration of all embryos of a construct
        AvgAmpAllAP(cc).AmpSD=nanmean(ConAmpSD);
        AvgAmpAllAP(cc).AmpSE=nanmean(ConAmpSE);
        AvgAmpAllAP(cc).AllAmps=[ConAmpAllAP];
        AvgAmpAllAP(cc).AllSD=nanstd([AvgAmpAllAP(cc).AllAmps]);
        AvgAmpAllAP(cc).AllSE=((AvgAmpAllAP(cc).AllSD)/(sqrt(length(AvgAmpAllAP(cc).AllAmps))));
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
%Plot mean of each embryo for each construct
DistalColor=[8 180 238] ./ 255;
DistalEmptyColor=[8 210 238] ./ 255; 
Distal32CColor=[118 180 238] ./ 255;
DoubDistColor=[1 17 181] ./ 255;
ProxColor=[251 220 50] ./ 255;
ProxEmptyColor=[251 250 50] ./255;
DoubProxColor=[251 190 100] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[12 195 82] ./ 255;
BothEmptyColor=[12 250 150] ./ 255;
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
Colors(10).Color=Distal32CColor;
Colors(11).Color=BothSep32CColor;
Colors(12).Color=Both32CColor;
Colors(13).Color=DoubProx32CColor;

fontsize=18;
fontname='Helvetica';

Egglength=APbinID .*100;
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for ee=2:2:length(AvgAmpAllAP(cc).EmbryoAmp)
        APBoxing(ee,cc)=AvgAmpAllAP(cc).EmbryoAmp(ee).MeanAmp(APtoUse);
    end
end
APBoxing(APBoxing==0)=nan;
figure 

    for ee=2:2:length(AvgAmpAllAP(1).EmbryoAmp)
    plot(1,AvgAmpAllAP(1).EmbryoAmp(ee).MeanAmp(APtoUse),'o','LineWidth',1.5,'Color',DistalColor)
    hold on
    end
    for ee=2:2:length(AvgAmpAllAP(2).EmbryoAmp)
    plot(2,AvgAmpAllAP(2).EmbryoAmp(ee).MeanAmp(APtoUse),'o','LineWidth',1.5,'Color',ProxColor)
        
    hold on
    end
    
    for ee=2:2:length(AvgAmpAllAP(3).EmbryoAmp)
    plot(3,AvgAmpAllAP(3).EmbryoAmp(ee).MeanAmp(APtoUse),'o','LineWidth',1.5,'Color',BothSepColor)
    hold on
    end
    for ee=2:2:length(AvgAmpAllAP(4).EmbryoAmp)
    plot(4,AvgAmpAllAP(4).EmbryoAmp(ee).MeanAmp(APtoUse),'o','LineWidth',1.5,'Color',Colors(4).Color)
    hold on
    end
    for ee=2:2:length(AvgAmpAllAP(5).EmbryoAmp)
    plot(5,AvgAmpAllAP(5).EmbryoAmp(ee).MeanAmp(APtoUse),'o','LineWidth',1.5,'Color',Colors(5).Color)
    hold on
    end
    for ee=2:2:length(AvgAmpAllAP(6).EmbryoAmp)
    plot(6,AvgAmpAllAP(6).EmbryoAmp(ee).MeanAmp(APtoUse),'o','LineWidth',1.5,'Color',Colors(6).Color);
    %errorbar(6,AvgAmpAllAP(6).EmbryoAmp(ee).MeanProd(APtoUse),AvgAmpAllAP(6).EmbryoAmp(ee).SE(APtoUse),'o','LineWidth',1.5,'Color',BothColor)
    hold on
    end
     for ee=2:2:length(AvgAmpAllAP(7).EmbryoAmp)
    plot(7,AvgAmpAllAP(7).EmbryoAmp(ee).MeanAmp(APtoUse),'o','LineWidth',1.5,'Color',Colors(7).Color)
    hold on
     end
    for ee=2:2:length(AvgAmpAllAP(8).EmbryoAmp)
    plot(8,AvgAmpAllAP(8).EmbryoAmp(ee).MeanAmp(APtoUse),'o','LineWidth',1.5,'Color',Colors(8).Color)
    hold on
    end
    for ee=2:2:length(AvgAmpAllAP(9).EmbryoAmp)
    plot(9,AvgAmpAllAP(9).EmbryoAmp(ee).MeanAmp(APtoUse),'o','LineWidth',1.5,'Color',Colors(9).Color)
    hold on
    end
    boxplot(APBoxing,'Colors','k')
xlim([0 10]);
xlabel('Construct');
xticks([1:9]);
xticklabels({'Dist', 'Prox', 'Both Sep','1x Dist', '1x Prox', '2x Dist', '2x Prox','Both','1x Both'});
ylabel('Fluorescence intensity (AU)');
title(['Mean burst amplitude',' ' ,num2str(EgglengthUse),'% egg length'])

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
xticks([1:11]);
xticklabels({'Dist','Prox','BothSep','1xDist','1xProx','2xDist','2xProx','Both','1xBoth','Dist32','BothSep32'});
ylabel('Fluorescent intensity');
title(['Avg burst amplitude AP bin',' ',num2str(APtoUse)]);

%Plot construct means vs AP position
figure 
for cc=1:length(ConstructList)
    plot(Egglength,AvgAmpAllAP(cc).AvgAmp,'Color',Colors(cc).Color,'LineWidth',1.5);
    hold on 
end
legend('Dist', 'Prox', 'Both Sep', '1x Dist','1x Prox', '2x Dist', '2x Prox', 'Both','1x Both');
xlabel('% Egg length');
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);

%doubles vs SE
figure 
%errorbar(Egglength, AvgAmpAllAP(6).AvgAmp, AvgAmpAllAP(6).All95Conf, 'Color', Colors(6).Color, 'LineWidth', 1.5);
hold on
errorbar(Egglength, AvgAmpAllAP(7).AvgAmp, AvgAmpAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth', 2.5);
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
legend('2x Prox', 'Both');
xlabel('% Egg length');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(9).AvgAmp, AvgAmpAllAP(9).All95Conf, 'Color', Colors(9).Color, 'LineWidth', 2.5);
legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

%Singles vs each other or SE
figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
%legend('Distal', 'Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(4).AvgAmp, AvgAmpAllAP(4).All95Conf, 'Color', Colors(4).Color, 'LineWidth', 2.5);
legend('Distal', '1x Distal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(5).AvgAmp, AvgAmpAllAP(5).All95Conf, 'Color', Colors(5).Color, 'LineWidth', 2.5);
legend('Proximal', '1x Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 1.5);
hold on
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 1.5);
errorbar(Egglength, AvgAmpAllAP(3).AvgAmp, AvgAmpAllAP(3).All95Conf, 'Color', Colors(3).Color, 'LineWidth', 1.5);
legend('Distal', 'Proximal','Both Sep');
xlabel('% Egg length');
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(2).AvgAmp, AvgAmpAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth', 2.5);
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
legend('Distal', 'Proximal','Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

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
legend('Prox', '2x Prox');
xlabel('% Egg length');
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(3).AvgAmp, AvgAmpAllAP(3).AmpSE, 'Color', Colors(3).Color, 'LineWidth', 1.5);
hold on
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).AmpSE, 'Color', Colors(8).Color, 'LineWidth', 1.5);
legend('Both Sep', 'Both');
xlabel('% Egg length');
ylabel('Fluorescence intensity');
title('Mean burst amplitude');
xlim([0 100]);

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

% 32C vs room temp
figure 
errorbar(Egglength, AvgAmpAllAP(1).AvgAmp, AvgAmpAllAP(1).All95Conf, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(10).AvgAmp, AvgAmpAllAP(10).All95Conf, 'Color', Colors(10).Color, 'LineWidth', 2.5);
legend('Distal', 'Distal 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(3).AvgAmp, AvgAmpAllAP(3).All95Conf, 'Color', Colors(3).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(11).AvgAmp, AvgAmpAllAP(11).All95Conf, 'Color', Colors(11).Color, 'LineWidth', 2.5);
legend('Both Sep', 'Both Sep 32C','Location','Best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
%title('Mean burst amplitude');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(8).AvgAmp, AvgAmpAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(12).AvgAmp, AvgAmpAllAP(12).All95Conf, 'Color', Colors(12).Color, 'LineWidth', 2.5);
legend('Both', 'Both 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
xlim([0 100]);

figure 
errorbar(Egglength, AvgAmpAllAP(7).AvgAmp, AvgAmpAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth', 2.5);
hold on
errorbar(Egglength, AvgAmpAllAP(13).AvgAmp, AvgAmpAllAP(13).All95Conf, 'Color', Colors(13).Color, 'LineWidth', 2.5);
legend('2x Proximal', '2xProximal 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Fluorescence intensity (AU)');
xlim([0 100]);

%Save amp info 
save([DropboxFolder filesep 'Constructs' filesep 'BurstAmplitude'],'AvgAmpAllAP');
%plot constructs against one another as bar graphs at a single AP position
% figure
% for ii=1:length(AvgAmpAllAP)
%     AvgAmpatAP(ii)=AvgAmpAllAP(ii).AvgAmp(APtoUse);
% 
%     errorbar(ii,AvgAmpatAP(ii), AvgAmpAllAP(ii).AmpSE(APtoUse),'o');
%     hold on
% end
% bar(AvgAmpatAP,'EdgeColor','k','LineWidth',1.5)
% xticks([1:6]);
% xticklabels({'Dist','Prox','BothSep','2xDist','2xProx','Both'});
% xlim([0 7])
% xlabel('Construct')
% ylabel('Fluorescent intensity (AU)');
% ylim([0 50000]);
% title(['Mean burst amplitude',' ',num2str(EgglengthUse),'% egg length']);

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
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'Kr2xProxEmpty';'KrDist17C';'KrBoth17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
SlopeUse=input('Want to use slope calculations?','s');
Halves=input('Do by halves calculations?','s');
%Count time of fluorescence turns on for each construct
if Halves ~= 'y'
AvgFirstSpot=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgFirstSpotCon=[];
    firsttime=1;
    ConFirSpotAllAP=[];
    ConFirSpotSE=[];
    ConFirSpotSD=[];
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
        FirSpotAllAP=[];
        for aa=1:length(APbinID)
            FirSpotAP=[];
            FirstSpotAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(FirstSpotAP)
                FirSpotAP=[FirSpotAP; nan];
            else
            for bb=1:length(FirstSpotAP)
                if ~isempty(BurstProperties(FirstSpotAP(bb)).Duration)
                FirSpotAP=[FirSpotAP;[BurstProperties(FirstSpotAP(bb)).FirstTimeOn]'];  %put all durations at a given AP value in a column going down
                else
                    %FirSpotAP=[FirSpotAP; 0];
                    FirSpotAP=[FirSpotAP; nan];
                end
               
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(FirSpotAP)
                FirSpotAllAP(bb,aa,ee)=FirSpotAP(bb);
            end
            FirSpotAllAP(FirSpotAllAP==0)=nan;
            
            FirstSpotSD(ee,aa,cc)=nanstd(FirSpotAllAP(:,aa,ee));
        %FirstSpotSE(ee,aa,cc)=FirstSpotSD(ee,aa,cc)/sqrt(length(FirstSpotAP));
         FirstSpotSE(ee,aa,cc)=FirstSpotSD(ee,aa,cc)/sqrt(length(FirSpotAP(~isnan(FirSpotAP))));             

            clear FirstSpotAP FirSpotAP
        end
        AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot=nanmean(FirSpotAllAP(:,:,ee));
        AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).SD=FirstSpotSD(ee,:,cc);
        AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).SE=FirstSpotSE(ee,:,cc);
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(FirSpotAllAP,3)
            ConFirSpotAllAP=[ConFirSpotAllAP; FirSpotAllAP(:,:,bb)];
        end
        for bb=1:size(FirstSpotSD,3)
            ConFirSpotSD=[ConFirSpotSD; FirstSpotSD(:,:,bb)];
        end
        for bb=1:size(FirstSpotSE,3)
            ConFirSpotSE=[ConFirSpotSE;FirstSpotSE(:,:,bb)];
        end
        
    end
        AvgFirSpotAllAP(cc).AvgFirSpot=nanmean(ConFirSpotAllAP);  %Avg duration of all embryos of a construct
        AvgFirSpotAllAP(cc).ConSD=nanmean(ConFirSpotSD);
        AvgFirSpotAllAP(cc).ConSE=nanmean(ConFirSpotSE);
        AvgFirSpotAllAP(cc).AllFirSpots=[ConFirSpotAllAP];
        AvgFirSpotAllAP(cc).AllSD=nanstd(AvgFirSpotAllAP(cc).AllFirSpots);
        %AvgFirSpotAllAP(cc).AllSE=((AvgFirSpotAllAP(cc).AllSD)/sqrt(length(AvgFirSpotAllAP(cc).AllFirSpots)));
        for aa=1:length(APbinID)
        AvgFirSpotAllAP(cc).AllSE(aa)=(AvgFirSpotAllAP(cc).AllSD(aa))/sqrt(length(AvgFirSpotAllAP(cc).AllFirSpots(~isnan(AvgFirSpotAllAP(cc).AllFirSpots(:,aa)))));
        end
        AvgFirSpotAllAP(cc).All95Conf=(AvgFirSpotAllAP(cc).AllSE).*1.95;

end

%ignore AP bins with only one value 
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if length(AvgFirSpotAllAP(cc).AllFirSpots(~isnan(AvgFirSpotAllAP(cc).AllFirSpots(:,aa)))) ==1
            AvgFirSpotAllAP(cc).All95Conf(aa)=nan;
            AvgFirSpotAllAP(cc).AllSD(aa)=nan;
            AvgFirSpotAllAP(cc).AllSE(aa)=nan;
            AvgFirSpotAllAP(cc).AvgFirSpot(aa)=nan;
        end
    end
end
end
%%
if Halves ~= 'y'
%Plot mean of each embryo at specific AP position for each construct
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
Egglength=APbinID .*100;
for cc=1:length(ConstructList)
    for ee=1:length(AvgFirSpotAllAP(cc).EmbryosFirSpot)
        APBoxing(ee,cc)=AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse);
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


Colors(1).Color=DistalColor; Colors(2).Color=ProxColor; Colors(3).Color=BothSepColor;
Colors(4).Color=DistalEmptyColor;
Colors(5).Color=ProxEmptyColor;
Colors(6).Color=DoubDistColor; Colors(7).Color=DoubProxColor; Colors(8).Color=BothColor;
Colors(9).Color=BothEmptyColor;
Colors(10).Color=DistalColor;
Colors(11).Color=ProxColor;
Colors(12).Color=BothSepColor;
Colors(13).Color=BothColor;
Colors(14).Color=DoubProxColor;
Colors(15).Color=DoubProxColor;
Colors(16).Color=DistalColor;
Colors(17).Color=BothColor;

fontsize=18;
fontname='Helvetica';
FigDirect=[DropboxFolder filesep 'Figures' filesep 'Transcriptional dynamics' filesep 'Firstfluo'];
%%

CheckEmbryos=input('Check individual embryos?','s');
if CheckEmbryos=='y'
    figure 
    for cc=1:length(ConstructList)
        figure
        ConMean=[];
        EmbNumbs=[];
        for ee=1:length(AvgFirSpotAllAP(cc).EmbryosFirSpot)
            plot(AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot)
            hold on 
            ConMean=[ConMean; AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot];
            EmbNumbs{ee}=num2str(ee);
        end
        plot(nanmean(ConMean),'LineWidth',2,'Color','r');
        ConMeaned=nanmean(ConMean);
        legend(EmbNumbs);
        title(ConstructList{cc});
        Weirds=[];
        Counter=0;
        figure
        for ee=1:length(AvgFirSpotAllAP(cc).EmbryosFirSpot)
            Strange='n'
            for aa=1:length(APbinID)
            if abs(ConMeaned(aa)-(AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot(aa)))>=(0.45*ConMeaned(aa))
                Strange='y'
                break
            else
            end
            end
            if Strange=='y'
                plot(AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot)
                hold on 
                plot(ConMeaned,'LineWidth',2,'Color','r');
                Counter=Counter+1;
                Weirds{Counter}=num2str(ee);
            end
            end
            title([ConstructList{cc},'weirds']);
            legend(Weirds);
        end
    end

    for ee=1:length(AvgFirSpotAllAP(1).EmbryosFirSpot)
    plot(1,AvgFirSpotAllAP(1).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse),'o','LineWidth',1.5,'Color',DistalColor)
    hold on
    end
    for ee=1:length(AvgFirSpotAllAP(2).EmbryosFirSpot)
    plot(2,AvgFirSpotAllAP(2).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse),'o','LineWidth',1.5,'Color',ProxColor)
    hold on
    end
    for ee=1:length(AvgFirSpotAllAP(3).EmbryosFirSpot)
    plot(3,AvgFirSpotAllAP(3).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse),'o','LineWidth',1.5,'Color',BothSepColor)
    hold on
    end
    for ee=1:length(AvgFirSpotAllAP(4).EmbryosFirSpot)
    plot(4,AvgFirSpotAllAP(4).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse),'o','LineWidth',1.5,'Color',Colors(4).Color)
    hold on
    end
    for ee=1:length(AvgFirSpotAllAP(5).EmbryosFirSpot)
    plot(5,AvgFirSpotAllAP(5).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse),'o','LineWidth',1.5,'Color',Colors(5).Color)
    hold on
    end
    for ee=1:length(AvgFirSpotAllAP(6).EmbryosFirSpot)
    plot(6,AvgFirSpotAllAP(6).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse),'o','LineWidth',1.5,'Color',Colors(6).Color)
    hold on
    end
     for ee=1:length(AvgFirSpotAllAP(7).EmbryosFirSpot)
    plot(7,AvgFirSpotAllAP(7).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse),'o','LineWidth',1.5,'Color',Colors(7).Color)
    hold on
     end
    for ee=1:length(AvgFirSpotAllAP(8).EmbryosFirSpot)
    plot(8,AvgFirSpotAllAP(8).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse),'o','LineWidth',1.5,'Color',Colors(8).Color)
    hold on
    end
    for ee=1:length(AvgFirSpotAllAP(9).EmbryosFirSpot)
    plot(9,AvgFirSpotAllAP(9).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse),'o','LineWidth',1.5,'Color',Colors(9).Color)
    hold on
    end
    boxplot(APBoxing,'Colors','k');
xlim([0 10]);
xlabel('Construct');
xticks([1:9]);
xticklabels({'Dist', 'Prox', 'Both Sep', '1x Dist', '1x Prox', '2x Dist', '2x Prox', 'Both','1x Both'});
ylabel('Time into nc14 (min)');
title(['Mean time of first spot',' ' ,num2str(EgglengthUse),'% egg length'])

% Plot mean turn on time each construct vs AP position
figure 
for cc=1:length(ConstructList)
    plot(Egglength,AvgFirSpotAllAP(cc).AvgFirSpot,'Color',Colors(cc).Color,'LineWidth',1.5)
    hold on 
end
legend('Distal', 'Proximal', 'Both Sep','1x Distal','1x Proximal', '2x Distal', '2x Proximal', 'Both','1x Both','Dist32C','BothSep32C');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);
%%
% singles vs both sep
figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',1.5);
errorbar(Egglength, AvgFirSpotAllAP(3).AvgFirSpot, AvgFirSpotAllAP(3).All95Conf, 'Color',Colors(3).Color,'LineWidth',1.5);
legend('Distal', 'Proximal', 'Both Sep');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',2.5);
%legend('Distal', 'Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'SinglesFirstFluo'],'-dsvg');

figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',2.5);
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth',2.5);
%legend('Distal', 'Proximal','Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'SinglevBothFirstFluo'],'-dsvg');
%% single allele
figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(4).AvgFirSpot, AvgFirSpotAllAP(4).All95Conf, 'Color', Colors(1).Color, 'LineWidth',2.5,'LineStyle',':');
%legend('Distal', '1x Distal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'HemiDistFirstFluo'],'-dsvg');

figure
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color',Colors(2).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(5).AvgFirSpot, AvgFirSpotAllAP(5).All95Conf, 'Color', Colors(2).Color, 'LineWidth',2.5,'LineStyle',':');
%legend('Proximal', '1x Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'HemiProxFirstFluo'],'-dsvg');

% figure
% errorbar(Egglength, AvgFirSpotAllAP(4).AvgFirSpot, AvgFirSpotAllAP(4).All95Conf, 'Color',Colors(4).Color,'LineWidth',2.5);
% hold on 
% errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth',2.5);
% legend('1x Distal', 'Both');
% set(gca, 'FontSize', fontsize, 'FontName', fontname);
% ylabel('Time of first spot (min)');
% %ylabel('Time into nc14 (min)');
% xlabel('% egg length')
% %title('Mean time of first spot')
% xlim([0 100]);
% 
% figure
% errorbar(Egglength, AvgFirSpotAllAP(5).AvgFirSpot, AvgFirSpotAllAP(5).All95Conf, 'Color',Colors(5).Color,'LineWidth',2.5);
% hold on 
% errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth',2.5);
% legend('1x Proximal', 'Both');
% set(gca, 'FontSize', fontsize, 'FontName', fontname);
% ylabel('time of first spot (min)');
% %ylabel('Time into nc14 (min)');
% xlabel('% egg length')
% %title('Mean time of first spot')
% xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color',Colors(8).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(9).AvgFirSpot, AvgFirSpotAllAP(9).All95Conf, 'Color', Colors(8).Color, 'LineWidth',2.5,'LineStyle',':');
%legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'HemiBothFirstFluo'],'-dsvg');

figure
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).All95Conf, 'Color',Colors(7).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(15).AvgFirSpot, AvgFirSpotAllAP(15).All95Conf, 'Color', Colors(15).Color, 'LineWidth',2.5,'LineStyle',':');
%legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'Hemi2xProxFirstFluo'],'-dsvg');
%% temperature comparisions 
figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(10).AvgFirSpot, AvgFirSpotAllAP(10).All95Conf, 'Color', Colors(10).Color, 'LineWidth',2.5,'LineStyle','-.');
%legend('Distal', 'Distal32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);
print( [FigDirect filesep 'DistTCompFirstFluo'],'-dsvg');

figure
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color',Colors(2).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(11).AvgFirSpot, AvgFirSpotAllAP(11).All95Conf, 'Color', Colors(11).Color, 'LineWidth',2.5,'LineStyle','-.');
%legend('Distal', 'Distal32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);
print( [FigDirect filesep 'ProxTCompFirstFluo'],'-dsvg');

figure
errorbar(Egglength, AvgFirSpotAllAP(3).AvgFirSpot, AvgFirSpotAllAP(3).All95Conf, 'Color',Colors(3).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(12).AvgFirSpot, AvgFirSpotAllAP(12).All95Conf, 'Color', Colors(12).Color, 'LineWidth',2.5,'LineStyle','-.');
%legend('Both Sep', 'Both Sep32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color',Colors(8).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(13).AvgFirSpot, AvgFirSpotAllAP(13).All95Conf, 'Color', Colors(13).Color, 'LineWidth',2.5,'LineStyle','-.');
%legend('Both', 'Both 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);
print( [FigDirect filesep 'BothTCompFirstFluo'],'-dsvg');

figure
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).All95Conf, 'Color',Colors(7).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(14).AvgFirSpot, AvgFirSpotAllAP(14).All95Conf, 'Color', Colors(14).Color, 'LineWidth',2.5,'LineStyle','-.');
%legend('2x Proximal', '2x Proximal 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);
print( [FigDirect filesep '2xProxTCompFirstFluo'],'-dsvg');

%% 17C vs RT
figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(16).AvgFirSpot, AvgFirSpotAllAP(16).All95Conf, 'Color', Colors(16).Color, 'LineWidth',2.5,'LineStyle','-.');
%legend('Distal', 'Distal32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);
print( [FigDirect filesep 'Dist17CCompFirstFluo'],'-dsvg');

figure
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color',Colors(8).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(17).AvgFirSpot, AvgFirSpotAllAP(17).All95Conf, 'Color', Colors(17).Color, 'LineWidth',2.5,'LineStyle','-.');
%legend('Distal', 'Distal32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);
print( [FigDirect filesep 'Both17CCompFirstFluo'],'-dsvg');

%%
%doubles vs SE
figure
errorbar(Egglength, AvgFirSpotAllAP(6).AvgFirSpot, AvgFirSpotAllAP(6).All95Conf, 'Color',Colors(6).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth',1.5);
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color',Colors(8).Color,'LineWidth',1.5);
legend('2x Distal', '2x Proximal', 'Both');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
title('Mean time of first spot')
xlim([0 100]);

figure 
errorbar(Egglength, AvgFirSpotAllAP(6).AvgFirSpot, AvgFirSpotAllAP(6).ConSE, 'Color',Colors(6).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).ConSE, 'Color',Colors(8).Color,'LineWidth',1.5);
legend('2x Distal',  'Both');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).ConSE, 'Color',Colors(7).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).ConSE, 'Color',Colors(8).Color,'LineWidth',1.5);
legend('2x Proximal',  'Both');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

% singles vs doubles
figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).ConSE, 'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(6).AvgFirSpot, AvgFirSpotAllAP(6).ConSE, 'Color',Colors(6).Color,'LineWidth',1.5);
legend('Distal',  '2x Distal');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).ConSE, 'Color',Colors(2).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).ConSE, 'Color',Colors(7).Color,'LineWidth',1.5);
legend('Proximal',  '2x Proximal');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

%%
%ANOVA at specific AP position of constructs against one another

for cc=1:length(ConstructList)
    for bb=1:length([AvgFirSpotAllAP(cc).AllFirSpots])
    ConComp(bb,cc)=AvgFirSpotAllAP(cc).AllFirSpots(bb,APtoUse);
    end
    ConComp(ConComp==0)=nan;
end
[p,tbl,stats]=anova1(ConComp);
xlabel('Construct')
xticks([1:8]);
xticklabels({'Dist','Prox','BothSep','1x Dist', '1x Prox','2xDist','2xProx','Both'});
ylabel('Time into nc14 (min)');
title(['Avg time of first fluorescence',' ',num2str(APtoUse)]);


%save info 
save([DropboxFolder filesep 'Constructs' filesep 'FirstTimeON'],'AvgFirSpotAllAP');
end
%%
if Halves =='y'
    AvgFirstSpot=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgFirstSpotCon=[];
    firsttime=1;
    ConFirSpotAllAP=[];
    ConFirSpotSE=[];
    ConFirSpotSD=[];
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
        FirSpotAllAP=[];
        for aa=1:length(APbinID)
            FirSpotAP=[];
            FirstSpotAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(FirstSpotAP)
                FirSpotAP=[FirSpotAP; nan];
            else
            for bb=1:length(FirstSpotAP)
                if ~isempty(BurstProperties(FirstSpotAP(bb)).Duration)
                FirSpotAP=[FirSpotAP;[BurstProperties(FirstSpotAP(bb)).FirstTimeOn]'];  %put all durations at a given AP value in a column going down
                else
                    FirSpotAP=[FirSpotAP; 0];
                end
               
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(FirSpotAP)
                FirSpotAllAP(bb,aa,ee)=FirSpotAP(bb);
            end
            FirSpotAllAP(FirSpotAllAP==0)=nan;
            
            FirstSpotSD(ee,aa,cc)=nanstd(FirSpotAllAP(:,aa,ee));
        FirstSpotSE(ee,aa,cc)=FirstSpotSD(ee,aa,cc)/sqrt(length(FirstSpotAP));             
            clear FirstSpotAP
        end
        AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot=nanmean(FirSpotAllAP(:,:,ee));
        AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).SD=FirstSpotSD(ee,:,cc);
        AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).SE=FirstSpotSE(ee,:,cc);
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(FirSpotAllAP,3)
            ConFirSpotAllAP=[ConFirSpotAllAP; FirSpotAllAP(:,:,bb)];
        end
        for bb=1:size(FirstSpotSD,3)
            ConFirSpotSD=[ConFirSpotSD; FirstSpotSD(:,:,bb)];
        end
        for bb=1:size(FirstSpotSE,3)
            ConFirSpotSE=[ConFirSpotSE;FirstSpotSE(:,:,bb)];
        end
        
    end
        AvgFirSpotAllAP(cc).AvgFirSpot=nanmean(ConFirSpotAllAP);  %Avg duration of all embryos of a construct
        AvgFirSpotAllAP(cc).ConSD=nanmean(ConFirSpotSD);
        AvgFirSpotAllAP(cc).ConSE=nanmean(ConFirSpotSE);
        AvgFirSpotAllAP(cc).AllFirSpots=[ConFirSpotAllAP];
        AvgFirSpotAllAP(cc).AllSD=nanstd(AvgFirSpotAllAP(cc).AllFirSpots);
        AvgFirSpotAllAP(cc).AllSE=((AvgFirSpotAllAP(cc).AllSD)/sqrt(length(AvgFirSpotAllAP(cc).AllFirSpots)));
        AvgFirSpotAllAP(cc).All95Conf=(AvgFirSpotAllAP(cc).AllSE).*1.95;

end

%% plot 1st half 
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
Egglength=APbinID .*100;
for cc=1:length(ConstructList)
    for ee=1:2:length(AvgFirSpotAllAP(cc).EmbryosFirSpot)
        APBoxing(ee,cc)=AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse);
    end
end
APBoxing(APBoxing==0)=nan;
figure
DistalColor=[8 180 238] ./ 255;
DistalEmptyColor=[8 210 238] ./ 255; 
Distal32CColor=[118 180 238] ./ 255;
DoubDistColor=[1 17 181] ./ 255;
ProxColor=[251 230 60] ./ 255;
ProxEmptyColor=[251 250 50] ./255;
DoubProxColor=[251 190 80] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[12 195 82] ./ 255;
BothEmptyColor=[12 250 150] ./ 255;
Both32CColor=[120 195 82] ./ 255;
DoubProx32CColor=[200 150 100] ./ 255;


Colors(1).Color=DistalColor; Colors(2).Color=ProxColor; Colors(3).Color=BothSepColor;
Colors(4).Color=DistalEmptyColor;
Colors(5).Color=ProxEmptyColor;
Colors(6).Color=DoubDistColor; Colors(7).Color=DoubProxColor; Colors(8).Color=BothColor;
Colors(9).Color=BothEmptyColor;
Colors(10).Color=Distal32CColor;
Colors(11).Color=BothSep32CColor;
Colors(12).Color=Both32CColor;
Colors(13).Color=DoubProx32CColor;

fontsize=18;
fontname='Helvetica';
for cc=1:length(ConstructList)
    for ee=1:2:length(AvgFirSpotAllAP(cc).EmbryosFirSpot)
        plot(cc,AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse),'o','LineWidth',1.5,'Color',Colors(cc).Color)
    hold on 
    end
end
    boxplot(APBoxing,'Colors','k');
xlim([0 10]);
xlabel('Construct');
xticks([1:9]);
xticklabels({'Dist', 'Prox', 'Both Sep', '1x Dist', '1x Prox', '2x Dist', '2x Prox', 'Both','1x Both'});
ylabel('Time into nc14 (min)');
title(['Mean time of first spot',' ' ,num2str(EgglengthUse),'% egg length'])

% Plot mean turn on time each construct vs AP position
figure 
for cc=1:length(ConstructList)
    plot(Egglength,AvgFirSpotAllAP(cc).AvgFirSpot,'Color',Colors(cc).Color,'LineWidth',1.5)
    hold on 
end
legend('Distal', 'Proximal', 'Both Sep','1x Distal','1x Proximal', '2x Distal', '2x Proximal', 'Both','1x Both','Dist32C','BothSep32C');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

% singles vs both sep
figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',1.5);
errorbar(Egglength, AvgFirSpotAllAP(3).AvgFirSpot, AvgFirSpotAllAP(3).All95Conf, 'Color',Colors(3).Color,'LineWidth',1.5);
legend('Distal', 'Proximal', 'Both Sep');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',2.5);
%legend('Distal', 'Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',2.5);
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth',2.5);
legend('Distal', 'Proximal','Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(4).AvgFirSpot, AvgFirSpotAllAP(4).All95Conf, 'Color', Colors(4).Color, 'LineWidth',2.5);
legend('Distal', '1x Distal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color',Colors(2).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(5).AvgFirSpot, AvgFirSpotAllAP(5).All95Conf, 'Color', Colors(5).Color, 'LineWidth',2.5);
legend('Proximal', '1x Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(4).AvgFirSpot, AvgFirSpotAllAP(4).All95Conf, 'Color',Colors(4).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth',2.5);
legend('1x Distal', 'Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(5).AvgFirSpot, AvgFirSpotAllAP(5).All95Conf, 'Color',Colors(5).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth',2.5);
legend('1x Proximal', 'Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color',Colors(8).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(9).AvgFirSpot, AvgFirSpotAllAP(9).All95Conf, 'Color', Colors(9).Color, 'LineWidth',2.5);
legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(10).AvgFirSpot, AvgFirSpotAllAP(10).All95Conf, 'Color', Colors(10).Color, 'LineWidth',2.5);
legend('Distal', 'Distal32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(3).AvgFirSpot, AvgFirSpotAllAP(3).All95Conf, 'Color',Colors(3).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(11).AvgFirSpot, AvgFirSpotAllAP(11).All95Conf, 'Color', Colors(11).Color, 'LineWidth',2.5);
legend('Both Sep', 'Both Sep32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color',Colors(8).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(12).AvgFirSpot, AvgFirSpotAllAP(12).All95Conf, 'Color', Colors(12).Color, 'LineWidth',2.5);
legend('Both', 'Both 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).All95Conf, 'Color',Colors(7).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(13).AvgFirSpot, AvgFirSpotAllAP(13).All95Conf, 'Color', Colors(13).Color, 'LineWidth',2.5);
legend('2x Proximal', '2x Proximal 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);

%doubles vs SE
figure
errorbar(Egglength, AvgFirSpotAllAP(6).AvgFirSpot, AvgFirSpotAllAP(6).All95Conf, 'Color',Colors(6).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth',1.5);
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color',Colors(8).Color,'LineWidth',1.5);
legend('2x Distal', '2x Proximal', 'Both');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
title('Mean time of first spot')
xlim([0 100]);

figure 
errorbar(Egglength, AvgFirSpotAllAP(6).AvgFirSpot, AvgFirSpotAllAP(6).ConSE, 'Color',Colors(6).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).ConSE, 'Color',Colors(8).Color,'LineWidth',1.5);
legend('2x Distal',  'Both');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).ConSE, 'Color',Colors(7).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).ConSE, 'Color',Colors(8).Color,'LineWidth',1.5);
legend('2x Proximal',  'Both');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

% singles vs doubles
figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).ConSE, 'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(6).AvgFirSpot, AvgFirSpotAllAP(6).ConSE, 'Color',Colors(6).Color,'LineWidth',1.5);
legend('Distal',  '2x Distal');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).ConSE, 'Color',Colors(2).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).ConSE, 'Color',Colors(7).Color,'LineWidth',1.5);
legend('Proximal',  '2x Proximal');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);


%ANOVA at specific AP position of constructs against one another

for cc=1:length(ConstructList)
    for bb=1:length([AvgFirSpotAllAP(cc).AllFirSpots])
    ConComp(bb,cc)=AvgFirSpotAllAP(cc).AllFirSpots(bb,APtoUse);
    end
    ConComp(ConComp==0)=nan;
end
[p,tbl,stats]=anova1(ConComp);
xlabel('Construct')
xticks([1:8]);
xticklabels({'Dist','Prox','BothSep','1x Dist', '1x Prox','2xDist','2xProx','Both'});
ylabel('Time into nc14 (min)');
title(['Avg time of first fluorescence',' ',num2str(APtoUse)]);

%% run 2nd half 
AvgFirstSpot=[];
clear AvgFirSpotAllAP
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgFirstSpotCon=[];
    firsttime=1;
    ConFirSpotAllAP=[];
    ConFirSpotSE=[];
    ConFirSpotSD=[];
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
        FirSpotAllAP=[];
        for aa=1:length(APbinID)
            FirSpotAP=[];
            FirstSpotAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(FirstSpotAP)
                FirSpotAP=[FirSpotAP; nan];
            else
            for bb=1:length(FirstSpotAP)
                if ~isempty(BurstProperties(FirstSpotAP(bb)).Duration)
                FirSpotAP=[FirSpotAP;[BurstProperties(FirstSpotAP(bb)).FirstTimeOn]'];  %put all durations at a given AP value in a column going down
                else
                    FirSpotAP=[FirSpotAP; 0];
                end
               
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(FirSpotAP)
                FirSpotAllAP(bb,aa,ee)=FirSpotAP(bb);
            end
            FirSpotAllAP(FirSpotAllAP==0)=nan;
            
            FirstSpotSD(ee,aa,cc)=nanstd(FirSpotAllAP(:,aa,ee));
        FirstSpotSE(ee,aa,cc)=FirstSpotSD(ee,aa,cc)/sqrt(length(FirstSpotAP));             
            clear FirstSpotAP
        end
        AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot=nanmean(FirSpotAllAP(:,:,ee));
        AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).SD=FirstSpotSD(ee,:,cc);
        AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).SE=FirstSpotSE(ee,:,cc);
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(FirSpotAllAP,3)
            ConFirSpotAllAP=[ConFirSpotAllAP; FirSpotAllAP(:,:,bb)];
        end
        for bb=1:size(FirstSpotSD,3)
            ConFirSpotSD=[ConFirSpotSD; FirstSpotSD(:,:,bb)];
        end
        for bb=1:size(FirstSpotSE,3)
            ConFirSpotSE=[ConFirSpotSE;FirstSpotSE(:,:,bb)];
        end
        
    end
        AvgFirSpotAllAP(cc).AvgFirSpot=nanmean(ConFirSpotAllAP);  %Avg duration of all embryos of a construct
        AvgFirSpotAllAP(cc).ConSD=nanmean(ConFirSpotSD);
        AvgFirSpotAllAP(cc).ConSE=nanmean(ConFirSpotSE);
        AvgFirSpotAllAP(cc).AllFirSpots=[ConFirSpotAllAP];
        AvgFirSpotAllAP(cc).AllSD=nanstd(AvgFirSpotAllAP(cc).AllFirSpots);
        AvgFirSpotAllAP(cc).AllSE=((AvgFirSpotAllAP(cc).AllSD)/sqrt(length(AvgFirSpotAllAP(cc).AllFirSpots)));
        AvgFirSpotAllAP(cc).All95Conf=(AvgFirSpotAllAP(cc).AllSE).*1.95;

end

%% Plot 2nd half
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
Egglength=APbinID .*100;
for cc=1:length(ConstructList)
    for ee=2:2:length(AvgFirSpotAllAP(cc).EmbryosFirSpot)
        APBoxing(ee,cc)=AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse);
    end
end
APBoxing(APBoxing==0)=nan;
figure
DistalColor=[8 180 238] ./ 255;
DistalEmptyColor=[8 210 238] ./ 255; 
Distal32CColor=[118 180 238] ./ 255;
DoubDistColor=[1 17 181] ./ 255;
ProxColor=[251 230 60] ./ 255;
ProxEmptyColor=[251 250 50] ./255;
DoubProxColor=[251 190 80] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[12 195 82] ./ 255;
BothEmptyColor=[12 250 150] ./ 255;
Both32CColor=[120 195 82] ./ 255;
DoubProx32CColor=[200 150 100] ./ 255;


Colors(1).Color=DistalColor; Colors(2).Color=ProxColor; Colors(3).Color=BothSepColor;
Colors(4).Color=DistalEmptyColor;
Colors(5).Color=ProxEmptyColor;
Colors(6).Color=DoubDistColor; Colors(7).Color=DoubProxColor; Colors(8).Color=BothColor;
Colors(9).Color=BothEmptyColor;
Colors(10).Color=Distal32CColor;
Colors(11).Color=BothSep32CColor;
Colors(12).Color=Both32CColor;
Colors(13).Color=DoubProx32CColor;

fontsize=18;
fontname='Helvetica';
for cc=1:length(ConstructList)
    for ee=2:2:length(AvgFirSpotAllAP(cc).EmbryosFirSpot)
        %if ~isempty(AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot)
        plot(cc,AvgFirSpotAllAP(cc).EmbryosFirSpot(ee).MeanFirstSpot(APtoUse),'o','LineWidth',1.5,'Color',Colors(cc).Color)
        %end  
        hold on 
    end
end
    boxplot(APBoxing,'Colors','k');
xlim([0 10]);
xlabel('Construct');
xticks([1:9]);
xticklabels({'Dist', 'Prox', 'Both Sep', '1x Dist', '1x Prox', '2x Dist', '2x Prox', 'Both','1x Both'});
ylabel('Time into nc14 (min)');
title(['Mean time of first spot',' ' ,num2str(EgglengthUse),'% egg length'])

% Plot mean turn on time each construct vs AP position
figure 
for cc=1:length(ConstructList)
    plot(Egglength,AvgFirSpotAllAP(cc).AvgFirSpot,'Color',Colors(cc).Color,'LineWidth',1.5)
    hold on 
end
legend('Distal', 'Proximal', 'Both Sep','1x Distal','1x Proximal', '2x Distal', '2x Proximal', 'Both','1x Both','Dist32C','BothSep32C');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

% singles vs both sep
figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',1.5);
errorbar(Egglength, AvgFirSpotAllAP(3).AvgFirSpot, AvgFirSpotAllAP(3).All95Conf, 'Color',Colors(3).Color,'LineWidth',1.5);
legend('Distal', 'Proximal', 'Both Sep');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',2.5);
%legend('Distal', 'Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',2.5);
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth',2.5);
legend('Distal', 'Proximal','Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(4).AvgFirSpot, AvgFirSpotAllAP(4).All95Conf, 'Color', Colors(4).Color, 'LineWidth',2.5);
legend('Distal', '1x Distal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).All95Conf, 'Color',Colors(2).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(5).AvgFirSpot, AvgFirSpotAllAP(5).All95Conf, 'Color', Colors(5).Color, 'LineWidth',2.5);
legend('Proximal', '1x Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(4).AvgFirSpot, AvgFirSpotAllAP(4).All95Conf, 'Color',Colors(4).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth',2.5);
legend('1x Distal', 'Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(5).AvgFirSpot, AvgFirSpotAllAP(5).All95Conf, 'Color',Colors(5).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color', Colors(8).Color, 'LineWidth',2.5);
legend('1x Proximal', 'Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color',Colors(8).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(9).AvgFirSpot, AvgFirSpotAllAP(9).All95Conf, 'Color', Colors(9).Color, 'LineWidth',2.5);
legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
%ylabel('Time into nc14 (min)');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(10).AvgFirSpot, AvgFirSpotAllAP(10).All95Conf, 'Color', Colors(10).Color, 'LineWidth',2.5);
legend('Distal', 'Distal32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(3).AvgFirSpot, AvgFirSpotAllAP(3).All95Conf, 'Color',Colors(3).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(11).AvgFirSpot, AvgFirSpotAllAP(11).All95Conf, 'Color', Colors(11).Color, 'LineWidth',2.5);
legend('Both Sep', 'Both Sep32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color',Colors(8).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(12).AvgFirSpot, AvgFirSpotAllAP(12).All95Conf, 'Color', Colors(12).Color, 'LineWidth',2.5);
legend('Both', 'Both 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).All95Conf, 'Color',Colors(7).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(13).AvgFirSpot, AvgFirSpotAllAP(13).All95Conf, 'Color', Colors(13).Color, 'LineWidth',2.5);
legend('2x Proximal', '2x Proximal 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Time of first spot (min)');
xlabel('% egg length')
xlim([0 100]);

%doubles vs SE
figure
errorbar(Egglength, AvgFirSpotAllAP(6).AvgFirSpot, AvgFirSpotAllAP(6).All95Conf, 'Color',Colors(6).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).All95Conf, 'Color', Colors(7).Color, 'LineWidth',1.5);
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).All95Conf, 'Color',Colors(8).Color,'LineWidth',1.5);
legend('2x Distal', '2x Proximal', 'Both');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
title('Mean time of first spot')
xlim([0 100]);

figure 
errorbar(Egglength, AvgFirSpotAllAP(6).AvgFirSpot, AvgFirSpotAllAP(6).ConSE, 'Color',Colors(6).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).ConSE, 'Color',Colors(8).Color,'LineWidth',1.5);
legend('2x Distal',  'Both');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).ConSE, 'Color',Colors(7).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(8).AvgFirSpot, AvgFirSpotAllAP(8).ConSE, 'Color',Colors(8).Color,'LineWidth',1.5);
legend('2x Proximal',  'Both');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

% singles vs doubles
figure
errorbar(Egglength, AvgFirSpotAllAP(1).AvgFirSpot, AvgFirSpotAllAP(1).ConSE, 'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(6).AvgFirSpot, AvgFirSpotAllAP(6).ConSE, 'Color',Colors(6).Color,'LineWidth',1.5);
legend('Distal',  '2x Distal');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);

figure
errorbar(Egglength, AvgFirSpotAllAP(2).AvgFirSpot, AvgFirSpotAllAP(2).ConSE, 'Color',Colors(2).Color,'LineWidth',1.5);
hold on 
errorbar(Egglength, AvgFirSpotAllAP(7).AvgFirSpot, AvgFirSpotAllAP(7).ConSE, 'Color',Colors(7).Color,'LineWidth',1.5);
legend('Proximal',  '2x Proximal');
ylabel('Time into nc14 (min)');
xlabel('% egg length')
xlim([0 100]);


%ANOVA at specific AP position of constructs against one another

for cc=1:length(ConstructList)
    for bb=1:length([AvgFirSpotAllAP(cc).AllFirSpots])
    ConComp(bb,cc)=AvgFirSpotAllAP(cc).AllFirSpots(bb,APtoUse);
    end
    ConComp(ConComp==0)=nan;
end
[p,tbl,stats]=anova1(ConComp);
xlabel('Construct')
xticks([1:8]);
xticklabels({'Dist','Prox','BothSep','1x Dist', '1x Prox','2xDist','2xProx','Both'});
ylabel('Time into nc14 (min)');
title(['Avg time of first fluorescence',' ',num2str(APtoUse)]);
end
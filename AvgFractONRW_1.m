%% Fraction of nc14 constructs are ON 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'Kr2xProxEmpty';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'HbEmpty'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
SlopeUse=input('Want to use slope calculations?','s');
%Count time of fluorescence turns on for each construct
AvgFractON=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgFractONCon=[];
    firsttime=1;
    ConFractONAllAP=[];
    ConFractONSE=[];
    ConFractONSD=[];
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
        FractONAllAP=[];
        for aa=1:length(APbinID)
            FractONAP=[];
            FractionONAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(FractionONAP)
                FractONAP=[FractONAP; nan];
            else
            for bb=1:length(FractionONAP)
                if ~isempty(BurstProperties(FractionONAP(bb)).Duration)
                FractONAP=[FractONAP;[BurstProperties(FractionONAP(bb)).FractON]'];  %put all durations at a given AP value in a column going down
                else
                    FractONAP=[FractONAP; nan];
                end
               
                end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(FractONAP)
                FractONAllAP(bb,aa,ee)=FractONAP(bb);
            end
            FractONAllAP(FractONAllAP==0)=nan;
            
            FractionONSD(ee,aa,cc)=nanstd(FractONAllAP(:,aa,ee));
        %FractionONSE(ee,aa,cc)=FractionONSD(ee,aa,cc)/sqrt(length(FractionONAP));
        FractionONSE(ee,aa,cc)=FractionONSD(ee,aa,cc)/sqrt(sum(~isnan(FractONAP)));
            clear FractionONAP
        end
        AvgFractONAllAP(cc).EmbryosFirSpot(ee).MeanFractON=nanmean(FractONAllAP(:,:,ee));
        AvgFractONAllAP(cc).EmbryosFirSpot(ee).SD=FractionONSD(ee,:,cc);
        AvgFractONAllAP(cc).EmbryosFirSpot(ee).SE=FractionONSE(ee,:,cc);
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(FractONAllAP,3)
            ConFractONAllAP=[ConFractONAllAP; FractONAllAP(:,:,bb)];
        end
        for bb=1:size(FractionONSD,3)
            ConFractONSD=[ConFractONSD; FractionONSD(:,:,bb)];
        end
        for bb=1:size(FractionONSE,3)
            ConFractONSE=[ConFractONSE;FractionONSE(:,:,bb)];
        end
        
    end
        AvgFractONAllAP(cc).AvgFractON=nanmean(ConFractONAllAP);  %Avg duration of all embryos of a construct
        AvgFractONAllAP(cc).ConSD=nanmean(ConFractONSD);
        AvgFractONAllAP(cc).ConSE=nanmean(ConFractONSE);
        AvgFractONAllAP(cc).AllFractON=[ConFractONAllAP];
        AvgFractONAllAP(cc).AllSD=nanstd(AvgFractONAllAP(cc).AllFractON);
        %AvgFractONAllAP(cc).AllSE=((AvgFractONAllAP(cc).AllSD)/sqrt(length(AvgFractONAllAP(cc).AllFractON)));
        for aa=1:length(APbinID)
            AvgFractONAllAP(cc).AllSE(aa)=(AvgFractONAllAP(cc).AllSD(aa))/sqrt(sum(~isnan(AvgFractONAllAP(cc).AllFractON(:,aa))));
        end
        AvgFractONAllAP(cc).All95Conf=(AvgFractONAllAP(cc).AllSE).*1.95;

end

%% Plotting 
Egglength=APbinID .*100;

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
Colors(6).Color=DoubDistColor; Colors(7).Color=DoubProxColor; Colors(8).Color=DoubProxEmptyColor;
Colors(9).Color=BothColor;
Colors(10).Color=BothEmptyColor;
Colors(11).Color=DistalColor;
Colors(12).Color=ProxColor
Colors(13).Color=BothSepColor;
Colors(14).Color=BothColor;
Colors(15).Color=DoubProxColor;

fontsize=18;
fontname='Helvetica';
FigDirect=[DropboxFolder filesep 'Figures' filesep 'Transcriptional dynamics' filesep 'FractON'];

%%
% singles vs both sep
figure
errorbar(Egglength, AvgFractONAllAP(1).AvgFractON, AvgFractONAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(2).AvgFractON, AvgFractONAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',2.5);
errorbar(Egglength, AvgFractONAllAP(3).AvgFractON, AvgFractONAllAP(3).All95Conf, 'Color',Colors(3).Color,'LineWidth',2.5);
legend('Distal', 'Proximal', 'Both Sep');
ylabel('Fraction of nc14 spent ON');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);

figure
errorbar(Egglength, AvgFractONAllAP(1).AvgFractON, AvgFractONAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(2).AvgFractON, AvgFractONAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',2.5);
%legend('Distal', 'Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('fraction of nc14 spent ON');
xlabel('% egg length')
xlim([0 100]);
print( [FigDirect filesep 'SinglesFractON'],'-dsvg');


figure
errorbar(Egglength, AvgFractONAllAP(1).AvgFractON, AvgFractONAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(2).AvgFractON, AvgFractONAllAP(2).All95Conf, 'Color', Colors(2).Color, 'LineWidth',2.5);
errorbar(Egglength, AvgFractONAllAP(9).AvgFractON, AvgFractONAllAP(9).All95Conf, 'Color', Colors(9).Color, 'LineWidth',2.5);
%legend('Distal', 'Proximal','Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('fraction of nc14 spent on');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'SinglesvBothFractON'],'-dsvg');

%% Duplications vs SE
figure
errorbar(Egglength, AvgFractONAllAP(4).AvgFractON, AvgFractONAllAP(4).All95Conf, 'Color',Colors(4).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(5).AvgFractON, AvgFractONAllAP(5).All95Conf, 'Color', Colors(5).Color, 'LineWidth',2.5);
errorbar(Egglength, AvgFractONAllAP(9).AvgFractON, AvgFractONAllAP(9).All95Conf, 'Color', Colors(9).Color, 'LineWidth',2.5);
%legend('Distal', 'Proximal','Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('fraction of nc14 spent on');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
%%
figure
errorbar(Egglength, AvgFractONAllAP(1).AvgFractON, AvgFractONAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(4).AvgFractON, AvgFractONAllAP(4).All95Conf, 'Color', Colors(1).Color, 'LineWidth',2.5,'LineStyle',':');
%legend('Distal', '1x Distal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('fraction of nc14 spent ON');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'HemiDistFractON'],'-dsvg');

figure
errorbar(Egglength, AvgFractONAllAP(2).AvgFractON, AvgFractONAllAP(2).All95Conf, 'Color',Colors(2).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(5).AvgFractON, AvgFractONAllAP(5).All95Conf, 'Color', Colors(2).Color, 'LineWidth',2.5,'LineStyle',':');
%legend('Proximal', '1x Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('fraction of nc14 spent ON');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'HemiProxFractON'],'-dsvg');

figure
errorbar(Egglength, AvgFractONAllAP(9).AvgFractON, AvgFractONAllAP(9).All95Conf, 'Color',Colors(9).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(10).AvgFractON, AvgFractONAllAP(10).All95Conf, 'Color', Colors(9).Color, 'LineWidth',2.5,'LineStyle',':');
%legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('fraction of nc14 spent ON');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'HemiBothFractON'],'-dsvg');

figure
errorbar(Egglength, AvgFractONAllAP(7).AvgFractON, AvgFractONAllAP(7).All95Conf, 'Color',Colors(7).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(8).AvgFractON, AvgFractONAllAP(8).All95Conf, 'Color', Colors(7).Color, 'LineWidth',2.5,'LineStyle',':');
%legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('fraction of nc14 spent ON');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'Hemi2xProxFractON'],'-dsvg');
%%
%Altered temperature 
figure
errorbar(Egglength, AvgFractONAllAP(1).AvgFractON, AvgFractONAllAP(1).All95Conf, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(11).AvgFractON, AvgFractONAllAP(11).All95Conf, 'Color', Colors(11).Color, 'LineWidth',2.5,'LineStyle','-.');
%legend('Distal', 'Distal 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('fraction of nc14 spent ON');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'DistTCompFractON'],'-dsvg');


figure
errorbar(Egglength, AvgFractONAllAP(3).AvgFractON, AvgFractONAllAP(3).All95Conf, 'Color',Colors(3).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(13).AvgFractON, AvgFractONAllAP(13).All95Conf, 'Color', Colors(13).Color, 'LineWidth',2.5,'LineStyle','-.');
%legend('Both Sep', 'Both Sep 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('fraction of nc14 spent ON');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'ProxTCompFractON'],'-dsvg');


figure
errorbar(Egglength, AvgFractONAllAP(9).AvgFractON, AvgFractONAllAP(9).All95Conf, 'Color',Colors(9).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(14).AvgFractON, AvgFractONAllAP(14).All95Conf, 'Color', Colors(14).Color, 'LineWidth',2.5,'LineStyle','-.');
%legend('Both', 'Both 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('fraction of nc14 spent ON');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep 'BothTCompFractON'],'-dsvg');


figure
errorbar(Egglength, AvgFractONAllAP(7).AvgFractON, AvgFractONAllAP(7).All95Conf, 'Color',Colors(7).Color,'LineWidth',2.5);
hold on 
errorbar(Egglength, AvgFractONAllAP(15).AvgFractON, AvgFractONAllAP(15).All95Conf, 'Color', Colors(15).Color, 'LineWidth',2.5,'LineStyle','-.');
%legend('2x Proximal', '2x Proximal 32C','Location','best');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('fraction of nc14 spent ON');
xlabel('% egg length')
%title('Mean time of first spot')
xlim([0 100]);
print( [FigDirect filesep '2xProxTCompFractON'],'-dsvg');


%Save Fraction ON info
save([DropboxFolder filesep 'Constructs' filesep 'FractON'],'AvgFractONAllAP');
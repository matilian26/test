%% calculate coefficient of variance for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'KrDistEmpty';'KrProxEmpty';'KrBothEmpty';'Kr2xProxEmpty';'KrDist17C';'KrBoth17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info   for right now just going to use nc14
ncUse=input('Want to only use nc14?','s');
NumberNuclei=nan(20,41,[length(ConstructList)]);
% go through each embryo of each construct
for cc=1:length(ConstructList)
     Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    APTable=[];
    ConTimeAvg=[];
    AllNucTimeTable=[];
    AllTotmRNA=[];
    AllNucStrengthTable=[];
    ConFluoMeanz=[];
    for ee=1:NEmbryos
        
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        else
            Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelation.mat'];
        end
        load(Filename);
        %if ncUse=='y'
%         nc_number=[CompiledParticles.nc];
%         CompiledParticles_14=CompiledParticles(nc_number==14);
        APstuff=[SpotDiff.APBin];
        %Spotstuff=[SpotDiff.SpotOne];
        %TimeTable=zeros(1000,41);
        for aa=1:length(APbinID)
            APsubset=[];
            APsubset=SpotDiff(APstuff==APbinID(aa));
            Spotstuff=[APsubset.SpotOne];
            APSpotsubset=APsubset(~isempty(Spotstuff));
            APNucs(aa)=length(APsubset);
        end
        MostAPNucs=max(APNucs);
        MostAPNucs=2*MostAPNucs;
        TimeTable=nan(MostAPNucs,41);
        TotmRNATable=nan(MostAPNucs,41);
        TotmRNATable=nan(MostAPNucs,41);
        StrengthTable=nan(MostAPNucs, 41);
        MeanzTable=nan(MostAPNucs,41);
        for aa=1:length(APbinID)
            APsubset=[];
%             VarPart=[];
%             SumVarPart=[];
            APsubset=SpotDiff(APstuff==APbinID(aa));
            if ~isempty(APsubset)
               
            for bb=1:length(APsubset)
                NumberNuclei(ee,aa,cc)=length(APsubset);
                TimeLength=[];
                TimeLength=sum(~isnan(APsubset(bb).SpotOne));   %Only want to look at frames where nucleus exisits (i.e 0 or number values)
                if TimeLength==1
                    AvgFluo=nan;
                else
                AvgFluo=nanmean([APsubset(bb).SpotOne]);
                end
                VarFluo=nanstd([APsubset(bb).SpotOne]);
                VarianceFluo=var([APsubset(bb).SpotOne],'omitnan');
                
            TimeTable(bb,aa)=(VarFluo/AvgFluo);
            if isempty(APsubset(bb).TotalmRNAOne)
                TotmRNATable(bb,aa)=nan;
            else
            TotmRNATable(bb,aa)=APsubset(bb).TotalmRNAOne;
            end
            MeanzTable(bb,aa)=AvgFluo;
            StrengthTable(bb,aa)=(VarianceFluo/AvgFluo);
            %FluoValues(cc).Embryo(ee).Nucleus(bb).FluoValue=FluoVals;
            end
            for bb=1:length(APsubset)
                if isfield(APsubset,'SpotTwo') & (~isempty(APsubset(bb).SpotTwo))
                    TimeLength2=[];
                TimeLength2=(sum(~isnan(APsubset(bb).SpotTwo)));
                if TimeLength2==1
                    AvgFluo2=nan;
                else
                AvgFluo2=nanmean([APsubset(bb).SpotTwo]);
                end
                VarFluo2=nanstd([APsubset(bb).SpotTwo]);
                VarianceFluo2=var([APsubset(bb).SpotTwo],'omitnan');
                TimeTable(bb+(length(APsubset)),aa)=(VarFluo2/AvgFluo2);
                if isempty(APsubset(bb).TotalmRNATwo)
                    TotmRNATable(bb+(length(APsubset)),aa)=nan;
                else
                TotmRNATable(bb+(length(APsubset)),aa)=APsubset(bb).TotalmRNATwo;
                end
                MeanzTable(bb+(length(APsubset)),aa)=AvgFluo2;
                StrengthTable(bb+(length(APsubset)),aa)=(VarianceFluo2/AvgFluo2);
                end
            end
            
            end
        end
        %TimeTable(TimeTable==0)=nan;
        WholeTimeTable(cc).Embryo(ee).TimeTable=TimeTable;
        WholeTimeTable(cc).Embryo(ee).TotalmRNATable=TotmRNATable;
        WholeTimeTable(cc).Embryo(ee).MeanFluo=MeanzTable;
        WholeTimeTable(cc).Embryo(ee).StrengthTable=StrengthTable;
        WholeTimeTable(cc).Embryo(ee).EmbryoAvg=nanmean(TimeTable);
        AllNucTimeTable=[AllNucTimeTable;TimeTable];
        AllTotmRNA=[AllTotmRNA;TotmRNATable];
        AllNucStrengthTable=[AllNucStrengthTable; StrengthTable];
        ConTimeAvg=[ConTimeAvg;(nanmean(TimeTable))];
        ConFluoMeanz=[ConFluoMeanz; MeanzTable];
    end
    WholeTimeTable(cc).ConstructAvg=nanmean(ConTimeAvg);
    WholeTimeTable(cc).AllAvgFluo=ConFluoMeanz;
    WholeTimeTable(cc).AllNucs=AllNucTimeTable;
    WholeTimeTable(cc).TotalmRNA=AllTotmRNA;
    WholeTimeTable(cc).AllNucsStrength=AllNucStrengthTable;
    WholeTimeTable(cc).AvgAllNucs=nanmean([WholeTimeTable(cc).AllNucs]);
    WholeTimeTable(cc).SDAllNucs=nanstd([WholeTimeTable(cc).AllNucs]);
    %WholeTimeTable(cc).SEAllNucs=(WholeTimeTable(cc).SDAllNucs)/(sqrt(length(WholeTimeTable(cc).AllNucs)));
     for aa=1:length(APbinID)
         WholeTimeTable(cc).SEAllNucs(aa)=(WholeTimeTable(cc).SDAllNucs(aa))/(sqrt(sum(~isnan(WholeTimeTable(cc).AllNucs(:,aa)))));
     end
     WholeTimeTable(cc).Conf95All=(WholeTimeTable(cc).SEAllNucs).*1.95; 
end

% Get rid of places where only have one data point for a whole AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(WholeTimeTable(cc).AllNucs(:,aa))) ==1
            WholeTimeTable(cc).AllNucs([find(~isnan(WholeTimeTable(cc).AllNucs(:,aa)))],aa)=nan;
            WholeTimeTable(cc).AvgAllNucs(aa)=nan;
        end
        NumberNuclei(cc,aa)=sum(~isnan(WholeTimeTable(cc).AllNucs(:,aa)));
    end
end
%% Plotting 
DistalColor=[1 64 172]./255;
Distal32CColor=[118 180 238] ./ 255;
DoubDistColor=[73 184 253] ./ 255;
ProxColor=[238 123 23]./255;
DoubProxColor=[215 183 58] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[52 119 71]./255;
Both32CColor=[120 195 82] ./ 255;
DoubProx32CColor=[200 150 100] ./ 255;

Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DoubDistColor;
Colors(5).Color=DoubProxColor;
Colors(6).Color=BothColor;
Colors(7).Color=DistalColor;
Colors(8).Color=ProxColor;
Colors(9).Color=BothSepColor;
Colors(10).Color=BothColor;
Colors(11).Color=DoubProxColor;
Colors(12).Color=DistalColor;
Colors(13).Color=ProxColor;
Colors(14).Color=BothColor;
Colors(15).Color=DoubProxColor;
Colors(16).Color=DistalColor;
Colors(17).Color=BothColor;

fontsize=18;
fontname='Helvetica';

WhichAP=input('Which AP bin to use?');
EggLength=APbinID .* 100;

FigDirect=[DropboxFolder filesep 'Figures'];
%% CV vs total mRNA production
for cc=1:length(ConstructList)
figure 
[sorted_x, sorted_x_index] = sort(WholeTimeTable(cc).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(cc).AllNucs(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(cc).Color,'LineWidth',1.5)
set(gca,'fontname',fontname,'fontsize',fontsize);
ylabel('temporal CV');
xlabel('total mRNA production');
title(ConstructList{cc});
end


%% CV vs Total mRNA 32C comparisons 
figure
[sorted_x, sorted_x_index] = sort(WholeTimeTable(1).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(1).AllNucs(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(1).Color,'LineWidth',1.5);

[sorted_x, sorted_x_index] = sort(WholeTimeTable(7).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(7).AllNucs(sorted_x_index),70);
plot(sorted_x, smooth_y, 'Color',Colors(7).Color,'LineWidth',1.5,'LineStyle','-.')
set(gca,'FontSize',fontsize,'FontName',fontname);

figure
[sorted_x, sorted_x_index] = sort(WholeTimeTable(6).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(6).AllNucs(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(6).Color,'LineWidth',1.5);

[sorted_x, sorted_x_index] = sort(WholeTimeTable(10).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(10).AllNucs(sorted_x_index),70);
plot(sorted_x, smooth_y, 'Color',Colors(10).Color,'LineWidth',1.5,'LineStyle','-.')

figure
[sorted_x, sorted_x_index] = sort(WholeTimeTable(6).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(6).AllNucs(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(6).Color,'LineWidth',1.5);
BothRTNoise=[sorted_x,smooth_y];

[sorted_x, sorted_x_index] = sort(WholeTimeTable(17).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(17).AllNucs(sorted_x_index),70);
plot(sorted_x, smooth_y, 'Color',Colors(17).Color,'LineWidth',1.5)
set(gca,'FontSize',fontsize,'FontName',fontname);
xlabel('total mRNA produced');
ylabel('temporal noise');
Both17CNoise=[sorted_x,smooth_y];
[H, pval,KStat]=kstest_2s_2d(BothRTNoise,Both17CNoise);
%%
figure
[sorted_x, sorted_x_index] = sort(WholeTimeTable(1).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(1).AllNucs(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(1).Color,'LineWidth',1.5);
DistNoise=[sorted_x,smooth_y];

[sorted_x, sorted_x_index] = sort(WholeTimeTable(6).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(6).AllNucs(sorted_x_index),70);
plot(sorted_x, smooth_y, 'Color',Colors(6).Color,'LineWidth',1.5)
set(gca,'FontSize',fontsize,'FontName',fontname);
xlabel('total mRNA produced');
ylabel('temporal noise');
BothNoise=[sorted_x,smooth_y];
[H, pval,KStat]=kstest_2s_2d(BothNoise,DistNoise);

figure
[sorted_x, sorted_x_index] = sort(WholeTimeTable(4).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(4).AllNucs(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(4).Color,'LineWidth',1.5);
DoubDistNoise=[sorted_x,smooth_y];

[sorted_x, sorted_x_index] = sort(WholeTimeTable(6).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(6).AllNucs(sorted_x_index),70);
plot(sorted_x, smooth_y, 'Color',Colors(6).Color,'LineWidth',1.5)
set(gca,'FontSize',fontsize,'FontName',fontname);
BothNoise=[sorted_x,smooth_y];
[H, pval,KStat]=kstest_2s_2d(BothNoise,DoubDistNoise);
%% CV vs Total mRNA Homo vs hemizygotes
figure
[sorted_x, sorted_x_index] = sort(WholeTimeTable(1).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(1).AllNucs(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(1).Color,'LineWidth',1.5);

[sorted_x, sorted_x_index] = sort(WholeTimeTable(12).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(12).AllNucs(sorted_x_index),70);
plot(sorted_x, smooth_y, 'Color',Colors(12).Color,'LineWidth',1.5,'LineStyle',':')
set(gca,'FontSize',fontsize,'FontName',fontname);
print('-painters',[FigDirect filesep 'HemiDistTotmRNAvTempCV'],'-dsvg');

figure
[sorted_x, sorted_x_index] = sort(WholeTimeTable(2).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(2).AllNucs(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(2).Color,'LineWidth',1.5);

[sorted_x, sorted_x_index] = sort(WholeTimeTable(13).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(13).AllNucs(sorted_x_index),70);
plot(sorted_x, smooth_y, 'Color',Colors(13).Color,'LineWidth',1.5,'LineStyle',':')
set(gca,'FontSize',fontsize,'FontName',fontname);
print('-painters',[FigDirect filesep 'HemiProxTotmRNAvTempCV'],'-dsvg');

figure
[sorted_x, sorted_x_index] = sort(WholeTimeTable(6).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(6).AllNucs(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(6).Color,'LineWidth',1.5);

[sorted_x, sorted_x_index] = sort(WholeTimeTable(14).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(14).AllNucs(sorted_x_index),70);
plot(sorted_x, smooth_y, 'Color',Colors(14).Color,'LineWidth',1.5,'LineStyle',':')
set(gca,'FontSize',fontsize,'FontName',fontname);
print('-painters',[FigDirect filesep 'HemiBothTotmRNAvTempCV'],'-dsvg');

figure
[sorted_x, sorted_x_index] = sort(WholeTimeTable(5).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(5).AllNucs(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(5).Color,'LineWidth',1.5);

[sorted_x, sorted_x_index] = sort(WholeTimeTable(15).TotalmRNA(:));
smooth_y = smooth(WholeTimeTable(15).AllNucs(sorted_x_index),70);
plot(sorted_x, smooth_y, 'Color',Colors(15).Color,'LineWidth',1.5,'LineStyle',':')
set(gca,'FontSize',fontsize,'FontName',fontname);
print('-painters',[FigDirect filesep 'Hemi2xProxTotmRNAvTempCV'],'-dsvg');
%%
%Look at noise vs mean fluorescence 
CheckNoisevFluo=input('plot noise vs mean?','s');
if CheckNoisevFluo=='y';
    for cc=1:length(ConstructList)
        for ee=1:length(WholeTimeTable(cc).Embryo)
            for aa=1:length(APbinID)
            if sum(~isnan(WholeTimeTable(cc).Embryo(ee).MeanFluo(:,aa)))==1
                WholeTimeTable(cc).Embryo(ee).TimeTable(:,aa)=nan;
                WholeTimeTable(cc).Embryo(ee).MeanFluo(:,aa)=nan;
                WholeTimeTable(cc).Embryo(ee).EmbryoAvg(:,aa)=nan;
            end
            end
        end
    end
    for cc=1:length(ConstructList)
        figure
        for ee=1:length(WholeTimeTable(cc).Embryo)
            scatter(WholeTimeTable(cc).Embryo(ee).MeanFluo(:), WholeTimeTable(cc).Embryo(ee).TimeTable(:));
            hold on 
        end
        xlabel('Avg fluorescence');
        title(ConstructList{cc});
        ylabel('Temporal noise');
        ylim([0 10]);
        xlim([0 70000]);
        
        figure
        for ee=1:length(WholeTimeTable(cc).Embryo)
            scatter(WholeTimeTable(cc).Embryo(ee).MeanFluo(:), WholeTimeTable(cc).Embryo(ee).StrengthTable(:));
            hold on
        end
        xlabel('Avg Fluorescence');
        ylabel('Noise strength');
        title(ConstructList{cc});
    end
    
    % Compare these noise plots across constructs
    figure
    for ee=1:length(WholeTimeTable(1).Embryo);
        scatter(WholeTimeTable(1).Embryo(ee).MeanFluo(:), WholeTimeTable(1).Embryo(ee).TimeTable(:),[], Colors(1).Color);
        hold on 
    end
    for ee=1:length(WholeTimeTable(2).Embryo);
        scatter(WholeTimeTable(2).Embryo(ee).MeanFluo(:), WholeTimeTable(2).Embryo(ee).TimeTable(:),[], Colors(2).Color);
    end 
    xlabel('Avg fluorescence');
    ylabel('Temporal noise');
    
end

LSLines=input('Do least squares line?', 's');
if LSLines=='y'
    for cc=1:length(ConstructList)
        EmbMeanProdVal=[];
        EmbCVVal=[];
        EmbNoiseStrengthVal=[];
        for ee=1:length(WholeTimeTable(cc).Embryo)
            EmbMeanProdVal=[EmbMeanProdVal; WholeTimeTable(cc).Embryo(ee).MeanFluo];
            EmbCVVal=[EmbCVVal; WholeTimeTable(cc).Embryo(ee).TimeTable];
            EmbNoiseStrengthVal=[EmbNoiseStrengthVal; WholeTimeTable(cc).Embryo(ee).StrengthTable];
        end
        figure
        scatter(EmbMeanProdVal(:), EmbCVVal(:));
        lsline;
        xlabel('Mean fluorescence');
        ylabel('Temporal noise');
        title(ConstructList{cc});
        figure 
        scatter(EmbMeanProdVal(:), EmbNoiseStrengthVal(:))
        lsline;
        xlabel('Mean fluorescence');
        ylabel('Noise strength');
        title(ConstructList{cc});
    end
end

%%
%Plot all 
figure 
for cc=1:length(ConstructList)
    plot(EggLength,WholeTimeTable(cc).AvgAllNucs,'Color',Colors(cc).Color,'LineWidth',1.5);
    hold on
end
title('Relative noise across time');
xlabel('% egg length');
ylabel('coefficient of variation');
legend('Distal','Proximal','Both Sep','2x Distal', '2x Proximal', 'Both');
%%
%Singles vs Both sep
figure 
errorbar(EggLength, WholeTimeTable(1).AvgAllNucs, WholeTimeTable(1).Conf95All, 'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
errorbar(EggLength, WholeTimeTable(2).AvgAllNucs, WholeTimeTable(2).Conf95All, 'Color',Colors(2).Color,'LineWidth',1.5);
errorbar(EggLength, WholeTimeTable(3).AvgAllNucs, WholeTimeTable(3).Conf95All, 'Color',Colors(3).Color,'LineWidth',1.5);
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylabel('coefficient of variation');

figure 
errorbar(EggLength, WholeTimeTable(1).AvgAllNucs, WholeTimeTable(1).Conf95All, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(2).AvgAllNucs, WholeTimeTable(2).Conf95All, 'Color',Colors(2).Color,'LineWidth',2.5);
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
%legend('Distal','Proximal');
print('-painters',[FigDirect filesep 'ProxDistTempCV'],'-dsvg');

figure 
errorbar(EggLength, WholeTimeTable(1).AvgAllNucs, WholeTimeTable(1).Conf95All, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(2).AvgAllNucs, WholeTimeTable(2).Conf95All, 'Color',Colors(2).Color,'LineWidth',2.5);
errorbar(EggLength, WholeTimeTable(6).AvgAllNucs, WholeTimeTable(6).Conf95All, 'Color',Colors(6).Color,'LineWidth',2.5);
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
%legend('Distal','Proximal', 'Both');
print('-painters',[FigDirect filesep 'SinglesBothTempCV'],'-dsvg');
%% try median values
figure 
errorbar(EggLength, nanmedian(WholeTimeTable(1).AllNucs), WholeTimeTable(1).Conf95All, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, nanmedian(WholeTimeTable(2).AllNucs), WholeTimeTable(2).Conf95All, 'Color',Colors(2).Color,'LineWidth',2.5);
errorbar(EggLength, nanmedian(WholeTimeTable(6).AllNucs), WholeTimeTable(6).Conf95All, 'Color',Colors(6).Color,'LineWidth',2.5);
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
%legend('Distal','Proximal', 'Both');

%% Hemizygotes
figure 
errorbar(EggLength, WholeTimeTable(1).AvgAllNucs, WholeTimeTable(1).Conf95All, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(12).AvgAllNucs, WholeTimeTable(12).Conf95All, 'Color',Colors(12).Color,'LineWidth',2.5,'LineStyle',':');
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
%legend('Distal','Proximal', 'Both');
print('-painters',[FigDirect filesep 'HemiDistTempCV'],'-dsvg');

figure 
errorbar(EggLength, WholeTimeTable(2).AvgAllNucs, WholeTimeTable(2).Conf95All, 'Color',Colors(2).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(13).AvgAllNucs, WholeTimeTable(13).Conf95All, 'Color',Colors(13).Color,'LineWidth',2.5,'LineStyle',':');
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
%legend('Distal','Proximal', 'Both');
print('-painters',[FigDirect filesep 'HemiProxTempCV'],'-dsvg');

figure 
errorbar(EggLength, WholeTimeTable(6).AvgAllNucs, WholeTimeTable(6).Conf95All, 'Color',Colors(6).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(14).AvgAllNucs, WholeTimeTable(14).Conf95All, 'Color',Colors(14).Color,'LineWidth',2.5,'LineStyle',':');
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
%legend('Distal','Proximal', 'Both');
print('-painters',[FigDirect filesep 'HemiBothTempCV'],'-dsvg');

figure 
errorbar(EggLength, WholeTimeTable(5).AvgAllNucs, WholeTimeTable(5).Conf95All, 'Color',Colors(5).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(15).AvgAllNucs, WholeTimeTable(15).Conf95All, 'Color',Colors(15).Color,'LineWidth',2.5,'LineStyle',':');
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
%legend('Distal','Proximal', 'Both');
print('-painters',[FigDirect filesep 'Hemi2xProxTempCV'],'-dsvg');
%%
%Duplicates vs SE
figure 
errorbar(EggLength, WholeTimeTable(1).AvgAllNucs, WholeTimeTable(1).Conf95All, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, WholeTimeTable(2).AvgAllNucs, WholeTimeTable(2).Conf95All, 'Color',Colors(2).Color,'LineWidth',2.5); 
errorbar(EggLength, WholeTimeTable(5).AvgAllNucs, WholeTimeTable(5).Conf95All, 'Color',Colors(5).Color,'LineWidth',2.5);
errorbar(EggLength, WholeTimeTable(6).AvgAllNucs, WholeTimeTable(6).Conf95All, 'Color',Colors(6).Color,'LineWidth',2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
print('-painters',[FigDirect filesep 'AllTempCV'],'-dsvg');


figure 
errorbar(EggLength, WholeTimeTable(4).AvgAllNucs, WholeTimeTable(4).Conf95All, 'Color',Colors(4).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(6).AvgAllNucs, WholeTimeTable(6).Conf95All, 'Color',Colors(6).Color,'LineWidth',2.5);
%title('Relative noise across time');
xlabel('% egg length');
xlim([0 100]);
set(gca,'fontsize', fontsize, 'fontname', fontname);
ylabel('coefficient of variation');
%legend('2x Distal', 'Both');

figure 
errorbar(EggLength, WholeTimeTable(4).AvgAllNucs, WholeTimeTable(4).Conf95All, 'Color',Colors(4).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, WholeTimeTable(5).AvgAllNucs, WholeTimeTable(5).Conf95All, 'Color',Colors(5).Color,'LineWidth',2.5);
errorbar(EggLength, WholeTimeTable(6).AvgAllNucs, WholeTimeTable(6).Conf95All, 'Color',Colors(6).Color,'LineWidth',2.5);
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
%legend('2x Proximal', 'Both');
print('-painters',[FigDirect filesep '2EnhancervBothTempCV'],'-dsvg');

%Zoom in to center 20%
figure 
errorbar(EggLength, WholeTimeTable(4).AvgAllNucs, WholeTimeTable(4).Conf95All, 'Color',Colors(4).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, WholeTimeTable(5).AvgAllNucs, WholeTimeTable(5).Conf95All, 'Color',Colors(5).Color,'LineWidth',2.5);
errorbar(EggLength, WholeTimeTable(6).AvgAllNucs, WholeTimeTable(6).Conf95All, 'Color',Colors(6).Color,'LineWidth',2.5);
%title('Relative noise across time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([40 60]);
ylim([0 5]);
ylabel('coefficient of variation');
%legend('2x Proximal', 'Both');

%%
%Singles vs doubles
figure 
errorbar(EggLength, WholeTimeTable(1).AvgAllNucs, WholeTimeTable(1).SEAllNucs, 'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
errorbar(EggLength, WholeTimeTable(4).AvgAllNucs, WholeTimeTable(4).SEAllNucs, 'Color',Colors(4).Color,'LineWidth',1.5);
%title('Relative noise across time');
xlabel('% egg length');
xlim([0 100]);
set(gca, 'fontsize',fontsize, 'fontname', fontname);
ylabel('coefficient of variation');
legend('Distal', '2x Distal');

figure 
errorbar(EggLength, WholeTimeTable(2).AvgAllNucs, WholeTimeTable(2).SEAllNucs, 'Color',Colors(2).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(5).AvgAllNucs, WholeTimeTable(5).SEAllNucs, 'Color',Colors(5).Color,'LineWidth',2.5);
%title('relative noise across time');
xlabel('% egg length');
xlim([0 100]);
ylabel('coefficient of variation');
%legend('Proximal', '2x Proximal');
set(gca, 'fontname', fontname, 'fontsize', fontsize);

figure 
errorbar(EggLength, WholeTimeTable(3).AvgAllNucs, WholeTimeTable(3).SEAllNucs, 'Color',Colors(3).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(6).AvgAllNucs, WholeTimeTable(6).SEAllNucs, 'Color',Colors(6).Color,'LineWidth',2.5);
%title('relative noise across time');
xlabel('% egg length');
xlim([0 100]);
ylabel('coefficient of variation');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
legend('Both Sep', 'Both');
%%
%Temperature alterations
figure 
errorbar(EggLength, WholeTimeTable(1).AvgAllNucs, WholeTimeTable(1).SEAllNucs, 'Color',Colors(1).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(7).AvgAllNucs, WholeTimeTable(7).SEAllNucs, 'Color',Colors(7).Color,'LineWidth',2.5,'LineStyle','-.');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
%legend('Distal', 'Distal 32C');
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
print('-painters',[FigDirect filesep 'DistTCompTempCV'],'-dsvg');

figure 
errorbar(EggLength, WholeTimeTable(2).AvgAllNucs, WholeTimeTable(2).SEAllNucs, 'Color',Colors(2).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(8).AvgAllNucs, WholeTimeTable(8).SEAllNucs, 'Color',Colors(8).Color,'LineWidth',2.5,'LineStyle','-.');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
print('-painters',[FigDirect filesep 'ProxTCompTempCV'],'-dsvg');

figure 
errorbar(EggLength, WholeTimeTable(3).AvgAllNucs, WholeTimeTable(3).SEAllNucs, 'Color',Colors(3).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(9).AvgAllNucs, WholeTimeTable(9).SEAllNucs, 'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-.');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
%legend('Both Sep', 'Both Sep 32C');
xlabel('% egg length');
xlim([0 100]);
ylabel('coefficient of variation');

figure 
errorbar(EggLength, WholeTimeTable(6).AvgAllNucs, WholeTimeTable(6).SEAllNucs, 'Color',Colors(6).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(10).AvgAllNucs, WholeTimeTable(10).SEAllNucs, 'Color',Colors(10).Color,'LineWidth',2.5,'LineStyle','-.');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
%legend('Both','Both 32C');
print('-painters',[FigDirect filesep 'BothTCompTempCV'],'-dsvg');


figure 
errorbar(EggLength, WholeTimeTable(5).AvgAllNucs, WholeTimeTable(5).SEAllNucs, 'Color',Colors(5).Color,'LineWidth',2.5);
hold on 
errorbar(EggLength, WholeTimeTable(11).AvgAllNucs, WholeTimeTable(11).SEAllNucs, 'Color',Colors(11).Color,'LineWidth',2.5,'LineStyle','-.');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
xlim([0 100]);
ylim([0 9]);
ylabel('coefficient of variation');
%legend('2x Proximal','2x Proximal 32C');
print('-painters',[FigDirect filesep '2ProxTCompTempCV'],'-dsvg');

%%
%Do single ANOVA of constructs at specified position
for cc=1:length(ConstructList)
    MaxAPtoUse=length(WholeTimeTable(cc).AllNucs);
end
MaxAPtoUse=max(MaxAPtoUse);
APbox=[];
for cc=1:length(ConstructList)
    for bb=1:length(WholeTimeTable(cc).AllNucs)
        APbox(bb,cc)=WholeTimeTable(cc).AllNucs(bb,WhichAP);
    end
    if length(WholeTimeTable(cc).AllNucs) < MaxAPtoUse
        APbox([end:MaxAPtoUse],cc)=nan;
    end
end
[p,tbl,stats]=anova1(APbox);
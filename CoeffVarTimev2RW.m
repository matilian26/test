%% calculate coefficient of variance for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info   for right now just going to use nc14
ncUse=input('Want to only use nc14?','s');

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
    for ee=1:NEmbryos
        
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        Filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            Filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end
        load(Filename);
        %if ncUse=='y'
%         nc_number=[CompiledParticles.nc];
%         CompiledParticles_14=CompiledParticles(nc_number==14);
        APstuff=[BurstProperties.APBin];
        TimeTable=zeros(1000,41);
        for aa=1:length(APbinID)
            APsubset=[];
            VarPart=[];
            SumVarPart=[];
            APsubset=BurstProperties(APstuff==APbinID(aa));
            if ~isempty(APsubset)
               
            for bb=1:length(APsubset)
                TimeLength=[];
                TimeLength=length(APsubset(bb).FluoFrames);
                AvgFluo=nanmean([APsubset(bb).Fluo]);
                for tt=1:TimeLength
                    VarPart(tt)=((APsubset(bb).Fluo(tt) - AvgFluo)^2);
                end
            SumVarPart=sum(VarPart);
            TimeTable(bb,aa)=(sqrt((1/TimeLength)*(SumVarPart))/AvgFluo);
            %FluoValues(cc).Embryo(ee).Nucleus(bb).FluoValue=FluoVals;
            end
            
            end
        end
        TimeTable(TimeTable==0)=nan;
        WholeTimeTable(cc).Embryo(ee).TimeTable=TimeTable;
        WholeTimeTable(cc).Embryo(ee).EmbryoAvg=nanmean(TimeTable);
        AllNucTimeTable=[AllNucTimeTable;TimeTable];
        ConTimeAvg=[ConTimeAvg;(nanmean(TimeTable))];
    end
    WholeTimeTable(cc).ConstructAvg=nanmean(ConTimeAvg);
    WholeTimeTable(cc).AllNucs=AllNucTimeTable;
    WholeTimeTable(cc).AvgAllNucs=nanmean([WholeTimeTable(cc).AllNucs]);
    WholeTimeTable(cc).SDAllNucs=nanstd([WholeTimeTable(cc).AllNucs]);
    WholeTimeTable(cc).SEAllNucs=(WholeTimeTable(cc).SDAllNucs)/(sqrt(length(WholeTimeTable(cc).AllNucs)));
end
        
%% Plotting 
DistalColor=[8 180 238] ./ 255;
DoubDistColor=[1 17 181] ./ 255;
ProxColor=[251 230 60] ./ 255;
DoubProxColor=[251 190 80] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothColor=[12 195 82] ./ 255;

Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DoubDistColor;
Colors(5).Color=DoubProxColor;
Colors(6).Color=BothColor;

fontsize=18;
fontname='Helvetica';

WhichAP=input('Which AP bin to use?');
EggLength=APbinID .* 100;

for cc=1:length(ConstructList)
    for ee=1:length(WholeTimeTable(cc).Embryo)
        APTimeBox(ee,cc)=nanmean(WholeTimeTable(cc).Embryo(ee).TimeTable(:,WhichAP));
    end
end
APTimeBox(APTimeBox==0)=nan;
figure
for cc=1:length(ConstructList)
    for ee=1:length(WholeTimeTable(cc).Embryo)
        AvgTimeatAP=[];
        AvgTimeatAP=nanmean(WholeTimeTable(cc).Embryo(ee).TimeTable(:,WhichAP));
        plot(cc,AvgTimeatAP,'o','Color',Colors(cc).Color);
        hold on 
    end
end
boxplot(APTimeBox,'Colors','k')
xticks([1:6]);
xticklabels({'Distal', 'Proximal', 'Both Sep', '2x Distal', '2x Proximal', 'Both'});
xlabel('Construct');
ylabel('Mean coefficient of variation');
title(['Noise across time',' ', num2str(EggLength(WhichAP)), '% egg length']);

%plot mean time variance of each construct vs AP position 
figure 
for cc=1:length(ConstructList)
    plot(EggLength, WholeTimeTable(cc).AvgAllNucs,'Color',Colors(cc).Color,'LineWidth',1.5);
    hold on 
end
legend('Distal','Proximal','Both Sep', '2x Distal', '2x Proximal', 'Both')
xlabel('% egg length')
xlim([0 100]);
ylabel('Coefficient of variation');
title('Avg noise across time');

%Singles vs SE
figure
errorbar(EggLength,WholeTimeTable(1).AvgAllNucs,WholeTimeTable(1).SEAllNucs,'Color',Colors(1).Color,'LineWidth',1.5);
hold on
errorbar(EggLength,WholeTimeTable(2).AvgAllNucs,WholeTimeTable(2).SEAllNucs,'Color',Colors(2).Color,'LineWidth',1.5);
errorbar(EggLength,WholeTimeTable(6).AvgAllNucs,WholeTimeTable(6).SEAllNucs,'Color',Colors(6).Color,'LineWidth',1.5);
%title('Relative expression noise in time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Coefficient of variation');
xlim([0 100]);
%legend('Distal', 'Proximal', 'Both');

% plot doubles compared to SE
figure
errorbar(EggLength,WholeTimeTable(4).AvgAllNucs,WholeTimeTable(4).SEAllNucs,'Color',Colors(4).Color,'LineWidth',1.5);
hold on
errorbar(EggLength,WholeTimeTable(5).AvgAllNucs,WholeTimeTable(5).SEAllNucs,'Color',Colors(5).Color,'LineWidth',1.5);
errorbar(EggLength,WholeTimeTable(6).AvgAllNucs,WholeTimeTable(6).SEAllNucs,'Color',Colors(6).Color,'LineWidth',1.5);
%title('Relative expression noise in time');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Coefficient of variation');
xlim([0 100]);
%legend('2x Dist', '2x Prox', 'Both');

figure 
errorbar(EggLength, WholeTimeTable(4).AvgAllNucs,WholeTimeTable(4).SEAllNucs,'Color',Colors(4).Color,'LineWidth',1.5);
hold on 
errorbar(EggLength, WholeTimeTable(6).AvgAllNucs, WholeTimeTable(6).SEAllNucs,'Color',Colors(6).Color,'LineWidth',1.5)
legend('2x Distal', 'Both');
xlabel('% egg length')
ylabel('Coefficient of variation');
title('Relative noise across time');
xlim([0 100]);


figure 
errorbar(EggLength, WholeTimeTable(5).AvgAllNucs,WholeTimeTable(5).SEAllNucs,'Color',Colors(5).Color,'LineWidth',1.5);
hold on 
errorbar(EggLength, WholeTimeTable(6).AvgAllNucs,WholeTimeTable(6).SEAllNucs,'Color',Colors(6).Color,'LineWidth',1.5)
%legend('2x Proximal', 'Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length')
ylabel('Coefficient of variation');
title('Relative noise across time');
xlim([0 100]);

%Doubles vs singles 
figure 
errorbar(EggLength, WholeTimeTable(1).AvgAllNucs,WholeTimeTable(1).SEAllNucs,'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
errorbar(EggLength, WholeTimeTable(4).AvgAllNucs,WholeTimeTable(4).SEAllNucs,'Color',Colors(4).Color,'LineWidth',1.5)
legend('Distal', '2x Distal');
xlabel('% egg length')
ylabel('Coefficient of variation');
title('Relative noise across time');
xlim([0 100]);

figure 
errorbar(EggLength, WholeTimeTable(2).AvgAllNucs,WholeTimeTable(2).SEAllNucs,'Color',Colors(2).Color,'LineWidth',1.5);
hold on 
plot(EggLength, WholeTimeTable(5).AvgAllNucs,WholeTimeTable(5).SEAllNucs,'Color',Colors(5).Color,'LineWidth',1.5)
legend('Proximal', '2x Proximal');
xlabel('% egg length')
ylabel('Coefficient of variation');
title('Relative noise across time');
xlim([0 100]);

% ANOVA comparing constructs at specified AP position 
for cc=1:length(ConstructList)
    ConAPTimes=[];
    for ee=1:length(WholeTimeTable(cc).Embryo);
        ConAPTimes=[ConAPTimes;WholeTimeTable(cc).Embryo(ee).TimeTable(:,WhichAP)];
    end
    LengthConTimes(cc)=size(ConAPTimes,1);
end
MaxTimeLength=max(LengthConTimes);
AllConsAPTimes=[];
for cc=1:length(ConstructList)
    ConAPTimes=[];
    for ee=1:length(WholeTimeTable(cc).Embryo)
        ConAPTimes=[ConAPTimes;WholeTimeTable(cc).Embryo(ee).TimeTable(:,WhichAP)];
    end
    if size(ConAPTimes,1) < MaxTimeLength
        ConAPTimes([end:MaxTimeLength],:)=nan;
    end
    AllConsAPTimes=[AllConsAPTimes, ConAPTimes];
end
[p,tbl,stats]=anova1(AllConsAPTimes);

%Saving temporal noise 
save([DropboxFolder filesep 'Constructs' filesep 'AllTemporalNoise'],'WholeTimeTable');

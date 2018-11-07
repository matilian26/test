%% Calculate construct "decay" and "onset" rates 

ConstructList= {'KrDist','KrProx','KrBothSep','KrDistEmpty','KrProxEmpty','KrDistDuplicN','KrProxDuplic','Kr2xProxEmpty','KrBoth','KrBothEmpty'} %{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
%Count for each construct
AvgDecay=[];
AvgOnset=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgDecayCon=[];
    AvgOnsetCon=[];
    firsttime=1;
    ConDecayAllAP=[];
    ConOnsetAllAP=[];
    ConDecaySE=[];
    ConOnsetSE=[];
    ConDecaySD=[];
    ConOnsetSD=[];
    ConOnsetSE=[];
    EmbsArray=[];
    %ConDurSD=[];
     DecayAllAP=[];
     OnsetAllAP=[];
     mRNAProdErrorAllAP=[];

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
       
        for aa=1:length(APbinID)
            DecayAP=[];
            OnsetAP=[];
            ProdErrorAP=[];
            DecayValAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(DecayValAP)
                DecayAP=[DecayAP; nan];
                OnsetAP=[OnsetAP;nan];
            else
            for bb=1:length(DecayValAP)
                DecayAP=[DecayAP;[BurstProperties(DecayValAP(bb)).Decay]'];  %put all mRNA outputs at a given AP value in a column going down
                %ProdErrorAP=[ProdErrorAP;BurstProperties(DecayValAP(bb)).TotalmRNAError];           
                OnsetAP=[OnsetAP;[BurstProperties(DecayValAP(bb)).Onset]'];
            end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(DecayAP)
                DecayAllAP(bb,aa,ee)=DecayAP(bb);
                 OnsetAllAP(bb,aa,ee)=OnsetAP(bb);
            end
            
%             for bb=1:length(ProdErrorAP)
%                 mRNAProdErrorAllAP(bb,aa,ee)=ProdErrorAP(bb);
%             end
            DecayAllAP(DecayAllAP==0)=nan;
            OnsetAllAP(OnsetAllAP==0)=nan;
            %mRNAProdErrorAllAP(mRNAProdErrorAllAP==0)=nan;
            
            DecaySD(ee,aa,cc)=nanstd(DecayAllAP(:,aa,ee));
            OnsetSD(ee,aa,cc)=nanstd(OnsetAllAP(:,aa,ee));
        DecaySE(ee,aa,cc)=DecaySD(ee,aa,cc)/sqrt(length(DecayValAP));   
        OnsetSE(ee,aa,cc)=OnsetSD(ee,aa,cc)/sqrt(length(DecayValAP));
            clear DecayValAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        
         EmbryoAvgDecay(ee).MeanProd=nanmean(DecayAllAP(:,:,ee));
         %EmbryoAvgDecay(ee).ProdSE=nanmean(mRNAProdErrorAllAP(:,:,ee));  %no idea really if this is right
         AvgDecayAllAP(cc).Embryos(ee).MeanDecay=nanmean(DecayAllAP(:,:,ee));
        AvgDecayAllAP(cc).Embryos(ee).SE=DecaySE(ee,:,cc);
         AvgDecayAllAP(cc).Embryos(ee).SD=DecaySD(ee,:,cc);
         AvgDecayAllAP(cc).Embryos(ee).MeanOnset=nanmean(OnsetAllAP(:,:,ee));
         AvgDecayAllAP(cc).Embryos(ee).SEOnset=OnsetSE(ee,:,cc);
         AvgDecayAllAP(cc).Embryos(ee).SDOnset=OnsetSD(ee,:,cc);
         
    end
    for bb=1:size(DecayAllAP,3)
            ConDecayAllAP=[ConDecayAllAP; DecayAllAP(:,:,bb)];
            ConOnsetAllAP=[ConOnsetAllAP; OnsetAllAP(:,:,bb)];
        end
        for bb=1:size(DecaySD,3)
            ConDecaySD=[ConDecaySD; DecaySD(:,:,bb)];
            ConOnsetSD=[ConOnsetSD; OnsetSD(:,:,bb)];
        end
        for bb=1:size(DecaySE,3)
            ConDecaySE=[ConDecaySE;DecaySE(:,:,bb)];
            ConOnsetSE=[ConOnsetSE; OnsetSE(:,:,bb)];
        end
%         EmbsArray=[];
%         for ee=1:NEmbryos
%             EmbsArray=[EmbsArray; [EmbryoAvgProd(ee).MeanProd]];  %3/16 10pm not working, not sure why, trying to put mean at each ap bin for each embryo in one array on top of one another 
%         end
        AvgDecayAllAP(cc).AvgDecay=nanmean(ConDecayAllAP,1);  %Avg mRNA production of all embryos of a construct by AP position
        AvgDecayAllAP(cc).ConSD=nanmean(ConDecaySD,1);
        AvgDecayAllAP(cc).ConSE=nanmean(ConDecaySE,1);
        AvgDecayAllAP(cc).AllDecay=[ConDecayAllAP];
        AvgDecayAllAP(cc).AllSD=nanstd([AvgDecayAllAP(cc).AllDecay]);
        AvgDecayAllAP(cc).AllSE=(([AvgDecayAllAP(cc).AllSD])/sqrt(length(AvgDecayAllAP(cc).AllDecay)));
        AvgDecayAllAP(cc).All95Conf=(AvgDecayAllAP(cc).AllSE).*1.95;
        AvgDecayAllAP(cc).AvgOnset=nanmean(ConOnsetAllAP,1);
        AvgDecayAllAP(cc).AllOnset=[ConOnsetAllAP];
        AvgDecayAllAP(cc).AllOnsetSD=nanstd([AvgDecayAllAP(cc).AllOnset]);
        AvgDecayAllAP(cc).AllOnsetSE=(([AvgDecayAllAP(cc).AllOnsetSD])/sqrt(length(AvgDecayAllAP(cc).AllOnset)));
        AvgDecayAllAP(cc).AllOnset95Conf=(AvgDecayAllAP(cc).AllOnsetSE).*1.95;
        %for ee=1:length(EmbryoAvgProd)
        AvgDecayAllAP(cc).EmbryoMeans=[EmbsArray];
        %end 
        clear DecayAllAP OnsetAllAP;
end

% Get rid of single values for an AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(AvgDecayAllAP(cc).AllDecay(:,aa)))==1
            AvgDecayAllAP(cc).AllSD(aa)=nan;
            AvgDecayAllAP(cc).AvgDecay(aa)=nan;
            AvgDecayAllAP(cc).AllSE(aa)=nan;
            AvgDecayAllAP(cc).All95Conf(aa)=nan;
        end
    end
end

%% Plotting
EggLength=APbinID.*100;

DistalColor=[8 180 238] ./ 255;
DistalEmptyColor=[8 210 238] ./ 255; 
DoubDistColor=[1 17 181] ./ 255;
ProxColor=[251 220 50] ./ 255;
ProxEmptyColor=[251 250 50] ./255;
DoubProxColor=[251 190 100] ./ 255;
DoubProxEmptyColor=[251 220 50] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothColor=[12 195 82] ./ 255;
BothEmptyColor=[12 250 100] ./ 255;

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

fontsize=18;
fontname='Helvetica';

figure
for cc=1:length(ConstructList)
    plot(EggLength, nanmean(abs(AvgDecayAllAP(cc).AllDecay)),'Color',Colors(cc).Color,'LineWidth',1.5)
    hold on
end
xlabel('% Egg length');
ylabel('Decay time (min)');

figure
for cc=1:length(ConstructList)
    plot(EggLength, nanmean(abs(AvgDecayAllAP(cc).AllDecay)),'Color',Colors(cc).Color,'LineWidth',1.5)
    hold on
end
xlabel('% Egg length');
ylabel('Onset time (min)');


% Singles
figure
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(1).AllDecay)),AvgDecayAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(2).AllDecay)),AvgDecayAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
xlabel('% egg length');
xlim([0 100]);
ylabel('Fluorescence decay time (min)');
set(gca, 'FontSize', fontsize, 'FontName', fontname);


figure
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(1).AllOnset)),AvgDecayAllAP(1).AllOnset95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(2).AllOnset)),AvgDecayAllAP(2).AllOnset95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
xlabel('% egg length');
xlim([0 100]);
ylabel('Fluorescence onset time (min)');
set(gca, 'FontSize', fontsize, 'FontName', fontname);

figure
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(1).AllDecay)),AvgDecayAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(2).AllDecay)),AvgDecayAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(3).AllDecay)),AvgDecayAllAP(3).All95Conf,'Color',Colors(3).Color,'LineWidth',2.5);
xlabel('% egg length');
xlim([0 100]);
ylabel('Fluorescence decay time (min)');
set(gca, 'FontSize', fontsize, 'FontName', fontname);

figure
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(1).AllDecay)),AvgDecayAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(2).AllDecay)),AvgDecayAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(9).AllDecay)),AvgDecayAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
xlabel('% egg length');
xlim([0 100]);
ylabel('Fluorescence decay time (min)');
set(gca, 'FontSize', fontsize, 'FontName', fontname);

figure
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(1).AllDecay)),AvgDecayAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(4).AllDecay)),AvgDecayAllAP(4).All95Conf,'Color',Colors(4).Color,'LineWidth',2.5);
xlabel('% egg length');
xlim([0 100]);
legend('Distal','1x Distal');
ylabel('Fluorescence decay time (min)');
set(gca, 'FontSize', fontsize, 'FontName', fontname);

figure
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(2).AllDecay)),AvgDecayAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(5).AllDecay)),AvgDecayAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',2.5);
xlabel('% egg length');
xlim([0 100]);
label('Proximal','1x Proximal');
ylabel('Fluorescence decay time (min)');
set(gca, 'FontSize', fontsize, 'FontName', fontname);

% Doubles
figure
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(6).AllDecay)),AvgDecayAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(7).AllDecay)),AvgDecayAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
errorbar(EggLength, nanmean(abs(AvgDecayAllAP(9).AllDecay)),AvgDecayAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
xlabel('% egg length');
xlim([0 100]);
ylabel('Fluorescence decay time (min)');
set(gca, 'FontSize', fontsize, 'FontName', fontname);

figure
errorbar(EggLength, (1./nanmean(abs(AvgDecayAllAP(7).AllDecay))),AvgDecayAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
hold on
errorbar(EggLength, (1./nanmean(abs(AvgDecayAllAP(8).AllDecay))),AvgDecayAllAP(8).All95Conf,'Color',Colors(8).Color,'LineWidth',2.5);
xlabel('% egg length');
xlim([0 100]);
legend('2x Proximal', 'Single 2x Proximal');
ylabel('Fluorescence decay time (min)');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
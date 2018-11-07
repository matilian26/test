% % %% calculate average mRNA produced per nucleus for each construct 
% %load constructs
% ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
%     %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
%     
% [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
%  Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
%  ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
% 
% %Ask if only want nc14 info
% ncUse=input('Want to only use nc14?','s');
% 
% %First make sure have done SmoothingStuffnc14RW for all embryos
% for cc=1:length(ConstructList)
%     Data= LoadMS2SetsCS(ConstructList{cc});
%     Datalength(cc)=length(Data);
%     NEmbryos = length(Data);
%     for ee=1:NEmbryos
%         PrefixName=Data(ee).Prefix;
%         
%         firstfilename=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
%         load(firstfilename);
%         if ncUse=='y'
%             run('SmoothingStuffnc14RW.m')
%             save([DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14'],'BurstProperties')
%         %else
%             %run SmoothingStuffRW.m
%         %end
%     end
%     %clear BurstProperties
%     end
% end


%%
ConstructList= {'KrDist','KrProx','KrBothSep','KrDistEmpty','KrProxEmpty','KrDistDuplicN','KrProxDuplic','Kr2xProxEmpty','KrBoth','KrBothEmpty','KrDist32C','KrProx32C','KrBothSep32C','KrBoth32C','Kr2xProx32C','HbEmpty','KrDistDuplicN','KrDist17C','Kr2xDist32C','KrBoth17C'} %{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
SlopeUse=input('Want to use Slope calculations?','s');
%Count for each construct
AvgmRNAProd=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgmRNAProdCon=[];
    firsttime=1;
    Timez=[];
    ConmRNAProdAllAP=[];
    ConmRNAProdSE=[];
    ConProdSD=[];
    EmbsArray=[];
    %ConDurSD=[];
     mRNAProdAllAP=[];
     mRNAProdErrorAllAP=[];

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
        CompPars=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
        load(CompPars);
        
        NumberBursts(ee,cc)=length([BurstProperties.Duration]);
        %seperate out by AP bin
       
        for aa=1:length(APbinID)
            ProdAP=[];
            ProdErrorAP=[];
            ProductionAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(ProductionAP)
                ProdAP=[ProdAP; nan];
            else
            for bb=1:length(ProductionAP)
                ProdAP=[ProdAP;BurstProperties(ProductionAP(bb)).TotalmRNA];  %put all mRNA outputs at a given AP value in a column going down
                ProdErrorAP=[ProdErrorAP;BurstProperties(ProductionAP(bb)).TotalmRNAError];           
            end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(ProdAP)
                mRNAProdAllAP(bb,aa,ee)=ProdAP(bb);
            end
            for bb=1:length(ProdErrorAP)
                mRNAProdErrorAllAP(bb,aa,ee)=ProdErrorAP(bb);
            end
            mRNAProdAllAP(mRNAProdAllAP==0)=nan;
            mRNAProdErrorAllAP(mRNAProdErrorAllAP==0)=nan;
            
            ProductionSD(ee,aa,cc)=nanstd(mRNAProdAllAP(:,aa,ee));
        %ProductionSE(ee,aa,cc)=ProductionSD(ee,aa,cc)/sqrt(length(ProductionAP)); 
        ProductionSE(ee,aa,cc)=ProductionSD(ee,aa,cc)/sqrt(sum(~isnan(ProdAP)));
            clear ProductionAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        AvgProdAllAP(cc).nc14Time(ee)=ElapsedTime(end)-ElapsedTime(nc14);
        
         EmbryoAvgProd(ee).MeanProd=nanmean(mRNAProdAllAP(:,:,ee));
         EmbryoAvgProd(ee).ProdSE=nanmean(mRNAProdErrorAllAP(:,:,ee));  %no idea really if this is right
         AvgProdAllAP(cc).EmbryosProd(ee).MeanProd=nanmean(mRNAProdAllAP(:,:,ee));
        AvgProdAllAP(cc).EmbryosProd(ee).SE=ProductionSE(ee,:,cc);
         AvgProdAllAP(cc).EmbryosProd(ee).SD=ProductionSD(ee,:,cc);
         Timez=[Timez,(ElapsedTime(end)-ElapsedTime(nc14))];
    end
    for bb=1:size(mRNAProdAllAP,3)
            ConmRNAProdAllAP=[ConmRNAProdAllAP; mRNAProdAllAP(:,:,bb)];
        end
        for bb=1:size(ProductionSD,3)
            ConProdSD=[ConProdSD; ProductionSD(:,:,bb)];
        end
        for bb=1:size(ProductionSE,3)
            ConmRNAProdSE=[ConmRNAProdSE;ProductionSE(:,:,bb)];
        end
%         EmbsArray=[];
%         for ee=1:NEmbryos
%             EmbsArray=[EmbsArray; [EmbryoAvgProd(ee).MeanProd]];  %3/16 10pm not working, not sure why, trying to put mean at each ap bin for each embryo in one array on top of one another 
%         end
        AvgProdAllAP(cc).nc14Time=mean(Timez);
        AvgProdAllAP(cc).AvgProd=nanmean(ConmRNAProdAllAP,1);  %Avg mRNA production of all embryos of a construct by AP position
        AvgProdAllAP(cc).ConSD=nanmean(ConProdSD,1);
        AvgProdAllAP(cc).ConSE=nanmean(ConmRNAProdSE,1);
        AvgProdAllAP(cc).AllProds=[ConmRNAProdAllAP];
        AvgProdAllAP(cc).AllSD=nanstd([AvgProdAllAP(cc).AllProds]);
        %AvgProdAllAP(cc).AllSE=(([AvgProdAllAP(cc).AllSD])/sqrt(length(AvgProdAllAP(cc).AllProds)));
        for aa=1:length(APbinID)
            AvgProdAllAP(cc).AllSE(aa)=(([AvgProdAllAP(cc).AllSD(aa)])./sqrt(sum(~isnan(AvgProdAllAP(cc).AllProds(:,aa)))));
        end
        AvgProdAllAP(cc).All95Conf=(AvgProdAllAP(cc).AllSE).*1.95;
        %for ee=1:length(EmbryoAvgProd)
        AvgProdAllAP(cc).EmbryoMeans=[EmbsArray];
        %end 
        clear mRNAProdAllAP mRNAProdErrorAllAP;
end

% Get rid of single values for an AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(AvgProdAllAP(cc).AllProds(:,aa)))==1
            AvgProdAllAP(cc).AllSD(aa)=nan;
            AvgProdAllAP(cc).AvgProd(aa)=nan;
            AvgProdAllAP(cc).AllSE(aa)=nan;
            AvgProdAllAP(cc).All95Conf(aa)=nan;
        end
    end
end

%% Coefficent of Variance
% for cc=1:length(ConstructList)
% CoeffVar(cc).Means=nanmean([AvgProdAllAP(cc).AllProds]);
% CoeffVar(cc).SD=nanstd([AvgProdAllAP(cc).AllProds],1);
% CoeffVar(cc).CoeffVar=((CoeffVar(cc).SD)./(CoeffVar(cc).Means));
% 
% CoeffVar(cc).Error=((CoeffVar(cc).SD)/(sqrt(length(AvgProdAllAP(cc).AllProds))));
% end
% EggLength=APbinID .* 100;
% DistalColor=[8 180 238] ./ 255;
% DoubDistColor=[1 17 181] ./ 255;
% ProxColor=[251 230 60] ./ 255;
% DoubProxColor=[251 190 80] ./ 255;
% BothSepColor=[94 250 81] ./ 255;
% BothColor=[12 195 82] ./ 255;
% 
% Colors(1).Color=DistalColor;
% Colors(2).Color=ProxColor;
% Colors(3).Color=BothSepColor;
% Colors(4).Color=DoubDistColor;
% Colors(5).Color=DoubProxColor;
% Colors(6).Color=BothColor;
% 
% figure 
% for cc=1:length(ConstructList)
%     errorbar(EggLength,CoeffVar(cc).CoeffVar,CoeffVar(cc).Error,'Color',Colors(cc).Color,'LineWidth',1.5);
%     hold on 
% end
% figure
% for cc=1:length(ConstructList)
%     plot(EggLength,CoeffVar(cc).CoeffVar,'Color',Colors(cc).Color, 'LineWidth',1.5);
%     hold on 
% end
    
    
    

%% Plotting 
%Plot mean of each embryo for each construct
EggLength=APbinID.*100;
FractInfo=[DropboxFolder filesep 'Constructs' filesep 'FractON.mat'];
% load(FractInfo);
% AvgTimenc14=(377*(mean([AvgProdAllAP(1:14).nc14Time])));
% for cc=1:length(ConstructList)
%     AvgProdAllAP(cc).MS2Conv=(377*(AvgTimenc14*AvgFractONAllAP(cc).AvgFractON));
% end
% MS2Convnc13=(377*15);%(377*15);
% MS2Convnc14=(377*45);
MS2Conversion=(377*((1.301+5.630)/1.5)); %Frnap * Telongation ((length of MS2 +rest of transcript) / elongation rate)
HbMS2Conversion=(377*((1.275+4.021)/1.5));

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
Colors(4).Color=DistalColor;
Colors(5).Color=ProxColor;
Colors(6).Color=DoubDistColor; 
Colors(7).Color=DoubProxColor;
Colors(8).Color=DoubProxColor;
Colors(9).Color=BothColor;
Colors(10).Color=BothColor;
Colors(11).Color=DistalColor;
Colors(12).Color=ProxColor;
Colors(13).Color=BothSepColor;
Colors(14).Color=BothColor;
Colors(15).Color=DoubProxColor;
Colors(16).Color='k';
Colors(17).Color=DoubDistColor;
Colors(18).Color=DistalColor;
Colors(19).Color=DoubDistColor;
Colors(20).Color=BothColor;

fontsize=18;
fontname='Helvetica';
FigDirect=[DropboxFolder filesep 'Figures'];
%%
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for ee=1:length(AvgProdAllAP(cc).EmbryosProd)
        APBoxing(ee,cc)=AvgProdAllAP(cc).EmbryosProd(ee).MeanProd(APtoUse);
    end
end
APBoxing(APBoxing==0)=nan;
figure 
MaxFluoVal=0;
for cc=1:length(ConstructList)
    
    for ee=1:length(AvgProdAllAP(cc).EmbryosProd);
        yyaxis left
        plot(cc,AvgProdAllAP(cc).EmbryosProd(ee).MeanProd(APtoUse),'o','LineWidth',1.5,'Color',Colors(cc).Color);
        hold on   
        yyaxis right
        plot(cc, (AvgProdAllAP(cc).EmbryosProd(ee).MeanProd(APtoUse)/MS2Conversion),'o','LineWidth',1.5,'Color',Colors(cc).Color);
    if AvgProdAllAP(cc).EmbryosProd(ee).MeanProd(APtoUse)>MaxFluoVal
        MaxFluoVal=AvgProdAllAP(cc).EmbryosProd(ee).MeanProd(APtoUse);
    end
    end
end
yyaxis left
    boxplot(APBoxing,'Colors','k')
xlim([0 15]);
xlabel('Construct');
xticks([1:14]);
xticklabels({'Dist', 'Prox', 'Both Sep','1x Dist','1x Prox', '2x Dist', '2x Prox','sing 2xProx','Both','1x Both'});
ylabel('Total nascent mRNA (AU)');
title(['Mean mRNA production',' ' ,num2str(EgglengthUse),'% egg length'])
ylim([0 MaxFluoVal+1500]);
yyaxis right
ylabel('transcripts produced');
ylim([0 ((MaxFluoVal+1500)/MS2Conversion)])

%Do nc13 calc for Hb data

    for ee=1:length(AvgProdAllAP(16).EmbryosProd)
        APBoxing13(ee,1)=AvgProdAllAP(16).EmbryosProd(ee).MeanProd(13);
    end
APBoxing13(APBoxing13==0)=nan;
figure 
MaxFluoVal=0;
for ee=1:size(APBoxing13,2)
    yyaxis left
    plot(AvgProdAllAP(16).EmbryosProd(ee).MeanProd(13),'o','LineWidth',1.5,'Color','k');
    hold on 
    yyaxis right 
    plot((AvgProdAllAP(16).EmbryosProd(ee).MeanProd(13)/HbMS2Conversion),'o','LineWidth',1.5,'Color','k');
    if AvgProdAllAP(16).EmbryosProd(ee).MeanProd(13) >MaxFluoVal
        MaxFluoVal=AvgProdAllAP(16).EmbryosProd(ee).MeanProd(13);
    end
end
yyaxis left
ylabel('integrated fluorescence (AU)');
boxplot(APBoxing13,'Colors','k');
ylim([0 (MaxFluoVal+20000)]);
yyaxis right 
ylabel('transcripts produced');
ylim([0 ((MaxFluoVal+20000)/HbMS2Conversion)]);


%ANOVA at specific AP position of constructs against one another
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for bb=1:size([AvgProdAllAP(cc).AllProds],1)
    ConComp(bb,cc)=AvgProdAllAP(cc).AllProds(bb,APtoUse);
    end
    ConComp(ConComp==0)=nan;
end
[p,tbl,stats]=anova1(ConComp);
xlabel('Construct')
xticks([1:10]);
xticklabels({'Dist','Prox','BothSep','1x Dist','1x Prox','2xDist','2xProx','sing2xProx','Both','1xboth'});
ylabel('Total nascent mRNA');
title(['Avg mRNA production AP bin',' ',num2str(APtoUse)]);

%While we're here, is BothSep acting additively? 
ExpectedBoth=[([AvgProdAllAP(1).AvgProd(:)]./2) + ([AvgProdAllAP(2).AvgProd(:)]./2)];
ExpectedBothSE=[([AvgProdAllAP(1).ConSE(:)]./2) + ([AvgProdAllAP(2).ConSE(:)]./2)];
figure 
%errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).ConSE,'Color',Colors(1).Color);
hold on 
%errorbar(EggLength,AvgProdAllAP(2).AvgProd,AvgProdAllAP(2).ConSE,'Color',Colors(2).Color);
errorbar(EggLength,ExpectedBoth, ExpectedBothSE,'LineWidth',1.5);
errorbar(EggLength,AvgProdAllAP(3).AvgProd, AvgProdAllAP(3).ConSE,'Color',Colors(3).Color,'LineWidth',1.5);
xlabel('% egg length');
xlim([0 100])
ylabel('Avg total nascent mRNA');
title('Additivity of Both Sep');
legend('Expected', 'Both Sep');

% Additivity of Distals
ExpectedDoubDist=[([AvgProdAllAP(1).AvgProd(:)]) + ([AvgProdAllAP(1).AvgProd(:)])];
ExpectedDoubDistSE=AvgProdAllAP(1).ConSE(:);
figure 
hold on 
errorbar(EggLength,ExpectedDoubDist, ExpectedDoubDistSE,'Color',Colors(1).Color,'LineWidth',1.5);
errorbar(EggLength,AvgProdAllAP(17).AvgProd, AvgProdAllAP(17).ConSE,'Color',Colors(17).Color,'LineWidth',1.5);
xlabel('% egg length');
xlim([0 100])
ylabel('Avg total nascent mRNA');
title('Additivity of Distal');
legend('Expected', '2x Dist');

 %Additivity of Proximals
 ExpectedDoubProx=([AvgProdAllAP(2).AvgProd(:)]) + ([AvgProdAllAP(2).AvgProd(:)]);
 ExpectedDoubProxSE=AvgProdAllAP(2).ConSE(:);
 figure 
   hold on 
 errorbar(EggLength,ExpectedDoubProx, ExpectedDoubProxSE,'Color',Colors(2).Color,'LineWidth',1.5);
 errorbar(EggLength,AvgProdAllAP(7).AvgProd, AvgProdAllAP(7).ConSE,'Color',Colors(7).Color,'LineWidth',1.5);
 xlabel('% egg length');
 xlim([0 100])
 ylabel('Avg total nascent mRNA');
 title('Additivity of Proximal');
 
 % Plot mean total mRNA of each construct against AP position
 figure 
 for cc=1:length(ConstructList)-1
     plot(AvgProdAllAP(cc).AvgProd,'Color',Colors(cc).Color,'LineWidth',1.5);
     hold on 
 end
 %%
 %Compare single and doubles
 figure 
 errorbar(EggLength, AvgProdAllAP(1).AvgProd, AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength, AvgProdAllAP(17).AvgProd, AvgProdAllAP(17).All95Conf,'Color',Colors(17).Color,'LineWidth',2.5);
 legend('Distal', '2x Distal');
 xlabel('% Egg length')
 ylabel('integrated fluorescence')
 title('Avg mRNA production per nucleus')
 set(gca,'FontSize',fontsize,'FontName',fontname);
 ylim([0 700000]);
 xlim([0 100]);
 
 figure
 errorbar(EggLength,AvgProdAllAP(2).AvgProd,AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',1.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(7).AvgProd, AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',1.5);
 legend('Proximal',  '2x Proximal');
 xlabel('% Egg length')
 ylabel('integrated fluorescence')
 title('Avg mRNA production per nucleus') 
 ylim([0 700000]);
 xlim([0 100]);
 
 figure
 yyaxis left
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',1.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(2).AvgProd, AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',1.5);
 errorbar(EggLength,AvgProdAllAP(3).AvgProd, AvgProdAllAP(3).All95Conf,'Color',Colors(3).Color,'LineWidth',1.5);
 %legend('Distal','Proximal','Both Sep');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% Egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 700000]);
 xlim([0 100]);
 yyaxis right 
 errorbar(EggLength,(AvgProdAllAP(1).AvgProd./(MS2Conversion)),(AvgProdAllAP(1).All95Conf./(MS2Conversion)),'Color',Colors(1).Color,'LineWidth',1.5,'LineStyle','-');
 hold on 
 errorbar(EggLength,(AvgProdAllAP(2).AvgProd./(MS2Conversion)), (AvgProdAllAP(2).All95Conf./(MS2Conversion)),'Color',Colors(2).Color,'LineWidth',1.5,'LineStyle','-');
 errorbar(EggLength,(AvgProdAllAP(3).AvgProd./(MS2Conversion)), (AvgProdAllAP(3).All95Conf./(MS2Conversion)),'Color',Colors(3).Color,'LineWidth',1.5,'LineStyle','-');
 ylabel('Total mRNA molecules');
 ylim([0 (700000/MS2Conversion)]);
 
  figure
  %yyaxis left
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(2).AvgProd, AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
 errorbar(EggLength,AvgProdAllAP(9).AvgProd, AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 yyaxis right
 plot(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion), 'Color', Colors(1).Color);
 plot(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion), 'Color', Colors(2).Color);
 plot(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion), 'Color', Colors(9).Color);
 set(gca, 'FontSize', fontsize, 'FontName', fontname,'YColor','k');
 xlabel('% egg length')
 ylabel('total transcripts produced')
 %title('Avg mRNA production per nucleus') 
 ylim([0 900000./MS2Conversion]);
 xlim([0 100]);
 yyaxis left 
 ylabel('integrated fluorescence (AU)');
 ylim([0 900000]);
print( [FigDirect filesep 'SinglesvBothTotalmRNA'],'-dsvg');
 
  figure
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(2).AvgProd, AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
 %legend('Distal','Proximal','Both Sep');
 yyaxis right
 plot(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion), 'Color', Colors(1).Color);
 plot(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion), 'Color', Colors(2).Color);
 ylabel('total transcripts produced');
 set(gca, 'FontSize', fontsize, 'FontName', fontname,'YColor','k');
 ylim([0 900000./MS2Conversion]);
 xlim([0 100]);
 yyaxis left
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 900000]);
 print( '-painters',[FigDirect filesep 'DistProxTotalmRNA'],'-dsvg');

 %% single allele 
  figure
 errorbar(EggLength,AvgProdAllAP(3).AvgProd,AvgProdAllAP(3).All95Conf,'Color',Colors(3).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(4).AvgProd, AvgProdAllAP(4).All95Conf,'Color',Colors(4).Color,'LineWidth',2.5,'LineStyle',':');
  errorbar(EggLength,AvgProdAllAP(5).AvgProd, AvgProdAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
 legend('Both Sep','1x Distal','1x Proximal');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1100000]);
 xlim([0 100]);
 
 figure
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(4).AvgProd, AvgProdAllAP(4).All95Conf,'Color',Colors(4).Color,'LineWidth',2.5,'LineStyle',':');
 %legend('Distal','1x Distal');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1100000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion),(AvgProdAllAP(1).All95Conf./MS2Conversion),'Color',Colors(1).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(4).AvgProd./MS2Conversion),(AvgProdAllAP(4).All95Conf./MS2Conversion),'Color',Colors(4).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1100000/MS2Conversion)]);
set(gca, 'YColor','k');
print( '-painters',[FigDirect filesep 'HemiDistTotalmRNA'],'-dsvg');
 
  figure
 errorbar(EggLength,AvgProdAllAP(2).AvgProd,AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(5).AvgProd, AvgProdAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
 %legend('Proximal','1x Proximal');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1100000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(2).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(5).AvgProd./MS2Conversion),(AvgProdAllAP(5).All95Conf./MS2Conversion),'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1100000/MS2Conversion)]);
set(gca,'YColor','k');
print( '-painters',[FigDirect filesep 'HemiProxTotalmRNA'],'-dsvg');
 
  figure
 errorbar(EggLength,AvgProdAllAP(4).AvgProd,AvgProdAllAP(4).All95Conf,'Color',Colors(4).Color,'LineWidth',2.5,'LineStyle',':');
 hold on 
 errorbar(EggLength,AvgProdAllAP(9).AvgProd, AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 legend('1x Distal','Both');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1100000]);
 xlim([0 100]);
 
 figure
 errorbar(EggLength,AvgProdAllAP(5).AvgProd,AvgProdAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
 hold on 
 errorbar(EggLength,AvgProdAllAP(9).AvgProd, AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 legend('1x Proximal','Both');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1100000]);
 xlim([0 100]);
 
  figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(10).AvgProd, AvgProdAllAP(10).All95Conf,'Color',Colors(10).Color,'LineWidth',2.5,'LineStyle',':');
 %legend('Both','1x Both');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1100000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(10).AvgProd./MS2Conversion),(AvgProdAllAP(10).All95Conf./MS2Conversion),'Color',Colors(10).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1100000/MS2Conversion)]);
set(gca,'YColor','k');
print( '-painters',[FigDirect filesep 'HemiBothTotalmRNA'],'-dsvg');
 
  figure
 errorbar(EggLength,AvgProdAllAP(7).AvgProd,AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(8).AvgProd,AvgProdAllAP(8).All95Conf,'Color',Colors(8).Color,'LineWidth',2.5,'LineStyle',':');
 %legend('2xProximal','Single 2xProximal');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 ylim([0 1100000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(8).AvgProd./MS2Conversion),(AvgProdAllAP(8).All95Conf./MS2Conversion),'Color',Colors(8).Color,'LineWidth',2.5,'LineStyle',':');
ylabel('total transcripts produced');
ylim([0 (1100000/MS2Conversion)]);
set(gca,'YColor','k');
print( '-painters',[FigDirect filesep 'Hemi2xProxTotalmRNA'],'-dsvg');
 
  figure
 errorbar(EggLength,AvgProdAllAP(5).AvgProd,AvgProdAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
 hold on 
 errorbar(EggLength,AvgProdAllAP(8).AvgProd,AvgProdAllAP(8).All95Conf,'Color',Colors(8).Color,'LineWidth',2.5,'LineStyle',':');
 legend('1x Proximal','Single 2xProximal');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 ylim([0 1100000]);
 xlim([0 100]);
 
   figure
 errorbar(EggLength,AvgProdAllAP(4).AvgProd,AvgProdAllAP(4).All95Conf,'Color',Colors(4).Color,'LineWidth',2.5,'LineStyle',':');
 hold on 
 errorbar(EggLength,AvgProdAllAP(5).AvgProd,AvgProdAllAP(5).All95Conf,'Color',Colors(5).Color,'LineWidth',2.5,'LineStyle',':');
 errorbar(EggLength,AvgProdAllAP(10).AvgProd, AvgProdAllAP(10).All95Conf,'Color',Colors(10).Color,'LineWidth',2.5,'LineStyle',':');
 %legend('1x Distal','1x Proximal','1x Both');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 ylim([0 1100000]);
 xlim([0 100]);

 %%
 %Compare SE and duplicates
 figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(7).AvgProd, AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
errorbar(EggLength, AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 xlim([0 100]);
 ylim([0 900000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
ylabel('transcripts produced');
ylim([0 (900000/MS2Conversion)]);
set(gca,'FontSize',fontsize,'FontName',fontname);
 print( [FigDirect filesep '2ProxBothTotalmRNA'],'-dsvg');
%legend('2x Proximal', 'Both','Location','best');
 %title('Avg mRNA production per nucleus') 
 
  figure 
 hold on 
 errorbar(EggLength, AvgProdAllAP(17).AvgProd, AvgProdAllAP(17).All95Conf,'Color',Colors(17).Color,'LineWidth',2.5);
errorbar(EggLength, AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 xlim([0 100]);
 ylim([0 1500000]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(17).AvgProd./MS2Conversion),(AvgProdAllAP(17).All95Conf./MS2Conversion),'Color',Colors(17).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
ylabel('transcripts produced');
ylim([0 (1500000/MS2Conversion)]);
set(gca,'FontSize',fontsize,'FontName',fontname,'YColor','k');
 print( [FigDirect filesep '2DistBothTotalmRNA'],'-dsvg');
 %%
 %32C vs RT
 figure
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(11).AvgProd, AvgProdAllAP(11).All95Conf,'Color',Colors(11).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion),(AvgProdAllAP(1).All95Conf./MS2Conversion),'Color',Colors(1).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(11).AvgProd./MS2Conversion),(AvgProdAllAP(11).All95Conf./MS2Conversion),'Color',Colors(11).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (950000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep 'DistTempCompTotalmRNA'],'-dsvg');


figure
 errorbar(EggLength,AvgProdAllAP(2).AvgProd,AvgProdAllAP(2).All95Conf,'Color',Colors(2).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(12).AvgProd, AvgProdAllAP(12).All95Conf,'Color',Colors(12).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(2).AvgProd./MS2Conversion),(AvgProdAllAP(2).All95Conf./MS2Conversion),'Color',Colors(2).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(12).AvgProd./MS2Conversion),(AvgProdAllAP(12).All95Conf./MS2Conversion),'Color',Colors(12).Color,'LineWidth',2.5,'LineStyle','--');
ylabel('transcripts produced');
ylim([0 (950000/MS2Conversion)]);
%legend('Proximal', 'Proximal 32C','Location','best');
 print( [FigDirect filesep 'ProxTempTotalmRNA'],'-dsvg');


 
 figure
 errorbar(EggLength,AvgProdAllAP(3).AvgProd,AvgProdAllAP(3).All95Conf,'Color',Colors(3).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(13).AvgProd, AvgProdAllAP(12).All95Conf,'Color',Colors(13).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(3).AvgProd./MS2Conversion),(AvgProdAllAP(3).All95Conf./MS2Conversion),'Color',Colors(3).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(13).AvgProd./MS2Conversion),(AvgProdAllAP(12).All95Conf./MS2Conversion),'Color',Colors(13).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (950000/MS2Conversion)]);
%legend('Both Sep', 'Both Sep 32C','Location','best');
 print( [FigDirect filesep 'SepTempTotalmRNA'],'-dsvg');

 
figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(14).AvgProd, AvgProdAllAP(13).All95Conf,'Color',Colors(14).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(14).AvgProd./MS2Conversion),(AvgProdAllAP(14).All95Conf./MS2Conversion),'Color',Colors(14).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (950000/MS2Conversion)]);
%legend('Both', 'Both 32C','Location','best');
print( [FigDirect filesep 'BothTempTotalmRNA'],'-dsvg');

 
 figure
 errorbar(EggLength,AvgProdAllAP(7).AvgProd,AvgProdAllAP(7).All95Conf,'Color',Colors(7).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(15).AvgProd, AvgProdAllAP(15).All95Conf,'Color',Colors(15).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence (AU)')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(7).AvgProd./MS2Conversion),(AvgProdAllAP(7).All95Conf./MS2Conversion),'Color',Colors(7).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(15).AvgProd./MS2Conversion),(AvgProdAllAP(15).All95Conf./MS2Conversion),'Color',Colors(15).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (950000/MS2Conversion)]);
%legend('2x Prox', '2x Prox 32C','Location','best');
 print( [FigDirect filesep '2xProxTempTotalmRNA'],'-dsvg');
 
 figure
 errorbar(EggLength,AvgProdAllAP(6).AvgProd,AvgProdAllAP(6).All95Conf,'Color',Colors(6).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(19).AvgProd, AvgProdAllAP(19).All95Conf,'Color',Colors(19).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1150000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(6).AvgProd./MS2Conversion),(AvgProdAllAP(6).All95Conf./MS2Conversion),'Color',Colors(6).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(19).AvgProd./MS2Conversion),(AvgProdAllAP(19).All95Conf./MS2Conversion),'Color',Colors(19).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (1150000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep '2xDistTempCompTotalmRNA'],'-dsvg');
 
 %% 17C vs RT
 figure
 errorbar(EggLength,AvgProdAllAP(1).AvgProd,AvgProdAllAP(1).All95Conf,'Color',Colors(1).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(18).AvgProd, AvgProdAllAP(18).All95Conf,'Color',Colors(18).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 950000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(1).AvgProd./MS2Conversion),(AvgProdAllAP(1).All95Conf./MS2Conversion),'Color',Colors(1).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(18).AvgProd./MS2Conversion),(AvgProdAllAP(18).All95Conf./MS2Conversion),'Color',Colors(18).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced','Color','k');
set(gca,'ycolor','k');
ylim([0 (950000/MS2Conversion)]);
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep 'DistC17CompTotalmRNA'],'-dsvg');

figure
 errorbar(EggLength,AvgProdAllAP(9).AvgProd,AvgProdAllAP(9).All95Conf,'Color',Colors(9).Color,'LineWidth',2.5);
 hold on 
 errorbar(EggLength,AvgProdAllAP(20).AvgProd, AvgProdAllAP(20).All95Conf,'Color',Colors(20).Color,'LineWidth',2.5,'LineStyle','-.');
 set(gca, 'FontSize', fontsize, 'FontName', fontname);
 xlabel('% egg length')
 ylabel('integrated fluorescence')
 %title('Avg mRNA production per nucleus') 
 ylim([0 1100000]);
 xlim([0 100]);
 yyaxis right
 errorbar(EggLength, (AvgProdAllAP(9).AvgProd./MS2Conversion),(AvgProdAllAP(9).All95Conf./MS2Conversion),'Color',Colors(9).Color,'LineWidth',2.5,'LineStyle','-');
  errorbar(EggLength, (AvgProdAllAP(20).AvgProd./MS2Conversion),(AvgProdAllAP(20).All95Conf./MS2Conversion),'Color',Colors(20).Color,'LineWidth',2.5,'LineStyle','-.');
ylabel('transcripts produced');
ylim([0 (1100000/MS2Conversion)]);
set(gca, 'ycolor','k')
%legend('Distal', 'Distal 32C','Location','best');
print( [FigDirect filesep 'BothC17CompTotalmRNA'],'-dsvg');
 %%

 
 %Calculate area under curve for saying where expression %'s are 
 TrapzSE=[AvgProdAllAP(8).AvgProd]
 TrapzSE=TrapzSE(~isnan(TrapzSE));
 TotalAreaSE=trapz(TrapzSE);
 Area40to60EL=trapz(AvgProdAllAP(8).AvgProd(19:23));
 Percent40t060=Area40to60EL/TotalAreaSE; 
 
% Save the mRNA production info in the Constructs folder
save([DropboxFolder filesep 'Constructs' filesep 'AllTotalmRNAProd'],'AvgProdAllAP');
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

TotalProdManovaVect=[];
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
            if isempty(BurstProperties(nn).TotalmRNA)
                TotalProdManovaVect=[TotalProdManovaVect, nan];
            else
            TotalProdManovaVect=[TotalProdManovaVect, BurstProperties(nn).TotalmRNA];
            end
            APManovaVect=[APManovaVect, BurstProperties(nn).APBin];
            ConManovaVect=[ConManovaVect, cc];
        end
    end
end

[p,tbl,stats]=anovan(TotalProdManovaVect,{APManovaVect, ConManovaVect},'sstype',1,'varnames',{'AP bin','Construct'})
%multcompare(stats,'Dimension',2)






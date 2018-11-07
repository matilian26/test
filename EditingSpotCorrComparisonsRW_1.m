CheckAll=input('Run CompareSpots for all Data?','s')
if CheckAll=='y'
ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrDist32C';'KrBothSep32C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

% %Ask if only want nc14 info
% ncUse=input('Want to only use nc14?','s');

%First make sure have done SmoothingStuffnc14RW for all embryos
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        
        firstfilename=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
        load(firstfilename);
        secondfilename= [DropboxFolder filesep PrefixName filesep 'APDetection.mat'];
        load (secondfilename);
        thirdfilename= [DropboxFolder filesep PrefixName filesep PrefixName '_lin.mat'];
        load (thirdfilename);
        
            run('ComparingSpotsRWAdj.m')
            save([DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj'],'SpotDiff')
        %else
            %run SmoothingStuffRW.m
        %end
    end
    %clear BurstProperties
end
end


%%
%% %Actually doing the comparison
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDist32C';'KrProx32C';'KrBothSep32C';'Kr2xDist32C';'Kr2xProx32C';'KrBoth32C'; 'KrDist17C';'KrBoth17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

ConCorrMean=[];
MeanAPCorr=[];
MeanAllCorr=[];
APCorr=[];

for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgCorrCon=[];
    firsttime=1;
    ConCorrAllAP=[];
    ConCorrSE=[];
    ConCorrSD=[];
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        
        filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        load(filename);
        
        %seperate out by AP bin
        CorrAllAP=[];
        for aa=2:length(APbinID) %skip APbin 0
            CorrAP=[];
            CorrelationAP=find([SpotDiff.APBin]==APbinID(aa));
            if isempty(CorrelationAP)
                CorrAP=[CorrAP; nan];
            else
            for bb=1:length(CorrelationAP)  
                if ~isfield(SpotDiff, 'SpotCorr')
                    %break test out if can just set as nan if no corr spots
                    CorrAP=[CorrAP; nan];
                elseif isempty(SpotDiff(CorrelationAP(bb)).SpotCorr)
                    CorrAP=[CorrAP;nan];
                elseif length([SpotDiff(CorrelationAP(bb)).SpotCorr])==1
                    CorrAP=[CorrAP;SpotDiff(CorrelationAP(bb)).SpotCorr(1)];
                else
                CorrAP=[CorrAP;SpotDiff(CorrelationAP(bb)).SpotCorr(1,2)];  %put all durations at a given AP value in a column going down
                end
            end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(CorrAP)
                CorrAllAP(bb,aa,ee)=CorrAP(bb);
            end
            CorrAllAP(CorrAllAP==0)=nan;
            
            CorrelationSD(ee,aa,cc)=nanstd(CorrAllAP(:,aa,ee));
        %CorrelationSE(ee,aa,cc)=CorrelationSD(ee,aa,cc)/sqrt(length(CorrelationAP));
        CorrelationSE(ee,aa,cc)=CorrelationSD(ee,aa,cc)/sqrt(sum(~isnan(CorrAP)));
            clear CorrelationAP
        end
        %Compile all correlation data for a construct in one long column
        %to put in structure
        for bb=1:size(CorrAllAP,3)
            ConCorrAllAP=[ConCorrAllAP; CorrAllAP(:,:,bb)];
        end
        for bb=1:size(CorrelationSD,3)
            ConCorrSD=[ConCorrSD; CorrelationSD(:,:,bb)];
        end
        for bb=1:size(CorrelationSE,3)
            ConCorrSE=[ConCorrSE;CorrelationSE(:,:,bb)];
        end
        AvgCorrAllAP(cc).EmbryoCorr(ee).AvgCorr=nanmean(CorrAllAP(:,:,ee));
        AvgCorrAllAP(cc).EmbryoCorr(ee).SE=CorrelationSE(ee,:,cc);
        AvgCorrAllAP(cc).EmbryoCorr(ee).SD=CorrelationSD(ee,:,cc);
    end
        AvgCorrAllAP(cc).AvgCorr=nanmean(ConCorrAllAP);  %Avg duration of all embryos of a construct
        AvgCorrAllAP(cc).CorrSD=nanmean(ConCorrSD);
        AvgCorrAllAP(cc).CorrSE=nanmean(ConCorrSE);
        AvgCorrAllAP(cc).AllCorrs=[ConCorrAllAP];
        AvgCorrAllAP(cc).AllSD=nanstd(AvgCorrAllAP(cc).AllCorrs);
        %AvgCorrAllAP(cc).AllSE=((AvgCorrAllAP(cc).AllSD)/sqrt(length(AvgCorrAllAP(cc).AllCorrs)));
        for aa=1:length(APbinID)
        AvgCorrAllAP(cc).AllSE(aa)=(AvgCorrAllAP(cc).AllSD(aa))/sqrt(sum(~isnan(AvgCorrAllAP(cc).AllCorrs(:,aa))));
        AvgCorrAllAP(cc).NumNucUsed(aa)=sum(~isnan(AvgCorrAllAP(cc).AllCorrs(:,aa)));
        end
        AvgCorrAllAP(cc).All95Conf=(AvgCorrAllAP(cc).AllSE) .* 1.95;

end

%Remove data relying on only one data point for an AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
    if sum(~isnan(AvgCorrAllAP(cc).AllCorrs(:,aa)))==1
        AvgCorrAllAP(cc).AllSD(aa)=nan;
        AvgCorrAllAP(cc).AllSE(aa)=nan;
        AvgCorrAllAP(cc).All95Conf(aa)=nan;
        AvgCorrAllAP(cc).AvgCorr(aa)=nan;
    end
    end
end

%% Plotting 
%Plot mean of each embryo for each construct
DistalColor=[1 64 172]./255;
Distal32CColor=[118 180 238] ./ 255;
DistalEmptyColor=[8 210 238] ./ 255; 
DoubDistColor=[73 184 253] ./ 255;
ProxColor=[238 123 23]./255;
ProxEmptyColor=[251 250 50] ./255;
Proximal32CColor=[251 150 10] ./ 255;
DoubProxColor=[215 183 58] ./ 255;
DoubProxEmptyColor=[251 220 50] ./ 255;
BothSepColor=[149 188 114] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[52 119 71]./255;
Both32CColor=[149 188 114] ./ 255;
BothEmptyColor=[12 250 100] ./ 255;
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
Colors(10).Color=DoubDistColor;
Colors(11).Color=DoubProxColor;
Colors(12).Color=BothColor;
Colors(13).Color=DistalColor;
Colors(14).Color=BothColor;

fontsize=18;
fontname='Helvetica';
FigDirect=[DropboxFolder filesep 'Figures'];

EggLength=APbinID .* 100;
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;

%% 
for cc=1:length(ConstructList)
    for ee=1:length(AvgCorrAllAP(cc).EmbryoCorr)
        APBoxing(ee,cc)=AvgCorrAllAP(cc).EmbryoCorr(ee).AvgCorr(APtoUse);
    end
end
APBoxing(APBoxing==0)=nan;
figure 

    for ee=1:length(AvgCorrAllAP(1).EmbryoCorr)
    plot(1,AvgCorrAllAP(1).EmbryoCorr(ee).AvgCorr(APtoUse),'o','LineWidth',1.5,'Color',DistalColor)
    hold on
    end
    for ee=1:length(AvgCorrAllAP(2).EmbryoCorr)
    plot(2,AvgCorrAllAP(2).EmbryoCorr(ee).AvgCorr(APtoUse),'o','LineWidth',1.5,'Color',ProxColor)
        
    hold on
    end
    
    for ee=1:length(AvgCorrAllAP(3).EmbryoCorr)
    plot(3,AvgCorrAllAP(3).EmbryoCorr(ee).AvgCorr(APtoUse),'o','LineWidth',1.5,'Color',BothSepColor)
    hold on
    end
    for ee=1:length(AvgCorrAllAP(4).EmbryoCorr)
    plot(4,AvgCorrAllAP(4).EmbryoCorr(ee).AvgCorr(APtoUse),'o','LineWidth',1.5,'Color',DoubDistColor)
    hold on
    end
    for ee=1:length(AvgCorrAllAP(5).EmbryoCorr)
    plot(5,AvgCorrAllAP(5).EmbryoCorr(ee).AvgCorr(APtoUse),'o','LineWidth',1.5,'Color',DoubProxColor)
    hold on
    end
    for ee=1:length(AvgCorrAllAP(6).EmbryoCorr)
    plot(6,AvgCorrAllAP(6).EmbryoCorr(ee).AvgCorr(APtoUse),'o','LineWidth',1.5,'Color',BothColor);
    %errorbar(6,AvgCorrAllAP(6).EmbryoCorr(ee).MeanProd(APtoUse),AvgCorrAllAP(6).EmbryoCorr(ee).SE(APtoUse),'o','LineWidth',1.5,'Color',BothColor)
    hold on
    end
    boxplot(APBoxing,'Colors','k')
xlim([0 7]);
xlabel('Construct');
xticks([1:6]);
xticklabels({'Dist', 'Prox', 'Separated', '2x Dist', '2x Prox','Both'});
ylabel('Correlation');
title(['Mean correlation of allele activity',' ' ,num2str(EgglengthUse),'% egg length'])

%Plot mean correlation vs AP position for each construct 
figure 
for cc=1:length(ConstructList)
    plot(EggLength, AvgCorrAllAP(cc).AvgCorr,'Color',Colors(cc).Color, 'LineWidth', 1.5);
    hold on 
end
xlabel('% Egg length');
ylabel('Correlation');
title('Mean correlation of allele activity');
legend('Distal', 'Proximal', 'Separated', '2x Distal', '2x Proximal', 'Both');
%%
%Singles vs Both sep
figure 
errorbar(EggLength, AvgCorrAllAP(1).AvgCorr, AvgCorrAllAP(1).All95Conf, 'Color',Colors(1).Color, 'LineWidth',2.5); 
hold on 
errorbar(EggLength, AvgCorrAllAP(2).AvgCorr, AvgCorrAllAP(2).All95Conf, 'Color',Colors(2).Color, 'LineWidth',2.5); 
errorbar(EggLength, AvgCorrAllAP(3).AvgCorr, AvgCorrAllAP(3).All95Conf, 'Color',Colors(3).Color, 'LineWidth',2.5); 
%legend('Distal', 'Proximal', 'Separated');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('allele correlation');
xlim([0 100]);
print( [FigDirect filesep 'SinglesvSepCorrelation'],'-dsvg');

%% Duplicated enhancers
figure 
errorbar(EggLength, AvgCorrAllAP(4).AvgCorr, AvgCorrAllAP(4).All95Conf, 'Color',Colors(4).Color, 'LineWidth',2.5); 
hold on 
errorbar(EggLength, AvgCorrAllAP(5).AvgCorr, AvgCorrAllAP(5).All95Conf, 'Color',Colors(5).Color, 'LineWidth',2.5); 
errorbar(EggLength, AvgCorrAllAP(6).AvgCorr, AvgCorrAllAP(6).All95Conf, 'Color',Colors(6).Color, 'LineWidth',2.5); 
%legend('Distal', 'Proximal', 'Separated');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('allele correlation');
xlim([0 100]);


%%
%Temperature elevation 32C
figure 
errorbar(EggLength, AvgCorrAllAP(1).AvgCorr, AvgCorrAllAP(1).All95Conf, 'Color',Colors(1).Color, 'LineWidth',2.5); 
hold on 
errorbar(EggLength, AvgCorrAllAP(7).AvgCorr, AvgCorrAllAP(7).All95Conf, 'Color',Colors(7).Color, 'LineWidth',2.5,'LineStyle','-.'); 
legend('Distal', 'Distal 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Allele correlation');
xlim([0 100]);

figure 
errorbar(EggLength, AvgCorrAllAP(2).AvgCorr, AvgCorrAllAP(2).All95Conf, 'Color',Colors(2).Color, 'LineWidth',2.5); 
hold on 
errorbar(EggLength, AvgCorrAllAP(8).AvgCorr, AvgCorrAllAP(8).All95Conf, 'Color',Colors(8).Color, 'LineWidth',2.5,'LineStyle','-.'); 
%legend('Distal', 'Distal 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Allele correlation');
xlim([0 100]);

figure 
errorbar(EggLength, AvgCorrAllAP(3).AvgCorr, AvgCorrAllAP(3).All95Conf, 'Color',Colors(3).Color, 'LineWidth',2.5); 
hold on 
errorbar(EggLength, AvgCorrAllAP(9).AvgCorr, AvgCorrAllAP(9).All95Conf, 'Color',Colors(9).Color, 'LineWidth',2.5,'LineStyle','-.'); 
%legend('Both Sep', 'Both Sep 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Allele correlation');
xlim([0 100]);

figure 
errorbar(EggLength, AvgCorrAllAP(7).AvgCorr, AvgCorrAllAP(7).All95Conf, 'Color',Colors(7).Color, 'LineWidth',2.5,'LineStyle','--'); 
hold on 
errorbar(EggLength, AvgCorrAllAP(8).AvgCorr, AvgCorrAllAP(8).All95Conf, 'Color',Colors(8).Color, 'LineWidth',2.5,'LineStyle','--'); 
errorbar(EggLength, AvgCorrAllAP(9).AvgCorr, AvgCorrAllAP(9).All95Conf, 'Color',Colors(9).Color, 'LineWidth',2.5,'LineStyle','--'); 
%legend('Both Sep', 'Both Sep 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% Egg length');
ylabel('Allele correlation');
xlim([0 100]);

figure 
errorbar(EggLength, AvgCorrAllAP(10).AvgCorr, AvgCorrAllAP(10).All95Conf, 'Color',Colors(10).Color, 'LineWidth',2.5,'LineStyle','--'); 
hold on 
errorbar(EggLength, AvgCorrAllAP(11).AvgCorr, AvgCorrAllAP(11).All95Conf, 'Color',Colors(11).Color, 'LineWidth',2.5,'LineStyle','--'); 
errorbar(EggLength, AvgCorrAllAP(9).AvgCorr, AvgCorrAllAP(9).All95Conf, 'Color',Colors(9).Color, 'LineWidth',2.5,'LineStyle','--'); 
%legend('Distal', 'Proximal', 'Separated');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('allele correlation');
xlim([0 100]);
%%
%title('Mean correlation of allele activity');
%Plot bar graph of avg correlation at designated AP bin for each construct
APtoUse=input('Which AP bin to use?');
EggLengthAPbin=APbinID(APtoUse)*100;
figure
for cc=1:length(ConstructList)
    MeanCorratAP(1,cc)=AvgCorrAllAP(cc).AvgCorr(APtoUse);
    errorbar(cc,MeanCorratAP(1,cc), AvgCorrAllAP(cc).CorrSE(APtoUse),'o');
    hold on 
end
bar(MeanCorratAP,'EdgeColor','k','LineWidth',1.5)
xlabel('Construct');
xticks([1:6]);
xticklabels({'Dist','Prox','Separated','2xDist','2xProx','Both'});
ylabel('correlation');
title(['Mean correlation of allele activity',' ',num2str(EggLength),'% egg length']);

%Do ANOVA to compare constructs at designated AP bin
for cc=1:length(ConstructList)
    for bb=1:length([AvgCorrAllAP(cc).AllCorrs])
    NCorrComp(bb,cc)=AvgCorrAllAP(cc).AllCorrs(bb,APtoUse);
    end
    NCorrComp(NCorrComp==0)=nan;
end
[p,tbl,stats]=anova1(NCorrComp);

%% 2way anova of Duration vs AP position vs genotype 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info
%ncUse=input('Want to only use nc14?','s');

%Count for each construct

CorrManovaVect=[];
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
%         if ncUse=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
%         else
%             filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelation.mat'];
%         end 
        load(filename);
        for nn=1:length(SpotDiff)
             if isempty(SpotDiff(nn).SpotCorr)
                 CorrManovaVect=[CorrManovaVect, nan];
             elseif length(SpotDiff(nn).SpotCorr)==1
                 CorrManovaVect=[CorrManovaVect, SpotDiff(nn).SpotCorr(1)];
             else
             CorrManovaVect=[CorrManovaVect, SpotDiff(nn).SpotCorr(1,2)];
             end
            APManovaVect=[APManovaVect, SpotDiff(nn).APBin];
            ConManovaVect=[ConManovaVect, cc];
        end
    end
end

[p,tbl,stats]=anovan(CorrManovaVect,{APManovaVect, ConManovaVect},'sstype',1,'varnames',{'AP bin','Construct'})
%multcompare(stats,'Dimension',2)

%Do ANOVA within each construct to look at variation across AP bin 

% for cc=1:1 %length(ConstructList)
%     [p,tbl,statzanova1(AvgCorrAllAP(cc).AllCorrs);
%     xlabel('AP bin');
%     ylabel('Spot correlation');
%     title(ConstructList{cc});
%     end
%   
% end 

% Plot mean of each embryo at certain AP position for all constructs
% figure
% for cc=1:length(ConstructList)
% end
% %%
% for ii=1:length(ConstructList)
%     Data= LoadMS2SetsCS(ConstructList{ii});
%     
%     NEmbryos = length(Data);
%     Label = ConstructList(ii);
%     figure
%     APCorr=[];
%     
%     for jj=1:NEmbryos
%         PrefixName=Data(jj).Prefix;
%         filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
%         load(filename);
%         Bins=unique([SpotDiff.APBin]);%want to find mean correlation of spots at each APbin for fewer points/easier to see graph when combining
%         APBins=Data(1).APbinID;
%         APCorr=nan(1,length(APBins));
% for zz=1:length(APBins)
%     %temp3=find([APBins==Bins(zz)]);
%     temp2=find([SpotDiff.APBin]==APBins(zz));%issue with find I think has to do with empty values 
%     for qq=1:length(temp2)
%         if ~isempty(SpotDiff(temp2(qq)).SpotCorr)
%             if length(SpotDiff(temp2(qq)).SpotCorr)==1
%                 APCorr(qq,zz,jj)=[SpotDiff(temp2(qq)).SpotCorr(1)]; %rows are spots, columns APbin, page embryos
%             else
%     APCorr(qq,zz,jj)=[SpotDiff(temp2(qq)).SpotCorr(1,2)];
%             end
%         end
%     end 
% end
% APCorr(APCorr==0)=nan;
% 
% MeanAPCorr(jj,:,ii)=nanmean(APCorr(:,:,jj));%mean of each embryo is a row and 3rd dimension is construct
% MeanAllCorr(jj,ii)=nanmean(MeanAPCorr(jj,:,ii),2);
% %setting up a structure to compare correlation and AP position (within each
% %construct- working elsewhere on 2way ANOVA)
% 
% counter=0;
% for tt=1:length(SpotDiff)
%     APColumn=find([APBins==SpotDiff(tt).APBin]);
%     if ~isempty([SpotDiff(tt).SpotCorr])
%         counter=counter+1;
%     if length(SpotDiff(tt).SpotCorr)==1
%         SetupAP(counter,APColumn,jj)=SpotDiff(tt).SpotCorr(1)
%     else
%     SetupAP(counter,APColumn,jj)=SpotDiff(tt).SpotCorr(1,2);
%     end
%     %SetupAP(counter,APColumn,jj)=SpotDiff(tt).APBin;
%     end
% end
% 
% 
%     end
%     
%     ConstAPCor=[];
%     for yy=1:NEmbryos
%         ConstAPCor=vertcat(ConstAPCor, SetupAP(:,:,yy)); %for all embryos of a construct combine so columns are AP bins and going down are all correlation values 
%     end
% ConstAPCor(ConstAPCor==0)=nan;
% CorrValues(ii).AllCorrsbyAP=ConstAPCor;
% 
% end
% %%
% CompAPCorrs=[];
% for ii=1:length(ConstructList)
%     if ii==1
%     CompAPCorrs=[CompAPCorrs CorrValues(ii).AllCorrsbyAP(:,21)];
%     else
%         if length([CorrValues(ii).AllCorrsbyAP]) > length(CompAPCorrs)
%             CompAPCorrs(end:length([CorrValues(ii).AllCorrsbyAP]))=nan;
%             CompAPCorrs=[CompAPCorrs CorrValues(ii).AllCorrsbyAP(:,21)];
%         end
%         if length(CorrValues(ii).AllCorrsbyAP) < length(CompAPCorrs)
%             CorrValues(ii).AllCorrsbyAP(length(CompAPCorrs),1)=nan;
%             CompAPCorrs=[CompAPCorrs CorrValues(ii).AllCorrsbyAP(:,21)];
%         end
%         if length(CorrValues(ii).AllCorrsbyAP) == length(CompAPCorrs)
%             CompAPCorrs=[CompAPCorrs CorrValues(ii).AllCorrsbyAP(:,21)];
%         end
%     end
% end
% %%
%     
%     AllCorr(ii).AllValues=[APCorr(:)];
%     %AllCorr(:,ii)=vertcat(AllCorr,[APCorr(:)]);
% 
%     ConstructAPCorr=[];
%     for tt=1:NEmbryos
%     ConstructAPCorr=vertcat(ConstructAPCorr,APCorr(:,:,tt));
%     end
%     %
%     ConCorrMean(ii,:)=nanmean(ConstructAPCorr,1);
%     
%     h(ii)=plot(APBins, ConCorrMean(ii,:), 'LineWidth',2);
%     hold on
%     xlabel('Mean AP Position')
%     ylabel('Average fluorescence correlation')
%     Prefixname=ConstructList(ii)
%     %legend([h(ii)],Prefixname);
%     title(['2 Spot Correlation',Prefixname]);
%     ylim([0 1]);
%     xlim([0.4 0.75])
%     
%      
%     
%     clearvars -except APStats ConstructList DropboxFolder ConCorrMean AllCorr  APBins MeanAPCorr MeanAllCorr
% %end
% for jj=1:length(ConstructList)
%     NNucs(jj)=length(AllCorr(jj).AllValues);
%     MaxNucs=max(NNucs);
% end
% AllCorrMat=zeros(MaxNucs,4);
% for jj=1:length(ConstructList)
%     for qq=1:length(AllCorr(jj).AllValues)
%     AllCorrMat(qq,jj)=AllCorr(jj).AllValues(qq);
%     end
% end
% AllCorrMat(AllCorrMat==0)=nan;  %ANOVA of spot correlation at all AP points across whole embryo for all embryos of a construct with columns being each different construct
% [p,tbl,stats]=anova1(AllCorrMat);
% title('Spot Correlation across whole embryo')
% xlabel('construct')
% figure
% multcompare(stats);
% 
% 
% 
% 
% %MeanAPCorr(MeanAPCorr==0)=nan;
% % MeanAllCorr(MeanAllCorr==0)=nan;    %ANOVA of mean spot correlation across whole embryo with groups being the different constructs 
% % Constructs=[ConstructList{1},'', ConstructList{2},'', ConstructList{3},'',ConstructList{4}];
% % [p,tbl,stats]=anova1(MeanAllCorr)
% % title('Mean Spot Correlation across whole embryo')
% % xlabel('construct')
% % figure
% % multcompare(stats);
% 
% BothData=LoadMS2SetsCS(ConstructList{1});
% DistData=LoadMS2SetsCS(ConstructList{2});
% ProxData=LoadMS2SetsCS(ConstructList{3});
% BothSepData=LoadMS2SetsCS(ConstructList{4});
% DoubDistData=LoadMS2SetsCS(ConstructList{5});
% DoubProxData=LoadMS2SetsCS(ConstructList{6});
% maxNEmbryos=max(size(MeanAPCorr,1));
% MeanAPCorrAll=[];
% for ii=1:length(ConstructList)         %setting up array for 2way ANOVA
% MeanAPCorrAll=vertcat(MeanAPCorrAll,MeanAPCorr(:,:,ii));  %right now just look at 3 embryos (since have that for all constructs, until figure out unbalanced anova)
% end
% 
% 
% 
% figure
% C={'r','g','b','y','c','m'};
% for jj=1:length(ConstructList)  %want to plot construct means all on one graph for comparison 
%     h(jj)=plot(APBins,ConCorrMean(jj,:),'LineWidth',2,'Color',C{jj});%plot(APBins,h(jj),'LineWidth',2)
%     hold on
% end
% legend('Kr Both', 'Kr Dist', 'Kr Prox','Kr BothSep','Kr 2xDist','Kr 2xProx');
% title('Mean spot correlation all constructs')
% xlabel('Mean AP position')
% ylabel('Mean spot fluorescence correlation')
% ylim([-0.3 1]);
% xlim([0.3 0.75]);
%             
% figure 
% %NEmbryos=[];
% Data=[];
% 
% % for qq=1:length(ConstructList)
% %     Data(qq)= LoadMS2SetsCS(ConstructList{qq});
% % end
% BothEmbryos=length(BothData);
% for ii=1:BothEmbryos
%     q(ii)=plot(APBins,MeanAPCorr(ii,:,1),'o','LineWidth',2,'Color','r');
%     hold on
% end
% DistEmbryos=length(DistData);
% for qq=1:DistEmbryos
%     w(qq)=plot(APBins,MeanAPCorr(qq,:,2),'o','LineWidth',2,'Color','g');
% end
% ProxEmbryos=length(ProxData);
% for zz=1:ProxEmbryos
%     z(zz)=plot(APBins,MeanAPCorr(zz,:,3),'o','LineWidth',2,'Color','b');
% end
% BothSepEmbryos=length(BothSepData);
% for hh=1:BothSepEmbryos
%     h(hh)=plot(APBins,MeanAPCorr(hh,:,4), 'o', 'LineWidth', 2, 'Color', 'y');
% end
% DoubDistEmbryos=length(DoubDistData);
% for rr=1:DoubDistEmbryos
%     r(rr)=plot(APBins,MeanAPCorr(rr,:,5),'o','LineWidth',2,'Color','c');
% end
% DoubProxEmbryos=length(DoubProxData);
% for bb=1:DoubProxEmbryos
%     b(bb)=plot(APBins,MeanAPCorr(bb,:,5),'o','LineWidth',2,'Color','m');
% end
% legend([q(1), w(1), z(1), h(1),r(1),b(1)],{'Both','Dist','Prox','BothSep','2xDist','2xProx'});
%  xlabel('Mean AP position'); ylabel('Avg correlation of spot fluorescence');
%  title('Mean spot correlation all embryos');
%  ylim([0 1]);
% 
%  figure %plot the mean of correlation across the whole embryo for each embryo color-coded by construct 
%  for uu=1:BothEmbryos 
%      u(uu)=plot(1,MeanAllCorr(uu,1),'o','Color','r');
%      hold on 
%  end
%  for jj=1:DistEmbryos
%      t(jj)=plot(2,MeanAllCorr(jj,2),'o','Color','g');
%  end
%  for hh=1:ProxEmbryos
%      b(hh)=plot(3,MeanAllCorr(hh,3),'o','Color','b');
%  end
%  for h=1:BothSepEmbryos
%      as(h)=plot(4,MeanAllCorr(h,4),'o', 'Color', 'y');
%  end
%  for r=1:DoubDistEmbryos
%      ra(r)=plot(5,MeanAllCorr(r,5),'o','Color','c');
%  end
%  for v=1:DoubProxEmbryos
%      va(v)=plot(6,MeanAllCorr(v,6),'o','color','m');
%  end
%  legend([u(1),t(1),b(1),as(1),ra(1),va(1)],{'Both','Dist','Prox','BothSep','2xDist','2xProx'},'Location','best');
%  ylim([0 1]); xlim([0 7]);
%  ylabel('Mean spot correlation');
%  
%  title('Mean spot correlation across whole embryo');
%% calculate average mRNA produced per nucleus for each construct 
%load constructs
ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');

%First make sure have done SmoothingStuffnc14RW for all embryos
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        
        firstfilename=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
        load(firstfilename);
        if ncUse=='y'
            run('SmoothingStuffnc14RW.m')
            save([DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14'],'BurstProperties')
        %else
            %run SmoothingStuffRW.m
        %end
    end
    %clear BurstProperties
    end
end


%%
ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'}%{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
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
    ConmRNAProdAllAP=[];
    ConmRNAProdSE=[];
    ConProdSD=[];
    %ConDurSD=[];
     mRNAProdAllAP=[];
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
        ProductionSE(ee,aa,cc)=ProductionSD(ee,aa,cc)/sqrt(length(ProductionAP));             
            clear ProductionAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(mRNAProdAllAP,3)
            ConmRNAProdAllAP=[ConmRNAProdAllAP; mRNAProdAllAP(:,:,bb)];
        end
        for bb=1:size(ProductionSD,3)
            ConProdSD=[ConProdSD; ProductionSD(:,:,bb)];
        end
        for bb=1:size(ProductionSE,3)
            ConmRNAProdSE=[ConmRNAProdSE;ProductionSE(:,:,bb)];
        end
        EmbryoAvgProd(ee).MeanProd=nanmean(mRNAProdAllAP(:,:,ee));
        EmbryoAvgProd(ee).ProdSE=nanmean(mRNAProdErrorAllAP(:,:,ee));  %no idea really if this is right
%         AvgProdAllAP(cc).EmbryosProd(ee).MeanProd=nanmean(mRNAProdAllAP(:,:,ee));
%         AvgProdAllAP(cc).EmbryosProd(ee).SE=ProductionSE(ee,:,cc);
%         AvgProdAllAP(cc).EmbryosProd(ee).SD=ProductionSD(ee,:,cc);
%         
    end
        EmbsArray=[];
        for ee=1:length(EmbryoAvgProd)
            EmbsArray=[EmbsArray; [EmbryoAvgProd(ee).MeanProd]];  %3/16 10pm not working, not sure why, trying to put mean at each ap bin for each embryo in one array on top of one another 
        end
        AvgProdAllAP(cc).AvgProd=nanmean(ConmRNAProdAllAP,1);  %Avg mRNA production of all embryos of a construct by AP position
        AvgProdAllAP(cc).ConSD=nanmean(ConProdSD,1);
        AvgProdAllAP(cc).ConSE=nanmean(ConmRNAProdSE,1);
        AvgProdAllAP(cc).AllProds=[ConmRNAProdAllAP];
        for ee=1:length(EmbryoAvgProd)
        AvgProdAllAP(cc).EmbryoMeans=[EmbsArray];
        end 
        clear mRNAProdAllAP mRNAProdErrorAllAP;
end

%Plot mean of each embryo for each construct
% DistalColor=[8 180 238] ./ 255;
% DoubDistColor=[1 17 181] ./ 255;
% ProxColor=[251 220 50] ./ 255;
% DoubProxColor=[251 190 100] ./ 255;
% BothSepColor=[94 250 81] ./ 255;
% BothColor=[12 195 82] ./ 255;
% Colors=[BothColor; DistalColor; ProxColor; BothSepColor; DoubDistColor; DoubProxColor];
% APtoUse=input('Which AP bin to compare constructs?');
% EgglengthUse=APbinID(APtoUse)*100;
% for cc=1:length(ConstructList)
%     for ee=1:length(AvgProdAllAP(cc).EmbryosProd)
%         APBoxing(ee,cc)=AvgProdAllAP(cc).EmbryosProd(ee).MeanProd(APtoUse);
%     end
% end
% APBoxing(APBoxing==0)=nan;
% figure 
% 
%     for ee=1:length(AvgProdAllAP(1).EmbryosProd)
%     plot(1,AvgProdAllAP(1).EmbryosProd(ee).MeanProd(APtoUse),'o','LineWidth',1.5,'Color',DistalColor)
%     hold on
%     end
%     for ee=1:length(AvgProdAllAP(2).EmbryosProd)
%     plot(2,AvgProdAllAP(2).EmbryosProd(ee).MeanProd(APtoUse),'o','LineWidth',1.5,'Color',ProxColor)
%         
%     hold on
%     end
%     
%     for ee=1:length(AvgProdAllAP(3).EmbryosProd)
%     plot(3,AvgProdAllAP(3).EmbryosProd(ee).MeanProd(APtoUse),'o','LineWidth',1.5,'Color',BothSepColor)
%     hold on
%     end
%     for ee=1:length(AvgProdAllAP(4).EmbryosProd)
%     plot(4,AvgProdAllAP(4).EmbryosProd(ee).MeanProd(APtoUse),'o','LineWidth',1.5,'Color',DoubDistColor)
%     hold on
%     end
%     for ee=1:length(AvgProdAllAP(5).EmbryosProd)
%     plot(5,AvgProdAllAP(5).EmbryosProd(ee).MeanProd(APtoUse),'o','LineWidth',1.5,'Color',DoubProxColor)
%     hold on
%     end
%     for ee=1:length(AvgProdAllAP(6).EmbryosProd)
%     plot(6,AvgProdAllAP(6).EmbryosProd(ee).MeanProd(APtoUse),'o','LineWidth',1.5,'Color',BothColor);
%     %errorbar(6,AvgProdAllAP(6).EmbryosProd(ee).MeanProd(APtoUse),AvgProdAllAP(6).EmbryosProd(ee).SE(APtoUse),'o','LineWidth',1.5,'Color',BothColor)
%     hold on
%     end
%     boxplot(APBoxing,'Colors','byg')
% xlim([0 7]);
% xlabel('Construct');
% xticks([1:6]);
% xticklabels({'Dist', 'Prox', 'Both Sep', '2x Dist', '2x Prox', 'Both'});
% ylabel('Total nascent mRNA (AU)');
% title(['Mean mRNA production',' ' ,num2str(EgglengthUse),'% egg length'])


%ANOVA of each construct looking at variance with AP position
for cc=1:length(ConstructList)
    anova1(AvgProdAllAP(cc).AllProds);
    xlabel('AP bin')
    ylabel('Avg total nascent mRNA');
    title(ConstructList{cc});
end

%ANOVA at specific AP position of constructs against one another
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for bb=1:length([AvgProdAllAP(cc).AllProds])
    ConComp(bb,cc)=AvgProdAllAP(cc).AllProds(bb,APtoUse);
    end
    ConComp(ConComp==0)=nan;
end
[p,tbl,stats]=anova1(ConComp);
xlabel('Construct')
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
ylabel('Total nascent mRNA');
title(['Avg mRNA production AP bin',' ',num2str(APtoUse)]);

%plot the average across all embryos at each AP position for each construct
%as own graph
for cc=1:length(ConstructList)
    figure
    plot(1:length(APbinID),AvgProdAllAP(cc).AvgProd);
    hold on 
    errorbar(1:length(APbinID),AvgProdAllAP(cc).AvgProd,AvgProdAllAP(cc).ConSE,'o');
    xlabel('AP bin');
    ylabel('Total nascent mRNA');
    title(['Avg mRNA production', ConstructList{cc}]);
    ylim([0 1300000]);
    xlim([0 41]);
end

%plot constructs against one another as bar graphs at a single AP position
figure
for ii=1:length(AvgProdAllAP)
    AvgProdatAP(ii)=AvgProdAllAP(ii).AvgProd(APtoUse);

    errorbar(ii,AvgProdatAP(ii), AvgProdAllAP(ii).ConSE(APtoUse),'o');
    hold on
end
bar(AvgProdatAP,'EdgeColor','k','LineWidth',1.5)
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
xlim([0 7])
xlabel('Construct')
ylabel('Total nascent mRNA (AU)');
ylim([0 1200000]);
title(['Mean mRNA production',' ',num2str(EgglengthUse),'% egg length']);

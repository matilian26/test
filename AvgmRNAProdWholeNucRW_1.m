%load constructs
ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

%Ask if only want nc14 info
ncUse=input('Want to only use nc14?','s');
%Count for each construct
AvgmRNAwholenucProd=[];
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    AvgmRNAWholeNucProdCon=[];
    firsttime=1;
    ConmRNAWholeNucProdAllAP=[];
    ConmRNAWholeNucProdSE=[];
    ConWholeNucProdSD=[];
    %ConDurSD=[];
     mRNAWholeNucProdAllAP=[];
     mRNAWholeNucProdErrorAllAP=[];

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
            WholeNucProdAP=[];
            WholeNucProdErrorAP=[];
            WholeNucleusProductionAP=find([BurstProperties.APBin]==APbinID(aa));
            if isempty(WholeNucleusProductionAP)
                WholeNucProdAP=[WholeNucProdAP; nan];
            else
            for bb=1:length(WholeNucleusProductionAP)
                WholeNucProdAP=[WholeNucProdAP;BurstProperties(WholeNucleusProductionAP(bb)).BothAllelesProd];  %put all mRNA outputs at a given AP value in a column going down
                WholeNucProdErrorAP=[WholeNucProdErrorAP;BurstProperties(WholeNucleusProductionAP(bb)).BothAllelesError];           
            end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(WholeNucProdAP)
                mRNAWholeNucProdAllAP(bb,aa,ee)=WholeNucProdAP(bb);
            end
            for bb=1:length(WholeNucProdErrorAP)
                mRNAWholeNucProdErrorAllAP(bb,aa,ee)=WholeNucProdErrorAP(bb);
            end
            mRNAWholeNucProdAllAP(mRNAWholeNucProdAllAP==0)=nan;
            mRNAWholeNucProdErrorAllAP(mRNAWholeNucProdErrorAllAP==0)=nan;
            
            WholeNucleusProductionSD(ee,aa,cc)=nanstd(mRNAWholeNucProdAllAP(:,aa,ee));
        WholeNucleusProductionSE(ee,aa,cc)=WholeNucleusProductionSD(ee,aa,cc)/sqrt(length(WholeNucleusProductionAP));             
            clear WholeNucleusProductionAP
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        for bb=1:size(mRNAWholeNucProdAllAP,3)
            ConmRNAWholeNucProdAllAP=[ConmRNAWholeNucProdAllAP; mRNAWholeNucProdAllAP(:,:,bb)];
        end
        for bb=1:size(WholeNucleusProductionSD,3)
            ConWholeNucProdSD=[ConWholeNucProdSD; WholeNucleusProductionSD(:,:,bb)];
        end
        for bb=1:size(WholeNucleusProductionSE,3)
            ConmRNAWholeNucProdSE=[ConmRNAWholeNucProdSE;WholeNucleusProductionSE(:,:,bb)];
        end
        EmbryoAvgWholeNucProd(ee).MeanProd=nanmean(mRNAWholeNucProdAllAP(:,:,ee));
        EmbryoAvgWholeNucProd(ee).ProdSE=nanmean(mRNAWholeNucProdErrorAllAP(:,:,ee));  %no idea really if this is right
    end
        EmbsArray=[];
        for ee=1:length(EmbryoAvgWholeNucProd)
            EmbsArray=[EmbsArray; [EmbryoAvgWholeNucProd(ee).MeanProd]];  %3/16 10pm not working, not sure why, trying to put mean at each ap bin for each embryo in one array on top of one another 
        end
        AvgWholeNucProdAllAP(cc).AvgProd=nanmean(ConmRNAWholeNucProdAllAP,1);  %Avg mRNA production of all embryos of a construct by AP position
        AvgWholeNucProdAllAP(cc).ConSD=nanmean(ConWholeNucProdSD,1);
        AvgWholeNucProdAllAP(cc).ConSE=nanmean(ConmRNAWholeNucProdSE,1);
        AvgWholeNucProdAllAP(cc).AllProds=[ConmRNAWholeNucProdAllAP];
        for ee=1:length(EmbryoAvgWholeNucProd)
        AvgWholeNucProdAllAP(cc).EmbryoMeans=[EmbsArray];
        end 
        clear mRNAWholeNucProdAllAP mRNAWholeNucProdErrorAllAP;
end

%ANOVA of each construct looking at variance with AP position
for cc=1:length(ConstructList)
    anova1(AvgWholeNucProdAllAP(cc).AllProds);
    xlabel('AP bin')
    ylabel('Avg total nascent mRNA');
    title(ConstructList{cc});
end
%ANOVA at specific AP position of constructs against one another
APtoUse=input('Which AP bin to compare constructs?');
EgglengthUse=APbinID(APtoUse)*100;
for cc=1:length(ConstructList)
    for bb=1:length([AvgWholeNucProdAllAP(cc).AllProds])
    ConComp(bb,cc)=AvgWholeNucProdAllAP(cc).AllProds(bb,APtoUse);
    end
    ConComp(ConComp==0)=nan;
end
[p,tbl,stats]=anova1(ConComp);
xlabel('Construct')
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
ylabel('Total nascent mRNA per nucleus');
title(['Avg mRNA production AP bin',' ',num2str(APtoUse)]);

%plot the average across all embryos at each AP position for each construct
%as own graph
for cc=1:length(ConstructList)
    figure
    plot(1:length(APbinID),AvgWholeNucProdAllAP(cc).AvgProd);
    hold on 
    errorbar(1:length(APbinID),AvgWholeNucProdAllAP(cc).AvgProd,AvgWholeNucProdAllAP(cc).ConSE,'o');
    xlabel('AP bin');
    ylabel('Total nascent mRNA per nucleus');
    title(['Avg mRNA production', ConstructList{cc}]);
    ylim([0 1300000]);
    xlim([0 41]);
end

%plot constructs against one another as bar graphs at a single AP position
figure
for ii=1:length(AvgWholeNucProdAllAP)
    AvgWholeNucProdatAP(ii)=AvgWholeNucProdAllAP(ii).AvgProd(APtoUse);

    errorbar(ii,AvgWholeNucProdatAP(ii), AvgWholeNucProdAllAP(ii).ConSE(APtoUse),'o');
    hold on
end
bar(AvgWholeNucProdatAP,'EdgeColor','k','LineWidth',1.5)
xticks([1:6]);
xticklabels({'Both','Dist','Prox','BothSep','2xDist','2xProx'});
xlim([0 7])
xlabel('Construct')
ylabel('Total nascent mRNA (AU)');
ylim([0 1200000]);
title(['Mean mRNA production',' ',num2str(EgglengthUse),'% egg length']);
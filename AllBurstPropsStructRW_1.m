%% One structure for all burst properties 
ConstructList= {'KrDist','KrProx','KrBothSep','KrDistEmpty','KrProxEmpty','KrDistDuplicN','KrProxDuplic','KrBoth'} %{'KrDist','KrProx','KrBothSep', 'KrDistDuplicN', 'KrProxDuplic', 'KrBoth'};% %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
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
    EmbsArray=[];
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
                DurAP=[DurAP; nan];
                AmpAP=[AmpAP;nan];
                FreqAP=[FreqAP; nan]; 
            else
            for bb=1:length(ProductionAP)
                ProdAP=[ProdAP;BurstProperties(ProductionAP(bb)).TotalmRNA];  %put all mRNA outputs at a given AP value in a column going down
                ProdErrorAP=[ProdErrorAP;BurstProperties(ProductionAP(bb)).TotalmRNAError];  
                if ~isempty([BurstProperties(bb).Duration])
                DurAP=[DurAP; (BurstProperties(ProductionAP(bb)).Duration)'];
                else
                    DurAP=[DurAP; nan];
                end
                if ~isempty([BurstProperties(bb).BurstAmplitude]);
                    AmpAP=[AmpAP; ([BurstProperties(bb).BurstAmplitude])'];
                else
                    AmpAP=[AmpAP; nan];
                end
                FreqAP=[FreqAP; [BurstProperties(bb).Frequency]'];
            end
            end
            %make array with AP bins across as columns each page an embryo
            for bb=1:length(ProdAP)
                mRNAProdAllAP(bb,aa,ee)=ProdAP(bb);
            end
            for bb=1:length(ProdErrorAP)
                mRNAProdErrorAllAP(bb,aa,ee)=ProdErrorAP(bb);
            end
            for bb=1:length(DurAP)
                DurAllAP(bb,aa,ee)=DurAP(bb);
            end
            for bb=1:length(AmpAP)
                AmpAllAP(bb,aa,ee)=AmpAP(bb);
            end
            for bb=1:length(FreqAP)
                FreqAllAP(bb,aa,ee)=FreqAP(bb);
            end
            mRNAProdAllAP(mRNAProdAllAP==0)=nan;
            mRNAProdErrorAllAP(mRNAProdErrorAllAP==0)=nan;
            
            ProductionSD(ee,aa,cc)=nanstd(mRNAProdAllAP(:,aa,ee));
        ProductionSE(ee,aa,cc)=ProductionSD(ee,aa,cc)/sqrt(length(ProductionAP));
        DurationSD(ee,aa,cc)=nanstd(DurAllAP(:,aa,ee));
        DurationSE(ee,aa,cc)=DurationSD(ee,aa,cc)/(sqrt(length(ProductionAP));
        AmplitudeSD(ee,aa,cc)=nanstd(AmpAllAP(:,aa,ee));
        AmplitudeSE(ee,aa,cc)=AmplitudeSD(ee,aa,cc)/(sqrt(length(ProductionAP));
        FrequencySD(ee,aa,cc)=nanstd(FreqAllAP(:,aa,ee));
        FrequencySE(ee,aa,cc)=FrequencySD(ee,aa,cc)/(sqrt(length(ProductionAP)));
            clear ProductionAP
             
            
        end
        %Compile all raw duration data for a construct in one long column
        %to put in structure
        
         EmbryoAvg(ee).MeanProd=nanmean(mRNAProdAllAP(:,:,ee));
         EmbryoAvg(ee).ProdSE=nanmean(mRNAProdErrorAllAP(:,:,ee));  %no idea really if this is right
         for bb=1:size(DurAllAP,3)
             ConDurAllAP=[ConDurAllAP; DurAllAP(:,:,bb)];
             ConDurSD=[ConDurSD; DurationSD(:,:,bb)];
             ConDurSE=[ConDurSE; DurationSE(:,:,bb)];
             ConFreqAllAP=[ConFreqAllAP; FreqAllAP(:,:,bb)];
             ConAmpAllAP=[ConAmpAllAP; AmpAllAP(:,:,bb)];
         end
         AvgProdAllAP(cc).EmbryosProd(ee).MeanProd=nanmean(mRNAProdAllAP(:,:,ee));
        AvgProdAllAP(cc).EmbryosProd(ee).SE=ProductionSE(ee,:,cc);
         AvgProdAllAP(cc).EmbryosProd(ee).SD=ProductionSD(ee,:,cc);
         
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
        AvgProdAllAP(cc).AvgProd=nanmean(ConmRNAProdAllAP,1);  %Avg mRNA production of all embryos of a construct by AP position
        AvgProdAllAP(cc).ConSD=nanmean(ConProdSD,1);
        AvgProdAllAP(cc).ConSE=nanmean(ConmRNAProdSE,1);
        AvgProdAllAP(cc).AllProds=[ConmRNAProdAllAP];
        AvgProdAllAP(cc).AllSD=nanstd([AvgProdAllAP(cc).AllProds]);
        AvgProdAllAP(cc).AllSE=(([AvgProdAllAP(cc).AllSD])/sqrt(length(AvgProdAllAP(cc).AllProds)));
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
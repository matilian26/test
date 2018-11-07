%% calculate coefficient of variance for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDist32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
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
    AllIntraNoise=[];
    AllCoVarNoise=[];
    AllCorrSpots=[];
    AllTotalNoise=[];
    ConMeanTable=[];
    AllTotalmRNA=[];
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
        IntraNoise=nan(50, 41);
        CoVarNoise=nan(50,41);
        TotalNoiseVal=nan(50,41);
        CorrSpots=nan(50,41);
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
        MeanzTable=nan(MostAPNucs,41);
        TotalmRNA=nan(MostAPNucs,41);
        for aa=1:length(APbinID)
            APsubset=[];
%             VarPart=[];
%             SumVarPart=[];
            APsubset=SpotDiff(APstuff==APbinID(aa));
            DoneCounter=0;  % counting the number of nuclei actually used in comparisions
            if ~isempty(APsubset)
               
            for bb=1:length(APsubset)
                AvgSpotOne=[];
                AvgSpotTwo=[];
                DiffVal=[];
                SpotCorr=[];
                SquaredSum=[];
                MultVal=[];
                NumberNuclei(ee,aa,cc)=length(APsubset);
                TimeLength=[];
                TimeLength=sum(~isnan([APsubset(bb).SpotOne]));   %Only want to look at frames where nucleus exisits (i.e 0 or number values)
                AvgFluo=nanmean([APsubset(bb).SpotOne]);
                VarFluo=nanstd([APsubset(bb).SpotOne]);
                
            TimeTable(bb,aa)=(VarFluo/AvgFluo);
            MeanzTable(bb,aa)=AvgFluo;
            if ~isempty(APsubset(bb).TotalmRNAOne)
            TotalmRNA(bb,aa)=APsubset(bb).TotalmRNAOne;
            else
                TotalmRNA(bb,aa)=nan;
            end
            if isfield(APsubset,'SpotTwo') & (length(APsubset(bb).SpotOne)==length(APsubset(bb).SpotTwo))
                DoneCounter=DoneCounter+1;
            end
            for ss=1:length(APsubset(bb).SpotOne)
                if ~isfield(APsubset,'SpotTwo')
                    break
                else
                if length([APsubset(bb).SpotOne]) ~= length([APsubset(bb).SpotTwo])
                    break
                else
                    DiffVal(ss)=((APsubset(bb).SpotOne(ss) - APsubset(bb).SpotTwo(ss))^2);
                    SquaredSum(ss)=(((APsubset(bb).SpotOne(ss)^2))+(APsubset(bb).SpotTwo(ss))^2);
                    MultVal(ss)=(APsubset(bb).SpotOne(ss) * (APsubset(bb).SpotTwo(ss)));
                    SpotCorr(ss,1)=APsubset(bb).SpotOne(ss);
                    SpotCorr(ss,2)=APsubset(bb).SpotTwo(ss);
                end
                end
            end
            AvgDiffVal=nanmean(DiffVal);
            AvgSqrSum=nanmean(SquaredSum);
            AvgSpotOne=nanmean(APsubset(bb).SpotOne);
            VarSpotOne=var(APsubset(bb).SpotOne);
            AvgMultVal=nanmean(MultVal);
            IndSpotCorr=corrcoef(SpotCorr,'Rows','complete');
            if ~isfield(APsubset,'SpotTwo')
                break
%                 APsubset(bb).SpotTwo=APsubset(bb).SpotOne;
%                 APsubset(bb).SpotTwo(APsubset(bb).SpotOne>=1)=0;
            else
            %end (new addition)
                AvgSpotTwo=nanmean(APsubset(bb).SpotTwo);
                IntraNoise(bb,aa)=(AvgDiffVal/((2*(AvgSpotOne*AvgSpotTwo))));
                CoVarNoise(bb,aa)=(((AvgMultVal) - ((AvgSpotOne)*(AvgSpotTwo)))/((AvgSpotOne) * (AvgSpotTwo)));
                TotalNoiseVal(bb,aa)=((AvgSqrSum-(2*(AvgSpotOne)*(AvgSpotTwo)))/(2*(AvgSpotOne)*(AvgSpotTwo)));
            end
            if length(IndSpotCorr) >1
                CorrSpots(bb,aa)=IndSpotCorr(1,2);
            else
                CorrSpots(bb,aa)=nan;
            end
                %end
            
            %FluoValues(cc).Embryo(ee).Nucleus(bb).FluoValue=FluoVals;
            end
            EmbCounter(ee,aa,cc)=DoneCounter;
            EmbCounter(EmbCounter==0)=nan;
            for bb=1:length(APsubset)
                if ~isfield(APsubset,'SpotTwo')
                    break
                elseif ~isempty(APsubset(bb).SpotTwo)
                    TimeLength2=[];
                TimeLength2=(sum(~isnan(APsubset(bb).SpotTwo)));
                AvgFluo2=nanmean([APsubset(bb).SpotTwo]);
                VarFluo2=nanstd([APsubset(bb).SpotTwo]);
                
                TimeTable(bb+(length(APsubset)),aa)=(VarFluo2/AvgFluo2);
                if isfield(APsubset,'TotalmRNATwo') & (~isempty(APsubset(bb).TotalmRNATwo));
                TotalmRNA(bb+(length(APsubset)),aa)=APsubset(bb).TotalmRNATwo
                else
                    TotalmRNA(bb+(length(APsubset)),aa)=nan;
                end
                MeanzTable(bb+(length(APsubset)),aa)=AvgFluo2;
                end
            end
             TotalmRNA(TotalmRNA==0)=nan;
                IntraNoise(IntraNoise==0)=nan;
                IntraNoise(isinf(IntraNoise))=nan;
                CoVarNoise(CoVarNoise==0)=nan;
                CoVarNoise(isinf(CoVarNoise))=nan;
                TotalNoiseVal(TotalNoiseVal==0)=nan;
                TotalNoiseVal(isinf(TotalNoiseVal))=nan;
            end
        end
        %TimeTable(TimeTable==0)=nan;
        WholeNoise(cc).Embryo(ee).NoiseVals=IntraNoise;
        WholeNoise(cc).Embryo(ee).CoVarNoise=CoVarNoise;
        WholeNoise(cc).Embryo(ee).TotalNoiseCalc=TotalNoiseVal;
        WholeNoise(cc).Embryo(ee).SpotCorr=CorrSpots;
        WholeNoise(cc).Embryo(ee).TimeTable=TimeTable;
        WholeNoise(cc).Embryo(ee).TotalmRNA=TotalmRNA;
        WholeNoise(cc).Embryo(ee).MeanFluo=MeanzTable;
        WholeNoise(cc).Embryo(ee).EmbryoAvg=nanmean(TimeTable);
        AllNucTimeTable=[AllNucTimeTable;TimeTable];
        ConTimeAvg=[ConTimeAvg;(nanmean(TimeTable))];
        AllIntraNoise=[AllIntraNoise; IntraNoise];
        AllCoVarNoise=[AllCoVarNoise; CoVarNoise];
        AllTotalNoise=[AllTotalNoise;TotalNoiseVal];
        AllCorrSpots=[AllCorrSpots; CorrSpots];
        ConMeanTable=[ConMeanTable; MeanzTable];
        AllTotalmRNA=[AllTotalmRNA; TotalmRNA];
    end
    WholeNoise(cc).ConstructAvgNoise=nanmean(ConTimeAvg);
    WholeNoise(cc).AllNucsCV=AllNucTimeTable;
    WholeNoise(cc).AvgCVAllNucs=nanmean([WholeNoise(cc).AllNucsCV]);
    WholeNoise(cc).SDAllNucs=nanstd([WholeNoise(cc).AllNucsCV]);
    WholeNoise(cc).SEAllNucs=(WholeNoise(cc).SDAllNucs)/(sqrt(length(WholeNoise(cc).AllNucsCV)));
    WholeNoise(cc).AllIntraNoise=AllIntraNoise;
    WholeNoise(cc).SDAllIntraNoise=nanstd(AllIntraNoise);
    WholeNoise(cc).SEAllIntraNoise=(nanstd(AllIntraNoise)/(sqrt(length(AllIntraNoise(~isnan(AllIntraNoise))))));
    WholeNoise(cc).AllCoVarNoise=AllCoVarNoise;
    WholeNoise(cc).SDAllCoVarNoise=nanstd(AllCoVarNoise);
    WholeNoise(cc).SEAllCoVarNoise=(WholeNoise(cc).SDAllCoVarNoise)/(sqrt(length(AllCoVarNoise(~isnan(AllCoVarNoise)))));
    WholeNoise(cc).TotalNoise=AllTotalNoise;
    WholeNoise(cc).SDTotalNoise=nanstd(AllTotalNoise);
    WholeNoise(cc).SETotalNoise=(WholeNoise(cc).SDTotalNoise)/(sqrt(length(AllTotalNoise(~isnan(AllTotalNoise)))));
    WholeNoise(cc).AllCorrSpots=AllCorrSpots;
    WholeNoise(cc).SumSystNoise=[AllIntraNoise + AllCoVarNoise];
    WholeNoise(cc).AllAvgFluoVal=[ConMeanTable];
    WholeNoise(cc).AllTotalmRNA=[AllTotalmRNA];
    WholeNoise(cc).CountedNuclei=EmbCounter(:,:,cc);
end

% Get rid of places where only have one data point for a whole AP bin
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        if sum(~isnan(WholeNoise(cc).AllNucsCV(:,aa))) ==1
            WholeNoise(cc).AllNucsCV([find(~isnan(WholeNoise(cc).AllNucsCV(:,aa)))],aa)=nan;
            WholeNoise(cc).AvgCVAllNucs(aa)=nan;
        end
        NumberNuclei(cc,aa)=sum(~isnan(WholeNoise(cc).AllNucsCV(:,aa)));
    end
end

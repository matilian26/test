%% calculate coefficient of variance for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C';'KrDist17C';'Kr2xDist32C';'KrBoth17C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
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
    AllCorrTotmRNA=[];
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
        BothTotmRNA=nan(50,41);
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
            else
                if isempty(APsubset(bb).TotalmRNATwo) | (isempty(APsubset(bb).TotalmRNAOne))
                    BothTotmRNA(bb,aa)=nan;
                else
                    BothTotmRNA(bb,aa)=(APsubset(bb).TotalmRNAOne+APsubset(bb).TotalmRNATwo);
                end
                AvgSpotTwo=nanmean(APsubset(bb).SpotTwo);
                 IntraNoise(bb,aa)=(AvgDiffVal/((2*(AvgSpotOne*AvgSpotTwo))));
                 CoVarNoise(bb,aa)=(((AvgMultVal) - ((AvgSpotOne)*(AvgSpotTwo)))/((AvgSpotOne) * (AvgSpotTwo)));
                 TotalNoiseVal(bb,aa)=((AvgSqrSum-(2*(AvgSpotOne)*(AvgSpotTwo)))/(2*(AvgSpotOne)*(AvgSpotTwo)));
%                  IntraNoise(bb,aa)=sqrt(AvgDiffVal/((2*(AvgSpotOne*AvgSpotTwo))));
%                  CoVarNoise(bb,aa)=sqrt(((AvgMultVal) - ((AvgSpotOne)*(AvgSpotTwo)))/((AvgSpotOne) * (AvgSpotTwo)));
%                  TotalNoiseVal(bb,aa)=sqrt((AvgSqrSum-(2*(AvgSpotOne)*(AvgSpotTwo)))/(2*(AvgSpotOne)*(AvgSpotTwo))); 
              end
            if length(IndSpotCorr) >1
                CorrSpots(bb,aa)=IndSpotCorr(1,2);
            else
                CorrSpots(bb,aa)=nan;
            end
                %end
            
            %FluoValues(cc).Embryo(ee).Nucleus(bb).FluoValue=FluoVals;
            end
            for bb=1:length(APsubset)
                if ~isfield(APsubset,'SpotTwo')
                    break
                elseif ~isempty(APsubset(bb).SpotTwo)
                    TimeLength2=[];
                TimeLength2=(sum(~isnan(APsubset(bb).SpotTwo)));
                AvgFluo2=nanmean([APsubset(bb).SpotTwo]);
                VarFluo2=nanstd([APsubset(bb).SpotTwo]);
                
                TimeTable(bb+(length(APsubset)),aa)=(VarFluo2/AvgFluo2);
                if ~isempty(APsubset(bb).TotalmRNATwo);
                TotalmRNA(bb+(length(APsubset)),aa)=APsubset(bb).TotalmRNATwo
                else
                    TotalmRNA(bb+(length(APsubset)),aa)=nan;
                end
                MeanzTable(bb+(length(APsubset)),aa)=AvgFluo2;
                end
            end
             TotalmRNA(TotalmRNA==0)=nan;
                %IntraNoise(IntraNoise==0)=nan;
                IntraNoise(isinf(IntraNoise))=nan;
                %CoVarNoise(CoVarNoise==0)=nan;
                CoVarNoise(isinf(CoVarNoise))=nan;
                TotalNoiseVal(TotalNoiseVal==0)=nan;
                TotalNoiseVal(isinf(TotalNoiseVal))=nan;
            end
        end
        %TimeTable(TimeTable==0)=nan;
        CoVarNoise=real(CoVarNoise); %values were including tiny imaginary portions so get rid of those
        WholeNoise(cc).Embryo(ee).NoiseVals=IntraNoise;
        WholeNoise(cc).Embryo(ee).CoVarNoise=CoVarNoise;
        WholeNoise(cc).Embryo(ee).TotalNoiseCalc=TotalNoiseVal;
        WholeNoise(cc).Embryo(ee).SpotCorr=CorrSpots;
        WholeNoise(cc).Embryo(ee).TimeTable=TimeTable;
        WholeNoise(cc).Embryo(ee).TotalmRNA=TotalmRNA;
        WholeNoise(cc).Embryo(ee).MeanFluo=MeanzTable;
        WholeNoise(cc).Embryo(ee).EmbryoAvg=nanmean(TimeTable);
        
        AllCorrTotmRNA=[AllCorrTotmRNA; BothTotmRNA];
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
    WholeNoise(cc).CorrTotNoise=AllCorrTotmRNA;
    WholeNoise(cc).AllNucsCV=AllNucTimeTable;
    WholeNoise(cc).AvgCVAllNucs=nanmean([WholeNoise(cc).AllNucsCV]);
    WholeNoise(cc).SDAllNucs=nanstd([WholeNoise(cc).AllNucsCV]);
    WholeNoise(cc).SEAllNucs=(WholeNoise(cc).SDAllNucs)/(sqrt(length(WholeNoise(cc).AllNucsCV)));
    WholeNoise(cc).AllIntraNoise=AllIntraNoise;
    WholeNoise(cc).SDAllIntraNoise=nanstd(AllIntraNoise);
    WholeNoise(cc).SEAllIntraNoise=(nanstd(AllIntraNoise)/(sqrt(length(AllIntraNoise(~isnan(AllIntraNoise))))));
    WholeNoise(cc).Conf95InterNoise=(WholeNoise(cc).SEAllIntraNoise).*1.96;
    WholeNoise(cc).AllCoVarNoise=AllCoVarNoise;
    WholeNoise(cc).SDAllCoVarNoise=nanstd(AllCoVarNoise);
    WholeNoise(cc).SEAllCoVarNoise=(WholeNoise(cc).SDAllCoVarNoise)/(sqrt(length(AllCoVarNoise(~isnan(AllCoVarNoise)))));
    WholeNoise(cc).Conf95CoVar=(WholeNoise(cc).SEAllCoVarNoise).*1.96;
    WholeNoise(cc).TotalNoise=AllTotalNoise;
    WholeNoise(cc).SDTotalNoise=nanstd(AllTotalNoise);
    WholeNoise(cc).SETotalNoise=(WholeNoise(cc).SDTotalNoise)/(sqrt(length(AllTotalNoise(~isnan(AllTotalNoise)))));
    WholeNoise(cc).Conf95TotalNoise=(WholeNoise(cc).SETotalNoise).*1.96;
    WholeNoise(cc).AllCorrSpots=AllCorrSpots;
    WholeNoise(cc).SumSystNoise=[AllIntraNoise + AllCoVarNoise];
    WholeNoise(cc).AllAvgFluoVal=[ConMeanTable];
    WholeNoise(cc).AllTotalmRNA=[AllTotalmRNA];
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

%% Check if allele noise + co-variance equal total noise
for cc=1:length(ConstructList)
    WholeNoise(cc).IssuePoints=nan(length(WholeNoise(cc).TotalNoise),41);
    for aa=1:length(APbinID)
        for bb=1:length(WholeNoise(cc).TotalNoise)
            AlleleNoise=WholeNoise(cc).AllIntraNoise(bb,aa);
            CoVarVal=WholeNoise(cc).AllCoVarNoise(bb,aa);
            TotalNoiseVal=WholeNoise(cc).TotalNoise(bb,aa);
            if ~isnan(AlleleNoise) & (~isnan(CoVarVal))
            if (AlleleNoise + CoVarVal) ~= TotalNoiseVal
                WholeNoise(cc).IssuePoints(bb,aa)=(TotalNoiseVal-(AlleleNoise+CoVarVal));
            elseif (AlleleNoise + CoVarVal) == TotalNoiseVal
                WholeNoise(cc).IssuePoints(bb,aa)=1;
            end
            end
        end
    end
end

%% Plotting
DistalColor=[1 64 172]./255;
Distal32CColor=[118 180 238] ./ 255;
DistalEmptyColor=[8 210 238] ./ 255; 
DoubDistColor=[73 184 253] ./ 255;
ProxColor=[238 123 23]./255;
ProxEmptyColor=[251 250 50] ./255;
Proximal32CColor=[251 150 10] ./ 255;
DoubProxColor=[215 183 58] ./ 255;
DoubProxEmptyColor=[251 220 50] ./ 255;
BothSepColor=[94 250 81] ./ 255;
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
Colors(7).Color=Distal32CColor;
Colors(8).Color=Proximal32CColor;
Colors(9).Color=BothSep32CColor;
Colors(10).Color=Both32CColor;
Colors(11).Color=DoubProx32CColor;
Colors(12).Color=DistalColor;
Colors(13).Color=DoubDistColor;
Colors(14).Color=Both32CColor;

fontsize=18;
fontname='Helvetica';

EggLength=APbinID .* 100;

FigDirect=[DropboxFolder filesep 'Figures'];
%%
% Stacked bar graph of noise contributions at AP bin of max mRNA production
% of that construct
for cc=1:length(ConstructList)
    APbinToUse=find(nanmean(WholeNoise(cc).AllTotalmRNA)==max(nanmean(WholeNoise(cc).AllTotalmRNA)));
    NoiseContributions(cc,1)=nanmean(WholeNoise(cc).AllIntraNoise(:,APbinToUse));
    NoiseContributions(cc,2)=nanmean(WholeNoise(cc).AllCoVarNoise(:,APbinToUse));
    NoiseContributions(cc,3)=nanmean(WholeNoise(cc).TotalNoise(:,APbinToUse));
    NoiseContributions(cc,4)=APbinToUse;
    NoiseContributions(cc,5)=WholeNoise(cc).Conf95InterNoise(APbinToUse);
    NoiseContributions(cc,6)=WholeNoise(cc).Conf95CoVar(APbinToUse);
    NoiseContributions(cc,7)=nanmedian(WholeNoise(cc).TotalNoise(:,APbinToUse));
    NoiseContributions(cc,8)=((NoiseContributions(cc,3)-NoiseContributions(cc,7))/(NoiseContributions(cc,3)))*100; % % mean is larger than median
    %Median values for intra,co-variance,total noise
    NoiseContributions(cc,9)=nanmedian(WholeNoise(cc).AllIntraNoise(:,APbinToUse));
    NoiseContributions(cc,10)=nanmedian(WholeNoise(cc).AllCoVarNoise(:,APbinToUse));
    NoiseContributions(cc,11)=nanmedian(WholeNoise(cc).TotalNoise(:,APbinToUse));
end
APBoxing=[];
for cc=1:length(ConstructList)
    if cc>1
    if length(WholeNoise(cc).TotalNoise)>size(APBoxing,1)
        
        APOrig=size(APBoxing,1);
        APBoxing([APOrig:length(WholeNoise(cc).TotalNoise)],:)=nan;
    elseif length(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)))<size(APBoxing,1)
        OrigLength=length(WholeNoise(cc).TotalNoise);
        WholeNoise(cc).TotalNoise([OrigLength:size(APBoxing,1)],:)=nan;
        %added 8/2 to allow scatter plots
        WholeNoise(cc).CorrTotNoise([OrigLength:size(APBoxing,1)],:)=nan;
        
    end
    end
    APBoxing=[APBoxing, [WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))]];
end
%Make "box & whisker" type plots
figure
PlaceCounter=0;
for cc=[1,2,4,5,6]%length(ConstructList)
    PlaceCounter=PlaceCounter+1;
    PointsUsed=length(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    %plot((PlaceCounter.*ones(PointsUsed,1)),WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),'jitter
    %scatter((PlaceCounter.*ones(PointsUsed,1)),(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),[],Colors(cc).Color,'jitter','on','jitterAmount',0.3);
    scatter((PlaceCounter.*ones(PointsUsed,1)),((WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)))),[],Colors(cc).Color,'jitter','on','jitterAmount',0.3);
    SevPerc(cc)=prctile((WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),75);
    hold on 
    %errorbar(PlaceCounter, nanmedian((WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)))),(WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4))),'.','Color','k','LineWidth',2);
    %errorbar(PlaceCounter, nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4)),'.','Color','k','LineWidth',2);
    %y=nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    y=nanmedian((WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))));
    h=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); h.Color='k'; h.LineWidth=2;
end
MaxPerc=max(SevPerc);
%DataBox=boxplot(APBoxing(:,[1,2,5,6]),'Colors','k');
%ylim([0 round(MaxPerc)]);
ylim([0 MaxPerc]);
xlim([0 6]);
xticks(1:5);
xticklabels({'Dist','Prox','2x Dist','2xProx','Both'})%'Dist32','Sep32','Both32','2xProx32'});
%title('All nuclei');
ylabel('total noise');
set(gca, 'FontName', fontname, 'FontSize', fontsize);
print( [FigDirect filesep 'AllNucleiBoxplot'],'-dsvg');

ColorUsed={Colors(1).Color,Colors(2).Color,Colors(3).Color,Colors(4).Color,Colors(5).Color,Colors(6).Color,Colors(7).Color, Colors(8).Color,Colors(9).Color,Colors(10).Color};
% %Violin plot
% figure 
% distributionPlot(APBoxing,'color',ColorUsed);
% ylim([0 round(MaxPerc)]);
% xticklabels({'Dist','Prox','Sep','2xDist','2xProx','Both','Dist32','Sep32','Both32','2xProx32'});
% ylabel('Total noise');
% title('All nuclei');

% points of embryo's avgs 
ConBox=[];
for cc=1:length(ConstructList)
    EmbBox=[];
    for ee=1:length(WholeNoise(cc).Embryo)
        EmbBox=[EmbBox; nanmean(WholeNoise(cc).Embryo(ee).TotalNoiseCalc(:,NoiseContributions(cc,4)))];
    end
    if cc >1
        if ee > size(ConBox,1)
            ConBox([size(ConBox,1):ee],:)=nan;
        else if ee < size(ConBox,1)
                EmbBox([end+1:size(ConBox,1)],1)=nan;
            end
        end
    end
    ConBox=[ConBox,EmbBox];
end
figure
for cc=1:length(ConstructList)
    EmbPerc(cc)=prctile(ConBox(:,cc),75);
    for ee=1:length(WholeNoise(cc).Embryo)
        plot(cc,nanmean(WholeNoise(cc).Embryo(ee).TotalNoiseCalc(:,NoiseContributions(cc,4))),'o','Color',Colors(cc).Color,'LineWidth',2);
        hold on 
    end
end
MaxPrct=max(EmbPerc);
boxplot(ConBox,'Colors','k');
title('Embryo averages');
ylabel('Total noise');
xticklabels({'Dist','Prox','Sep','2xDist','2xProx','Both','Dist32','Prox32C','Sep32','Both32','2xProx32'});
ylim([0 (MaxPrct)]);

%% bar graphs of inter-allele and co-variance next to each other
figure
NoiseStuff=[];
NoiseError=[];
Counter=0;
for cc=[1,2,5,6]
    Counter=Counter+1;
    NoiseError=[NoiseError; [NoiseContributions(cc,6),NoiseContributions(cc,5)]];
    %NoiseStuff=[NoiseStuff;[NoiseContributions(cc,2),NoiseContributions(cc,1)]];
    %median values
    NoiseStuff=[NoiseStuff; [NoiseContributions(cc,10),NoiseContributions(cc,9)]];
end
b=bar(NoiseStuff,'BarWidth',1,'FaceColor','flat');
PlaceCounter=0;
for cc=[1,2,5,6]
    PlaceCounter=PlaceCounter+1;
    b(1).CData(PlaceCounter,:)=[Colors(cc).Color];
    b(2).CData(PlaceCounter,:)=[Colors(cc).Color];
    %y=nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    y=nanmedian(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    h=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); h.Color='k'; h.LineWidth=2;
end
ngroups=size(NoiseStuff,1);
groupwidth=min(0.8, 2/3.5);
hold on 
 for bb=1:2
     x=(1:ngroups)-groupwidth/2 + (2*bb-1) * groupwidth / 4;
     errorbar(x, NoiseStuff(:,bb), NoiseError(:,bb),'Color','k','LineWidth',2,'LineStyle','none');
 end
xticks([1:4]);
xticklabels({'Distal','Proximal','2x Proximal','Both'});
ylabel('Noise');
xlabel('Construct');
set(gca,'FontSize',fontsize,'FontName',fontname);
print( [FigDirect filesep 'NoiseBars'],'-dsvg');

%% Scatter plot style boxes of inter-allele & co-variance 
figure
PlaceCounter=1;
for cc=[1,2,4,5,6]%length(ConstructList)
    
    PointsUsed=length(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)));
    %plot((PlaceCounter.*ones(PointsUsed,1)),WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),'jitter
    scatter(((PlaceCounter-0.5).*ones(PointsUsed,1)),WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4)),20,Colors(cc).Color,'jitter','on','jitterAmount',0.3);
     hold on 
    scatter(((PlaceCounter+0.3).*ones(PointsUsed,1)),WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)),20,Colors(cc).Color,'x','jitter','on','jitterAmount',0.3);
    SevPerc(cc)=prctile(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),75);
    %errorbar((PlaceCounter-0.5), nanmedian(WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95CoVar(NoiseContributions(cc,4)),'.','Color','k','LineWidth',2);
    %errorbar((PlaceCounter+0.3), nanmedian(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95InterNoise(NoiseContributions(cc,4)),'.','Color','k','LineWidth',2);
    
    SevFivPercInter(cc)=prctile(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)),75);
    SevFivPercCoVar(cc)=prctile(WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4)),75);
    %errorbar(PlaceCounter, nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4)),'.','Color','k','LineWidth',2);
    %y=nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    y=nanmedian(WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4)));
    y2=nanmedian(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)));
    h=line([PlaceCounter-0.9,PlaceCounter-0.1],[y,y]); h.Color='k'; h.LineWidth=2;
    h2=line([PlaceCounter-0.1,PlaceCounter+0.7],[y2,y2]); h2.Color='k'; h2.LineWidth=2;
    PlaceCounter=PlaceCounter+2.4;
    %PlaceCounter=PlaceCounter+0.6;
end
SevFivPercInter=max(SevFivPercInter); SevFivPercCoVar=max(SevFivPercCoVar);
if SevFivPercInter > SevFivPercCoVar
    SevFivPerc=(0.75*SevFivPercInter);
    %SevFivPerc=(0.1*SevFivPercInter);
else
    SevFivPerc=(0.75*SevFivPercCoVar);
    %SevFivPerc=(0.1*SevFivPercCoVar);
end
set(gca,'FontSize',fontsize,'FontName',fontname);
ylabel('noise');
ylim([0 SevFivPerc])
xlim([0 12])
xticks([1,2.5,4,5.5,7]);
xticklabels({'Dist','Prox','2x Dist','2x Prox','Both'});
print( [FigDirect filesep 'NoiseBarsNuclei'],'-dsvg');

%% Temp compaisons scatter/bar graphs
figure
PlaceCounter=0;
for cc=[1,7,2,8,4,13,5,11,6,10]%length(ConstructList)
    PlaceCounter=PlaceCounter+1;
    PointsUsed=length(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    %plot((PlaceCounter.*ones(PointsUsed,1)),WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),'jitter
    %scatter((PlaceCounter.*ones(PointsUsed,1)),(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),[],Colors(cc).Color,'jitter','on','jitterAmount',0.3);
    scatter((PlaceCounter.*ones(PointsUsed,1)),((WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)))),[],Colors(cc).Color,'jitter','on','jitterAmount',0.3);
    SevPerc(cc)=prctile((WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),75);
    hold on 
    %errorbar(PlaceCounter, nanmedian((WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)))),(WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4))),'.','Color','k','LineWidth',2);
    %errorbar(PlaceCounter, nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4)),'.','Color','k','LineWidth',2);
    %y=nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    y=nanmedian((WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))));
    h=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); h.Color='k'; h.LineWidth=2;
end
MaxPerc=max(SevPerc);
%DataBox=boxplot(APBoxing(:,[1,2,5,6]),'Colors','k');
%ylim([0 round(MaxPerc)]);
ylim([0 MaxPerc]);
xlim([0 11]);
xticks(1:10);
xticklabels({'Dist','Prox','2x Dist','2xProx','Both'})%'Dist32','Sep32','Both32','2xProx32'});
%title('All nuclei');
ylabel('total noise');
set(gca, 'FontName', fontname, 'FontSize', fontsize);

figure
PlaceCounter=1;
for cc=[1,7,2,8,4,13,5,11,6,10]%length(ConstructList)
    
    PointsUsed=length(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)));
    %plot((PlaceCounter.*ones(PointsUsed,1)),WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),'jitter
    scatter(((PlaceCounter-0.5).*ones(PointsUsed,1)),WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4)),20,Colors(cc).Color,'jitter','on','jitterAmount',0.3);
     hold on 
    scatter(((PlaceCounter+0.3).*ones(PointsUsed,1)),WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)),20,Colors(cc).Color,'x','jitter','on','jitterAmount',0.3);
    SevPerc(cc)=prctile(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),75);
    %errorbar((PlaceCounter-0.5), nanmedian(WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95CoVar(NoiseContributions(cc,4)),'.','Color','k','LineWidth',2);
    %errorbar((PlaceCounter+0.3), nanmedian(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95InterNoise(NoiseContributions(cc,4)),'.','Color','k','LineWidth',2);
    
    SevFivPercInter(cc)=prctile(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)),75);
    SevFivPercCoVar(cc)=prctile(WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4)),75);
    %errorbar(PlaceCounter, nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4)),'.','Color','k','LineWidth',2);
    %y=nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    y=nanmedian(WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4)));
    y2=nanmedian(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)));
    h=line([PlaceCounter-0.7,PlaceCounter-0.3],[y,y]); h.Color='k'; h.LineWidth=2;
    h2=line([PlaceCounter+0.1,PlaceCounter+0.5],[y2,y2]); h2.Color='k'; h2.LineWidth=2;
    PlaceCounter=PlaceCounter+2.4;
    %PlaceCounter=PlaceCounter+0.6;
end
SevFivPercInter=max(SevFivPercInter); SevFivPercCoVar=max(SevFivPercCoVar);
if SevFivPercInter > SevFivPercCoVar
    SevFivPerc=(0.75*SevFivPercInter);
    %SevFivPerc=(0.1*SevFivPercInter);
else
    SevFivPerc=(0.75*SevFivPercCoVar);
    %SevFivPerc=(0.1*SevFivPercCoVar);
end
set(gca,'FontSize',fontsize,'FontName',fontname);
ylabel('noise');
ylim([0 SevFivPerc])
xlim([0 24])
xticks([1,2.5,4,5.5,7]);
xticklabels({'Dist','Prox','2x Dist','2x Prox','Both'});
%%

figure
bar(NoiseContributions([1:6],[1:2]),'stacked')
xticklabels({'Dist','Prox','Separated','2xDist','2xProx','Both'})%,'Dist32C','BothSep32C','Both32C'})
title('Total noise max mRNA production AP bin');
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
hold on
%errorbar(cumsum(NoiseContributions(:,[1:2])')',NoiseContributions(:,[5:6]), 'Color','k','LineStyle','none');
 for cc=1:6%length(ConstructList)
     APbinCon=NoiseContributions(cc,4);
     errorbar(cc,nanmean(WholeNoise(cc).AllIntraNoise(:,APbinCon)),WholeNoise(cc).SEAllIntraNoise(APbinCon),'Color','k','LineWidth',1.5);
     errorbar(cc,nanmean(WholeNoise(cc).AllCoVarNoise(:,APbinCon)),WholeNoise(cc).SEAllCoVarNoise(APbinCon),'Color','k','LineWidth',1.5);
 end

 figure
bar(NoiseContributions([1,3,5,6,7,8,9,11,10],[1:2]),'stacked')
xticklabels({'Dist','Sep','2xProx','Both','Dist32C','Prox32C','Sep32C','2xProx32C','Both32C'})
title('Total noise max mRNA production AP bin');
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
hold on
%errorbar(cumsum(NoiseContributions(:,[1:2])')',NoiseContributions(:,[5:6]), 'Color','k','LineStyle','none');
bb=0; 
for cc=[1,3,5,6,7,8,9,11,10]%length(ConstructList)
    bb=bb+1; 
     APbinCon=NoiseContributions(cc,4);
     errorbar(bb,nanmean(WholeNoise(cc).AllIntraNoise(:,APbinCon)),WholeNoise(cc).SEAllIntraNoise(APbinCon),'Color','k','LineWidth',1.5);
     errorbar(bb,nanmean(WholeNoise(cc).AllCoVarNoise(:,APbinCon)),WholeNoise(cc).SEAllCoVarNoise(APbinCon),'Color','k','LineWidth',1.5);
     
 end
 
for cc=1:length(ConstructList)
    APbinToUse=find(nanmean(WholeNoise(cc).AllTotalmRNA)==max(nanmean(WholeNoise(cc).AllTotalmRNA)));
    PercentNoiseContributions(cc,1)=(nanmean(WholeNoise(cc).AllIntraNoise(:,APbinToUse))/(nanmean(WholeNoise(cc).TotalNoise(:,APbinToUse)))).*100;
    PercentNoiseContributions(cc,2)=(nanmean(WholeNoise(cc).AllCoVarNoise(:,APbinToUse))/(nanmean(WholeNoise(cc).TotalNoise(:,APbinToUse)))).*100;
    PercentNoiseContributions(cc,3)=nanmean(WholeNoise(cc).TotalNoise(:,APbinToUse));
    PercentNoiseContributions(cc,4)=APbinToUse;
    NucsUsing=WholeNoise(cc).AllIntraNoise(:,APbinCon);
    NucsUsingCo=WholeNoise(cc).AllCoVarNoise(:,APbinCon);
    PercentNoiseContributions(cc,5)=sqrt(((PercentNoiseContributions(cc,1))/100)*(1-(PercentNoiseContributions(cc,1)/100))/(length(NucsUsing(~isnan(NucsUsing)))));
    PercentNoiseContributions(cc,6)=sqrt(((PercentNoiseContributions(cc,2))/100)*(1-(PercentNoiseContributions(cc,2)/100))/(length(NucsUsingCo(~isnan(NucsUsingCo)))));

end
figure
bar(PercentNoiseContributions([1:6],[1:2]),'stacked');
ylabel('% of total noise')
yticks([0:10:100]);
ylim([0 105]);
xticklabels({'Dist','Prox','BothSep','2xDist','2xProx','Both'})%'Dist32C','BothSep32C','Both32C','2xProx32C'})
title('Noise contribution at max mRNA production AP bin');
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
hold on 
%errorbar(cumsum(PercentNoiseContributions(:,[1:2])')',(PercentNoiseContributions(:,[5:6])*100),'Color','k','LineStyle','none');
for cc=1:6%length(ConstructList)
    errorbar(cc,PercentNoiseContributions(cc,1),(PercentNoiseContributions(cc,5))*100,'Color','k','LineWidth',1.5);
    errorbar(cc,PercentNoiseContributions(cc,2),(PercentNoiseContributions(cc,6))*100,'Color','k','LineWidth',1.5);
end
%%

%Compare Temps
figure
PlaceHolder=0;
for cc=[1,7]
    PlaceHolder=PlaceHolder+1;
    PointsUsed=length(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    scatter((PlaceHolder.*ones(PointsUsed,1)),WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),'jitter','on','jitteramount',0.3);
    hold on 
    errorbar(PlaceHolder, nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4)),'Color','k','LineWidth',2.5);
    h=line([PlaceHolder-0.5,PlaceHolder+0.5],[nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)))]);
    h.Color='k'; h.LineWidth=2.5;
end
set(gca,'FontSize',fontsize,'FontName',fontname);
ylim([0 6]);
xticks([1:2]);
xticklabels({'distal','distal 32C'});
ylabel('total noise');

figure
NoiseStuff=[];
NoiseError=[];
Counter=0;
for cc=[1,7]
    Counter=Counter+1;
    NoiseError=[NoiseError; [NoiseContributions(cc,6),NoiseContributions(cc,5)]];
    NoiseStuff=[NoiseStuff;[NoiseContributions(cc,2),NoiseContributions(cc,1)]];
end
b=bar(NoiseStuff,'BarWidth',1,'FaceColor','flat');
PlaceCounter=0;
for cc=[1,7]
    PlaceCounter=PlaceCounter+1;
    b(1).CData(PlaceCounter,:)=[Colors(cc).Color];
    b(2).CData(PlaceCounter,:)=[Colors(cc).Color];
    y=nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    h=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); h.Color='k'; h.LineWidth=2;
end
ngroups=size(NoiseStuff,1);
groupwidth=min(0.8, 2/3.5);
hold on 
 for bb=1:2
     x=(1:ngroups)-groupwidth/2 + (2*bb-1) * groupwidth / 4;
     errorbar(x, NoiseStuff(:,bb), NoiseError(:,bb),'Color','k','LineWidth',2,'LineStyle','none');
 end
 xticklabels({'distal', 'distal 32C'});
 set(gca,'FontSize',fontsize,'FontName',fontname);
ylabel('noise');

figure
bar(NoiseContributions([1,7],[1:2]),'stacked')
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
xticklabels({'Distal','Distal 32C'});
hold on
%errorbar(cumsum(NoiseContributions([1,7],[1:2])')',NoiseContributions([1,7],[5:6]), 'Color','k','LineStyle','none');
errorbar(1,NoiseContributions(1,1),NoiseContributions(1,5),'Color','k','LineWidth',1.5);
errorbar(1, NoiseContributions(1,2),NoiseContributions(1,6),'Color','k','LineWidth',1.5);
errorbar(2,NoiseContributions(7,1),NoiseContributions(7,5),'Color','k','LineWidth',1.5);
errorbar(2, NoiseContributions(7,2),NoiseContributions(7,6),'Color','k','LineWidth',1.5);

figure 
scatter(WholeNoise(1).CorrTotNoise(:), WholeNoise(1).TotalNoise(:), [], 'g');
hold on 
scatter(WholeNoise(7).CorrTotNoise(:), WholeNoise(7).TotalNoise(:), [], 'r');
ylabel('Total Noise');
xlabel('Total mRNA production');
legend('Distal', 'Distal 32C');
%%
figure 
[sorted_x, sorted_x_index] = sort(WholeNoise(1).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(1).TotalNoise(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'color',Colors(1).Color,'LineWidth',2.5)
RTNoisevmRNA=[sorted_x, smooth_y];

[sorted_x, sorted_x_index] = sort(WholeNoise(7).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(7).TotalNoise(sorted_x_index),70);
plot(sorted_x, smooth_y, 'color',Colors(1).Color,'LineWidth',1.5,'LineStyle','--')
ylabel('total noise');
xlabel('total mRNA production');
HotNoisevmRNA=[sorted_x,smooth_y];
%legend('Distal RT', 'Distal 32C');
set(gca, 'FontSize',fontsize, 'FontName', fontname);
print( [FigDirect filesep 'NoisevExpDistal'],'-dsvg');
[H, pval,KStat]=kstest_2s_2d(RTNoisevmRNA,HotNoisevmRNA);

%% Co-variance vs total mRNA expression 32C vs RT
figure 
[sorted_x, sorted_x_index] = sort(WholeNoise(1).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(1).AllCoVarNoise(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'color',Colors(1).Color,'LineWidth',2.5)
RTNoisevmRNA=[sorted_x, smooth_y];

[sorted_x, sorted_x_index] = sort(WholeNoise(7).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(7).AllCoVarNoise(sorted_x_index),70);
plot(sorted_x, smooth_y, 'color',Colors(1).Color,'LineWidth',1.5,'LineStyle','--')
ylabel('total noise');
xlabel('total mRNA production');
HotNoisevmRNA=[sorted_x,smooth_y];
%legend('Distal RT', 'Distal 32C');
set(gca, 'FontSize',fontsize, 'FontName', fontname);
%% 17C compare
figure 
RTNoisevmRNA=[]; HotNoisevmRNA=[];
[sorted_x, sorted_x_index] = sort(WholeNoise(1).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(1).TotalNoise(sorted_x_index),70);
RTNoisevmRNA=[sorted_x, smooth_y];
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(1).Color,'LineWidth',2.5)

[sorted_x, sorted_x_index] = sort(WholeNoise(12).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(12).TotalNoise(sorted_x_index),70);
HotNoisevmRNA=[sorted_x, smooth_y];
plot(sorted_x, smooth_y, 'Color',Colors(1).Color,'LineWidth',1.5,'LineStyle','--')
ylabel('total noise');
xlabel('total mRNA production');
%legend('Proximal RT', 'Proximal 32C');
set(gca, 'FontSize',fontsize, 'FontName', fontname);
%print( [FigDirect filesep 'NoisevExpProximal'],'-dsvg');
[H, pval,KStat]=kstest_2s_2d(RTNoisevmRNA,HotNoisevmRNA);

%%
figure
[ind1, ind2]=find(WholeNoise(1).AllCoVarNoise<0);
for aa=1:length(APbinID)
    NegFract(aa)=length(find(ind2==aa));
    NegFract(aa)=NegFract(aa)/(sum(~isnan(WholeNoise(1).AllCoVarNoise(:,aa))));
end
plot(EggLength, (NegFract.*100),'Color',Colors(1).Color, 'LineWidth',2.5);
ylabel('% of nuclei with negative co-variance');
xlabel('% egg length');
xlim([0 100])
hold on 
h=line([EggLength(NoiseContributions(1,4)) EggLength(NoiseContributions(1,4))],[0 20]); 
h.Color='g'; h.LineWidth=2;
[ind1, ind2]=find(WholeNoise(7).AllCoVarNoise<0);
for aa=1:length(APbinID)
    NegFract(aa)=length(find(ind2==aa));
    NegFract(aa)=NegFract(aa)/(sum(~isnan(WholeNoise(7).AllCoVarNoise(:,aa))));
end
plot(EggLength, NegFract.*100, 'Color', Colors(7).Color,'LineWidth',2.5,'LineStyle','-.');
h2=line([EggLength(NoiseContributions(7,4)) EggLength(NoiseContributions(7,4))],[0 20]);
h2.LineWidth=2; h2.Color='r';
%%
PercDistIntra=(nanmean(WholeNoise(7).AllIntraNoise))-(nanmean(WholeNoise(1).AllIntraNoise));
PercDistIntra=(PercDistIntra./(nanmean(WholeNoise(1).AllIntraNoise))).*100;
PercDistCoVar=(nanmean(WholeNoise(7).AllCoVarNoise))-(nanmean(WholeNoise(1).AllCoVarNoise));
PercDistCoVar=(PercDistCoVar./(nanmean(WholeNoise(1).AllCoVarNoise))).*100;
PercDistTot=(nanmean(WholeNoise(7).TotalNoise))-(nanmean(WholeNoise(1).TotalNoise));
PercDistTot=(PercDistTot./(nanmean(WholeNoise(1).TotalNoise))).*100;
figure
plot(EggLength,PercDistIntra,'LineWidth',2);
hold on 
plot(EggLength,PercDistCoVar,'LineWidth',2);
plot(EggLength,PercDistTot,'LineWidth',2);
CutOff=refline(0,0); CutOff.Color='k'; CutOff.LineWidth=1.5; CutOff.HandleVisibility='off'
OrigPeak=line([EggLength(NoiseContributions(1,4)) EggLength(NoiseContributions(1,4))],[-1 40]); OrigPeak.Color='g';
Dist32Peak=line([EggLength(NoiseContributions(7,4)) EggLength(NoiseContributions(7,4))],[-1 40]); Dist32Peak.Color='r';
legend('Inter-Allele noise', 'Co-Variance','Total Noise','RT peak','32C Peak','Location','best');
xlabel('% Egg length');
ylabel('% noise increase');
title('Kr Distal');

figure
errorbar(EggLength,nanmean(WholeNoise(1).TotalNoise),WholeNoise(1).SETotalNoise,'Color',Colors(1).Color,'LineWidth',2);
hold on 
errorbar(EggLength,nanmean(WholeNoise(7).TotalNoise),WholeNoise(7).SETotalNoise,'Color',Colors(1).Color,'LineWidth',2,'LineStyle','-.');
OrigPeak=line([EggLength(NoiseContributions(1,4)) EggLength(NoiseContributions(1,4))],[1 10]); OrigPeak.Color='g';
Dist32Peak=line([EggLength(NoiseContributions(7,4)) EggLength(NoiseContributions(7,4))],[1 10]); Dist32Peak.Color='r';
set(gca, 'FontSize',fontsize,'FontName',fontname);
xlim([0 100]);
xlabel('% egg length');
ylabel('total noise');
print( [FigDirect filesep 'DistTCompAPTotalNoise'],'-dsvg');

figure
c([1:16],1)=0; c([1:16],2)=0; c([1:16],3)=1;
c([17:25],1)=1; c([17:25],2)=0; c([17:25],3)=0;
c([26:41],1)=0; c([26:41],2)=1; c([26:41],3)=0;
c=c([1:41],[1:3]);
scatter(nanmean(WholeNoise(1).AllTotalmRNA),nanmean(WholeNoise(1).TotalNoise),[],c,'LineWidth',2);
hold on 
scatter(nanmean(WholeNoise(7).AllTotalmRNA), nanmean(WholeNoise(7).TotalNoise),[],c,'LineWidth',2,'Marker','d');
ylabel('Avg total noise');
xlabel('Avg total mRNA');
legend('Distal', 'Distal 32C');

plot((nanmean(WholeNoise(1).AllTotalmRNA)),nanmean(WholeNoise(1).TotalNoise),'o','Color',Colors(1).Color,'LineWidth',2);
hold on 
plot((nanmean(WholeNoise(7).AllTotalmRNA)),(nanmean(WholeNoise(1).TotalNoise)),'o','Color', Colors(7).Color,'LineWidth',2);
ylabel('Total Noise');
xlabel('Avg total mRNA');
%%
figure
PlaceHolder=0;
for cc=[2,8]
    PlaceHolder=PlaceHolder+1;
    PointsUsed=length(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    scatter((PlaceHolder.*ones(PointsUsed,1)),WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),'jitter','on','jitteramount',0.3);
    hold on 
    errorbar(PlaceHolder, nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4)),'Color','k','LineWidth',2.5);
    h=line([PlaceHolder-0.5,PlaceHolder+0.5],[nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)))]);
    h.Color='k'; h.LineWidth=2.5;
end
set(gca,'FontSize',fontsize,'FontName',fontname);
ylim([0 6]);
xticks([1:2]);
xticklabels({'proximal','proximal 32C'});
ylabel('total noise');

figure
NoiseStuff=[];
NoiseError=[];
Counter=0;
for cc=[2,8]
    Counter=Counter+1;
    NoiseError=[NoiseError; [NoiseContributions(cc,6),NoiseContributions(cc,5)]];
    NoiseStuff=[NoiseStuff;[NoiseContributions(cc,2),NoiseContributions(cc,1)]];
end
b=bar(NoiseStuff,'BarWidth',1,'FaceColor','flat');
PlaceCounter=0;
for cc=[2,8]
    PlaceCounter=PlaceCounter+1;
    b(1).CData(PlaceCounter,:)=[Colors(cc).Color];
    b(2).CData(PlaceCounter,:)=[Colors(cc).Color];
    y=nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    h=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); h.Color='k'; h.LineWidth=2;
end
ngroups=size(NoiseStuff,1);
groupwidth=min(0.8, 2/3.5);
hold on 
 for bb=1:2
     x=(1:ngroups)-groupwidth/2 + (2*bb-1) * groupwidth / 4;
     errorbar(x, NoiseStuff(:,bb), NoiseError(:,bb),'Color','k','LineWidth',2,'LineStyle','none');
 end
 xticklabels({'proximal', 'proximal 32C'});
 set(gca,'FontSize',fontsize,'FontName',fontname);
ylabel('noise');

figure
bar(NoiseContributions([2,8],[1:2]),'stacked')
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
xticklabels({'Proximal','Proximal 32C'});
hold on
%errorbar(cumsum(NoiseContributions([1,7],[1:2])')',NoiseContributions([1,7],[5:6]), 'Color','k','LineStyle','none');
errorbar(1,NoiseContributions(2,1),NoiseContributions(2,5),'Color','k','LineWidth',1.5);
errorbar(1, NoiseContributions(2,2),NoiseContributions(2,6),'Color','k','LineWidth',1.5);
errorbar(2,NoiseContributions(8,1),NoiseContributions(8,5),'Color','k','LineWidth',1.5);
errorbar(2, NoiseContributions(8,2),NoiseContributions(8,6),'Color','k','LineWidth',1.5);

figure 
scatter(WholeNoise(2).CorrTotNoise(:), WholeNoise(2).TotalNoise(:), [], 'g');
hold on 
scatter(WholeNoise(8).CorrTotNoise(:), WholeNoise(8).TotalNoise(:), [], 'r');
ylabel('Total Noise');
xlabel('Total mRNA production');
legend('Proximal', 'Proximal 32C');
%%
figure 
RTNoisevmRNA=[]; HotNoisevmRNA=[];
[sorted_x, sorted_x_index] = sort(WholeNoise(2).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(2).TotalNoise(sorted_x_index),70);
RTNoisevmRNA=[sorted_x, smooth_y];
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(2).Color,'LineWidth',2.5)

[sorted_x, sorted_x_index] = sort(WholeNoise(8).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(8).TotalNoise(sorted_x_index),70);
HotNoisevmRNA=[sorted_x, smooth_y];
plot(sorted_x, smooth_y, 'Color',Colors(2).Color,'LineWidth',1.5,'LineStyle','--')
ylabel('total noise');
xlabel('total mRNA production');
%legend('Proximal RT', 'Proximal 32C');
set(gca, 'FontSize',fontsize, 'FontName', fontname);
print( [FigDirect filesep 'NoisevExpProximal'],'-dsvg');
[H, pval,KStat]=kstest_2s_2d(RTNoisevmRNA,HotNoisevmRNA);
%%
figure
[ind1, ind2]=find(WholeNoise(2).AllCoVarNoise<0);
for aa=1:length(APbinID)
    NegFract(aa)=length(find(ind2==aa));
    NegFract(aa)=NegFract(aa)/(sum(~isnan(WholeNoise(2).AllCoVarNoise(:,aa))));
end
plot(EggLength, (NegFract.*100),'Color',Colors(2).Color, 'LineWidth',2.5);
ylabel('% of nuclei with negative co-variance');
xlabel('% egg length');
xlim([0 100])
hold on 
h=line([EggLength(NoiseContributions(2,4)) EggLength(NoiseContributions(2,4))], [0 20]);
h.Color='g', h.LineWidth=2;
set(gca, 'FontSize',fontsize,'FontName',fontname);

[ind1, ind2]=find(WholeNoise(8).AllCoVarNoise<0);
for aa=1:length(APbinID)
    NegFract(aa)=length(find(ind2==aa));
    NegFract(aa)=NegFract(aa)/(sum(~isnan(WholeNoise(8).AllCoVarNoise(:,aa))));
end
plot(EggLength, NegFract.*100, 'Color', Colors(8).Color,'LineWidth',2.5,'LineStyle','-.');
%%

PercProxIntra=(nanmean(WholeNoise(8).AllIntraNoise))-(nanmean(WholeNoise(2).AllIntraNoise));
PercProxIntra=(PercProxIntra./(nanmean(WholeNoise(2).AllIntraNoise))).*100;
PercProxCoVar=(nanmean(WholeNoise(8).AllCoVarNoise))-(nanmean(WholeNoise(2).AllCoVarNoise));
PercProxCoVar=(PercProxCoVar./(nanmean(WholeNoise(2).AllCoVarNoise))).*100;
PercProxTot=(nanmean(WholeNoise(8).TotalNoise))-(nanmean(WholeNoise(2).TotalNoise));
PercProxTot=(PercProxTot./(nanmean(WholeNoise(2).TotalNoise))).*100;
figure
plot(EggLength,PercProxIntra,'LineWidth',2);
hold on 
plot(EggLength,PercProxCoVar,'LineWidth',2);
plot(EggLength,PercProxTot,'LineWidth',2);
CutOff=refline(0,0); CutOff.Color='k'; CutOff.LineWidth=1.5; CutOff.HandleVisibility='off'
OrigPeak=line([EggLength(NoiseContributions(2,4)) EggLength(NoiseContributions(2,4))],[-1 40]); OrigPeak.Color='g';
Prox32Peak=line([EggLength(NoiseContributions(8,4)) EggLength(NoiseContributions(8,4))],[-1 40]); Prox32Peak.Color='r';
legend('Inter-Allele noise', 'Co-Variance','Total Noise','RT peak','32C Peak','Location','best');
xlabel('% egg length');
ylabel('% noise increase');
title('Kr Proximal');

figure
errorbar(EggLength,nanmean(WholeNoise(2).TotalNoise),WholeNoise(2).SETotalNoise,'Color',Colors(2).Color,'LineWidth',2);
hold on 
errorbar(EggLength,nanmean(WholeNoise(8).TotalNoise),WholeNoise(8).SETotalNoise,'Color',Colors(2).Color,'LineWidth',2,'LineStyle','-.');
OrigPeak=line([EggLength(NoiseContributions(2,4)) EggLength(NoiseContributions(2,4))],[1 10]); OrigPeak.Color='g';
Prox32Peak=line([EggLength(NoiseContributions(8,4)) EggLength(NoiseContributions(8,4))],[1 10]); Prox32Peak.Color='r';
set(gca, 'FontSize',fontsize,'FontName',fontname);
ylabel('total noise');
xlabel('% egg length');
print( [FigDirect filesep 'ProxTCompTotalNoiseAP'],'-dsvg');
%legend('Proximal', 'Proximal 32C');

%%
figure
NoiseStuff=[];
NoiseError=[];
Counter=0;
for cc=[3,9]
    Counter=Counter+1;
    NoiseError=[NoiseError; [NoiseContributions(cc,6),NoiseContributions(cc,5)]];
    NoiseStuff=[NoiseStuff;[NoiseContributions(cc,2),NoiseContributions(cc,1)]];
end
b=bar(NoiseStuff,'BarWidth',1,'FaceColor','flat');
PlaceCounter=0;
for cc=[3,9]
    PlaceCounter=PlaceCounter+1;
    b(1).CData(PlaceCounter,:)=[Colors(cc).Color];
    b(2).CData(PlaceCounter,:)=[Colors(cc).Color];
    y=nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    h=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); h.Color='k'; h.LineWidth=2;
end
ngroups=size(NoiseStuff,1);
groupwidth=min(0.8, 2/3.5);
hold on 
 for bb=1:2
     x=(1:ngroups)-groupwidth/2 + (2*bb-1) * groupwidth / 4;
     errorbar(x, NoiseStuff(:,bb), NoiseError(:,bb),'Color','k','LineWidth',2,'LineStyle','none');
 end
 xticklabels({'separate', 'separate 32C'});
 set(gca,'FontSize',fontsize,'FontName',fontname);
ylabel('noise');

figure
PlaceHolder=0;
for cc=[3,9]
    PlaceHolder=PlaceHolder+1;
    PointsUsed=length(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    scatter((PlaceHolder.*ones(PointsUsed,1)),WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),'jitter','on','jitteramount',0.3);
    hold on 
    errorbar(PlaceHolder, nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4)),'Color','k','LineWidth',2.5);
    h=line([PlaceHolder-0.5,PlaceHolder+0.5],[nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)))]);
    h.Color='k'; h.LineWidth=2.5;
end
set(gca,'FontSize',fontsize,'FontName',fontname);
ylim([0 6]);
xticks([1:2]);
xticklabels({'separate','separate 32C'});
ylabel('total noise');

figure
bar(NoiseContributions([3,9],[1:2]),'stacked')
title('Total noise max mRNA production AP bin');
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
xticklabels({'Separated','Separated 32C'});
hold on
%errorbar(cumsum(NoiseContributions([3,8],[1:2])')',NoiseContributions([3,8],[5:6]), 'Color','k','LineStyle','none');
errorbar(1,NoiseContributions(3,1),NoiseContributions(3,5),'Color','k','LineWidth',1.5);
errorbar(1, NoiseContributions(3,2),NoiseContributions(3,6),'Color','k','LineWidth',1.5);
errorbar(2,NoiseContributions(9,1),NoiseContributions(9,5),'Color','k','LineWidth',1.5);
errorbar(2, NoiseContributions(9,2),NoiseContributions(9,6),'Color','k','LineWidth',1.5);

figure 
scatter(WholeNoise(3).CorrTotNoise(:), WholeNoise(3).TotalNoise(:), [], 'g');
hold on 
scatter(WholeNoise(9).CorrTotNoise(:), WholeNoise(9).TotalNoise(:), [], 'r');
ylabel('Total Noise');
xlabel('Total mRNA production');
legend('Separate', 'Separate 32C');

figure 
[sorted_x, sorted_x_index] = sort(WholeNoise(3).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(3).TotalNoise(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'g','LineWidth',1.5)

[sorted_x, sorted_x_index] = sort(WholeNoise(9).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(9).TotalNoise(sorted_x_index),70);
plot(sorted_x, smooth_y, 'r','LineWidth',1.5)
ylabel('total noise');
xlabel('total mRNA production');
%legend('Separated RT', 'Separated 32C');
set(gca, 'FontSize',fontsize, 'FontName', fontname);
print( [FigDirect filesep 'NoisevExpSep'],'-dsvg');


PercSepIntra=(nanmean(WholeNoise(9).AllIntraNoise))-(nanmean(WholeNoise(3).AllIntraNoise))
PercSepIntra=(PercSepIntra./(nanmean(WholeNoise(3).AllIntraNoise))).*100;
PercSepCoVar=(nanmean(WholeNoise(9).AllCoVarNoise))-(nanmean(WholeNoise(3).AllCoVarNoise));
PercSepCoVar=(PercSepCoVar./(nanmean(WholeNoise(3).AllCoVarNoise))).*100;
PercSepTot=(nanmean(WholeNoise(9).TotalNoise))-(nanmean(WholeNoise(3).TotalNoise));
PercSepTot=(PercSepTot./(nanmean(WholeNoise(3).TotalNoise))).*100;
figure
plot(EggLength,PercSepIntra,'LineWidth',2);
hold on 
plot(EggLength,PercSepCoVar,'LineWidth',2);
plot(EggLength,PercSepTot, 'LineWidth',2);
CutOff=refline(0,0); CutOff.Color='k'; CutOff.LineWidth=1.5; CutOff.HandleVisibility='off'
OrigSepPeak=line([EggLength(NoiseContributions(3,4)) EggLength(NoiseContributions(3,4))],[-1 10]); OrigSepPeak.Color='g';
Sep32Peak=line([EggLength(NoiseContributions(9,4)) EggLength(NoiseContributions(9,4))],[-1 10]); Sep32Peak.Color='r';
legend('Inter-Allele noise', 'Co-Variance','Total noise','RT peak','32C Peak','Location','best');
xlabel('% Egg length');
ylabel('% noise increase');
title('Kr Separate');

figure
errorbar(EggLength,nanmean(WholeNoise(3).TotalNoise),WholeNoise(3).SETotalNoise,'Color',Colors(3).Color,'LineWidth',2);
hold on 
errorbar(EggLength,nanmean(WholeNoise(9).TotalNoise),WholeNoise(9).SETotalNoise,'Color',Colors(3).Color,'LineWidth',2,'LineStyle','-.');
%OrigSepPeak=line([EggLength(NoiseContributions(3,4)) EggLength(NoiseContributions(3,4))],[1 10]); OrigSepPeak.Color='g';
%Sep32Peak=line([EggLength(NoiseContributions(9,4)) EggLength(NoiseContributions(9,4))],[1 10]); Sep32Peak.Color='r';
%legend('Separate', 'Separate 32C');
set(gca, 'FontSize',fontsize,'FontName',fontname);
xlim([0 100]);
xlabel('% egg length');
ylabel('total noise');
%%
figure
PlaceHolder=0;
for cc=[6,10]
    PlaceHolder=PlaceHolder+1;
    PointsUsed=length(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    scatter((PlaceHolder.*ones(PointsUsed,1)),WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),'jitter','on','jitteramount',0.3);
    hold on 
    errorbar(PlaceHolder, nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4)),'Color','k','LineWidth',2.5);
    h=line([PlaceHolder-0.5,PlaceHolder+0.5],[nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)))]);
    h.Color='k'; h.LineWidth=2.5;
end
set(gca,'FontSize',fontsize,'FontName',fontname);
ylim([0 6]);
xticks([1:2]);
xticklabels({'both','both 32C'});
ylabel('total noise');

figure
NoiseStuff=[];
NoiseError=[];
Counter=0;
for cc=[6,10]
    Counter=Counter+1;
    NoiseError=[NoiseError; [NoiseContributions(cc,6),NoiseContributions(cc,5)]];
    NoiseStuff=[NoiseStuff;[NoiseContributions(cc,2),NoiseContributions(cc,1)]];
end
b=bar(NoiseStuff,'BarWidth',1,'FaceColor','flat');
PlaceCounter=0;
for cc=[6,10]
    PlaceCounter=PlaceCounter+1;
    b(1).CData(PlaceCounter,:)=[Colors(cc).Color];
    b(2).CData(PlaceCounter,:)=[Colors(cc).Color];
    y=nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    h=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); h.Color='k'; h.LineWidth=2;
end
ngroups=size(NoiseStuff,1);
groupwidth=min(0.8, 2/3.5);
hold on 
 for bb=1:2
     x=(1:ngroups)-groupwidth/2 + (2*bb-1) * groupwidth / 4;
     errorbar(x, NoiseStuff(:,bb), NoiseError(:,bb),'Color','k','LineWidth',2,'LineStyle','none');
 end
 xticklabels({'both', 'both 32C'});
 set(gca,'FontSize',fontsize,'FontName',fontname);
ylabel('noise');

figure
bar(NoiseContributions([6,10],[1:2]),'stacked')
title('Total noise max mRNA production AP bin');
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
xticklabels({'Both','Both 32C'});
hold on
%errorbar(cumsum(NoiseContributions([6,9],[1:2])')',NoiseContributions([6,9],[5:6]), 'Color','k','LineStyle','none');
errorbar(1,NoiseContributions(6,1),NoiseContributions(6,5),'Color','k','LineWidth',1.5);
errorbar(1, NoiseContributions(6,2),NoiseContributions(6,6),'Color','k','LineWidth',1.5);
errorbar(2,NoiseContributions(10,1),NoiseContributions(10,5),'Color','k','LineWidth',1.5);
errorbar(2, NoiseContributions(10,2),NoiseContributions(10,6),'Color','k','LineWidth',1.5);

figure 
scatter(WholeNoise(6).CorrTotNoise(:), WholeNoise(6).TotalNoise(:), [], 'g');
hold on 
scatter(WholeNoise(10).CorrTotNoise(:), WholeNoise(10).TotalNoise(:), [], 'r');
ylabel('Total Noise');
xlabel('Total mRNA production');
legend('Both', 'Both 32C');
%% Total noise vs exp 32C/RT
figure 
RTNoisevmRNA=[]; HotNoisevmRNA=[];
[sorted_x, sorted_x_index] = sort(WholeNoise(6).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(6).TotalNoise(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(6).Color,'LineWidth',2.5)
RTNoisevmRNA=[sorted_x, smooth_y];

[sorted_x, sorted_x_index] = sort(WholeNoise(10).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(10).TotalNoise(sorted_x_index),70);
HotNoisevmRNA=[sorted_x, smooth_y];
plot(sorted_x, smooth_y, 'Color',Colors(6).Color,'LineWidth',1.5,'LineStyle','--')
ylabel('total noise');
xlabel('total mRNA production');
%legend('Both RT', 'Both 32C');
set(gca, 'FontSize',fontsize, 'FontName', fontname);
print( [FigDirect filesep 'NoisevExpBoth'],'-dsvg');
[H, pval,KStat]=kstest_2s_2d(RTNoisevmRNA,HotNoisevmRNA);

%% Total noise v exp 17C/RT
figure 
[sorted_x, sorted_x_index] = sort(WholeNoise(6).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(6).TotalNoise(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'color',Colors(6).Color,'LineWidth',2.5)
RTNoisevmRNA=[sorted_x, smooth_y];

[sorted_x, sorted_x_index] = sort(WholeNoise(14).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(14).TotalNoise(sorted_x_index),70);
plot(sorted_x, smooth_y, 'color',Colors(6).Color,'LineWidth',1.5,'LineStyle','--')
ylabel('total noise');
xlabel('total mRNA production');
HotNoisevmRNA=[sorted_x,smooth_y];
%legend('Distal RT', 'Distal 32C');
set(gca, 'FontSize',fontsize, 'FontName', fontname);
print( [FigDirect filesep 'NoisevExp17CBoth'],'-dsvg');
[H, pval,KStat]=kstest_2s_2d(RTNoisevmRNA,HotNoisevmRNA);

%%
figure
PlaceCounter=1;
for cc=[6,10,14]%length(ConstructList)
    
    PointsUsed=length(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)));
    %plot((PlaceCounter.*ones(PointsUsed,1)),WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),'jitter
    scatter(((PlaceCounter-0.5).*ones(PointsUsed,1)),WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4)),20,Colors(cc).Color,'jitter','on','jitterAmount',0.3);
     hold on 
    scatter(((PlaceCounter+0.3).*ones(PointsUsed,1)),WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)),20,Colors(cc).Color,'x','jitter','on','jitterAmount',0.3);
    SevPerc(cc)=prctile(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),75);
    %errorbar((PlaceCounter-0.5), nanmedian(WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95CoVar(NoiseContributions(cc,4)),'.','Color','k','LineWidth',2);
    %errorbar((PlaceCounter+0.3), nanmedian(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95InterNoise(NoiseContributions(cc,4)),'.','Color','k','LineWidth',2);
    
    SevFivPercInter(cc)=prctile(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)),75);
    SevFivPercCoVar(cc)=prctile(WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4)),75);
    %errorbar(PlaceCounter, nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4)),'.','Color','k','LineWidth',2);
    %y=nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    y=nanmedian(WholeNoise(cc).AllCoVarNoise(:,NoiseContributions(cc,4)));
    y2=nanmedian(WholeNoise(cc).AllIntraNoise(:,NoiseContributions(cc,4)));
    h=line([PlaceCounter-0.7,PlaceCounter-0.3],[y,y]); h.Color='k'; h.LineWidth=2;
    h2=line([PlaceCounter+0.1,PlaceCounter+0.5],[y2,y2]); h2.Color='k'; h2.LineWidth=2;
    PlaceCounter=PlaceCounter+2.4;
    %PlaceCounter=PlaceCounter+0.6;
end
SevFivPercInter=max(SevFivPercInter); SevFivPercCoVar=max(SevFivPercCoVar);
if SevFivPercInter > SevFivPercCoVar
    SevFivPerc=(0.75*SevFivPercInter);
    %SevFivPerc=(0.1*SevFivPercInter);
else
    SevFivPerc=(0.75*SevFivPercCoVar);
    %SevFivPerc=(0.1*SevFivPercCoVar);
end
set(gca,'FontSize',fontsize,'FontName',fontname);
ylabel('noise');
ylim([0 SevFivPerc])
xlim([0 7])
xticks([1,3,6]);
xticklabels({'RT','32C','17C'});
%% Check fraction nuclei w negative co-variance 
figure
[ind1, ind2]=find(WholeNoise(6).AllCoVarNoise<0);
for aa=1:length(APbinID)
    NegFract(aa)=length(find(ind2==aa));
    NegFract(aa)=NegFract(aa)/(sum(~isnan(WholeNoise(6).AllCoVarNoise(:,aa))));
end
plot(EggLength, (NegFract.*100),'Color',Colors(6).Color, 'LineWidth',2.5);
ylabel('% of nuclei with negative co-variance');
xlabel('% egg length');
xlim([0 100])
hold on 
h=line([EggLength(NoiseContributions(6,4)) EggLength(NoiseContributions(6,4))], [0 20]);
h.LineWidth=2; h.Color='g';
set(gca, 'FontName',fontname,'FontSize',fontsize);

[ind1, ind2]=find(WholeNoise(10).AllCoVarNoise<0);
for aa=1:length(APbinID)
    NegFract(aa)=length(find(ind2==aa));
    NegFract(aa)=NegFract(aa)/(sum(~isnan(WholeNoise(10).AllCoVarNoise(:,aa))));
end
plot(EggLength, NegFract.*100, 'Color', Colors(6).Color,'LineWidth',2.5,'LineStyle','-.');
%%
PercBothIntra=(nanmean(WholeNoise(10).AllIntraNoise))-(nanmean(WholeNoise(6).AllIntraNoise))
PercBothIntra=(PercBothIntra./(nanmean(WholeNoise(6).AllIntraNoise))).*100;
PercBothCoVar=(nanmean(WholeNoise(10).AllCoVarNoise))-(nanmean(WholeNoise(6).AllCoVarNoise));
PercBothCoVar=(PercBothCoVar./(nanmean(WholeNoise(6).AllCoVarNoise))).*100;
PercBothTot=(nanmean(WholeNoise(10).TotalNoise))-(nanmean(WholeNoise(6).TotalNoise));
PercBothTot=(PercBothTot./(nanmean(WholeNoise(6).TotalNoise))).*100;
figure
plot(EggLength,PercBothIntra,'LineWidth',2);
hold on 
plot(EggLength,PercBothCoVar,'LineWidth',2);
plot(EggLength,PercBothTot,'LineWidth',2);
CutOff=refline(0,0); CutOff.Color='k'; CutOff.LineWidth=1.5; CutOff.HandleVisibility='off'
OrigBothPeak=line([EggLength(NoiseContributions(6,4)) EggLength(NoiseContributions(6,4))],[-1 10]); OrigBothPeak.Color='g';
Both32Peak=line([EggLength(NoiseContributions(10,4)) EggLength(NoiseContributions(10,4))],[-1 10]); Both32Peak.Color='r';
legend('Inter-Allele noise', 'Co-Variance','Total noise','RT peak','32C Peak','Location','best');
xlabel('% Egg length');
ylabel('% noise increase');
title('Kr Both');

figure
errorbar(EggLength,nanmean(WholeNoise(6).TotalNoise),WholeNoise(6).SETotalNoise,'Color',Colors(6).Color,'LineWidth',2);
hold on 
errorbar(EggLength,nanmean(WholeNoise(10).TotalNoise),WholeNoise(10).SETotalNoise,'Color',Colors(6).Color,'LineWidth',2,'LineStyle','-.');
OrigBothPeak=line([EggLength(NoiseContributions(6,4)) EggLength(NoiseContributions(6,4))],[1 10]); OrigBothPeak.Color='g';
Both32Peak=line([EggLength(NoiseContributions(10,4)) EggLength(NoiseContributions(10,4))],[1 10]); Both32Peak.Color='r';
%legend('Both RT', '32C');
set(gca, 'FontSize',fontsize,'FontName',fontname);
xlim([0 100]);
xlabel('% egg length');
ylabel('total noise');
print( [FigDirect filesep 'BothTCompTotalNoiseAP'],'-dsvg');
%%
figure
PlaceHolder=0;
for cc=[5,11]
    PlaceHolder=PlaceHolder+1;
    PointsUsed=length(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    scatter((PlaceHolder.*ones(PointsUsed,1)),WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)),'jitter','on','jitteramount',0.3);
    hold on 
    errorbar(PlaceHolder, nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),WholeNoise(cc).Conf95TotalNoise(NoiseContributions(cc,4)),'Color','k','LineWidth',2.5);
    h=line([PlaceHolder-0.5,PlaceHolder+0.5],[nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4))),nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)))]);
    h.Color='k'; h.LineWidth=2.5;
end
set(gca,'FontSize',fontsize,'FontName',fontname);
ylim([0 6]);
xticks([1:2]);
xticklabels({'2x proximal','2x proximal 32C'});
ylabel('total noise');

figure
NoiseStuff=[];
NoiseError=[];
Counter=0;
for cc=[5,11]
    Counter=Counter+1;
    NoiseError=[NoiseError; [NoiseContributions(cc,6),NoiseContributions(cc,5)]];
    NoiseStuff=[NoiseStuff;[NoiseContributions(cc,2),NoiseContributions(cc,1)]];
end
b=bar(NoiseStuff,'BarWidth',1,'FaceColor','flat');
PlaceCounter=0;
for cc=[5,11]
    PlaceCounter=PlaceCounter+1;
    b(1).CData(PlaceCounter,:)=[Colors(cc).Color];
    b(2).CData(PlaceCounter,:)=[Colors(cc).Color];
    y=nanmean(WholeNoise(cc).TotalNoise(:,NoiseContributions(cc,4)));
    h=line([PlaceCounter-0.5,PlaceCounter+0.5],[y,y]); h.Color='k'; h.LineWidth=2;
end
ngroups=size(NoiseStuff,1);
groupwidth=min(0.8, 2/3.5);
hold on 
 for bb=1:2
     x=(1:ngroups)-groupwidth/2 + (2*bb-1) * groupwidth / 4;
     errorbar(x, NoiseStuff(:,bb), NoiseError(:,bb),'Color','k','LineWidth',2,'LineStyle','none');
 end
 xticklabels({'2x proximal', '2x proximal 32C'});
 set(gca,'FontSize',fontsize,'FontName',fontname);
ylabel('noise');

figure
bar(NoiseContributions([5,11],[1:2]),'stacked')
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
xticklabels({'2xProx','2x Prox 32C'});
hold on
errorbar(1,NoiseContributions(5,1),NoiseContributions(5,5),'Color','k','LineWidth',1.5);
errorbar(1, NoiseContributions(5,2),NoiseContributions(5,6),'Color','k','LineWidth',1.5);
errorbar(2,NoiseContributions(11,1),NoiseContributions(11,5),'Color','k','LineWidth',1.5);
errorbar(2, NoiseContributions(11,2),NoiseContributions(11,6),'Color','k','LineWidth',1.5);
%%
figure
RTCutoff=prctile(WholeNoise(5).TotalNoise(:),75);
RTIndPts=(WholeNoise(5).TotalNoise(:)<=RTCutoff);
RTXPts=WholeNoise(5).CorrTotNoise(RTIndPts);
RTYPts=WholeNoise(5).TotalNoise(RTIndPts);
scatter(RTXPts, RTYPts);
hold on 
[n,c]=hist3([RTXPts, RTYPts]);
[C,h1]=contour(c{1},c{2},n','LineWidth',2);

figure
HotCutoff=prctile(WholeNoise(11).TotalNoise(:),75);
HotIndPts=(WholeNoise(11).TotalNoise(:)<=HotCutoff);
HotXPts=WholeNoise(11).CorrTotNoise(HotIndPts);
HotYPts=WholeNoise(11).TotalNoise(HotIndPts);
scatter(HotXPts,HotYPts);
hold on 
[n2,c2]=hist3([HotXPts, HotYPts]);
[C2,h2]=contour(c2{1},c2{2},n2','LineWidth',2);
%%
figure
RTNoisevmRNA=[]; HotNoisevmRNA=[];
[sorted_x, sorted_x_index] = sort(WholeNoise(5).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(5).TotalNoise(sorted_x_index),70);
RTNoisevmRNA=[sorted_x, smooth_y];
hold on 
plot(sorted_x, smooth_y, 'Color',Colors(5).Color,'LineWidth',2.5)
RT75=prctile(WholeNoise(5).TotalNoise(:),75);
RT75Idx=(WholeNoise(5).TotalNoise<=RT75);
RTX=WholeNoise(5).CorrTotNoise(RT75Idx);
RTY=WholeNoise(5).TotalNoise(RT75Idx);

[sorted_x, sorted_x_index] = sort(WholeNoise(11).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(11).TotalNoise(sorted_x_index),70);
Hot75=prctile(WholeNoise(11).TotalNoise(:),75);
Hot75Idx=(WholeNoise(11).TotalNoise<=Hot75);
Hot75Y=WholeNoise(11).TotalNoise(Hot75Idx);
Hot75X=WholeNoise(11).CorrTotNoise(Hot75Idx);
HotNoisevmRNA=[sorted_x, smooth_y];
plot(sorted_x, smooth_y, 'Color',Colors(5).Color,'LineWidth',1.5,'LineStyle','--')

ylabel('total noise');
xlabel('total mRNA production');
%legend('2x Proximal', '2x Proximal 32C');
set(gca, 'FontName',fontname, 'FontSize',fontsize);
print( [FigDirect filesep 'NoisevExp2xProx'],'-dsvg');
[H, pval,KStat]=kstest_2s_2d(RTNoisevmRNA,HotNoisevmRNA);

%%
%% 2xDistal 
figure 
[sorted_x, sorted_x_index] = sort(WholeNoise(4).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(4).TotalNoise(sorted_x_index),70);
hold on 
plot(sorted_x, smooth_y, 'color',Colors(4).Color,'LineWidth',2.5)
RTNoisevmRNA=[sorted_x, smooth_y];

[sorted_x, sorted_x_index] = sort(WholeNoise(13).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(13).TotalNoise(sorted_x_index),70);
plot(sorted_x, smooth_y, 'color',Colors(13).Color,'LineWidth',1.5,'LineStyle','--')
ylabel('total noise');
xlabel('total mRNA production');
HotNoisevmRNA=[sorted_x,smooth_y];
%legend('Distal RT', 'Distal 32C');
set(gca, 'FontSize',fontsize, 'FontName', fontname);
print( [FigDirect filesep 'NoisevExp2xDistal'],'-dsvg');
[H, pval,KStat]=kstest_2s_2d(RTNoisevmRNA,HotNoisevmRNA);
%% Fraction of nuclei w negative co-variance by AP bin 
figure
[ind1, ind2]=find(WholeNoise(5).AllCoVarNoise<0);
for aa=1:length(APbinID)
    NegFract(aa)=length(find(ind2==aa));
    NegFract(aa)=NegFract(aa)/(sum(~isnan(WholeNoise(5).AllCoVarNoise(:,aa))));
end
plot(EggLength, (NegFract.*100),'Color',Colors(5).Color, 'LineWidth',2.5);
ylabel('% of nuclei with negative co-variance');
xlabel('% egg length');
xlim([0 100])
h=line([EggLength(NoiseContributions(5,4)) EggLength(NoiseContributions(5,4))], [0 20]);
h.LineWidth=2; h.Color='g';
set(gca, 'FontName',fontname,'FontSize',fontsize);

hold on 
[ind1, ind2]=find(WholeNoise(11).AllCoVarNoise<0);
for aa=1:length(APbinID)
    NegFract(aa)=length(find(ind2==aa));
    NegFract(aa)=NegFract(aa)/(sum(~isnan(WholeNoise(11).AllCoVarNoise(:,aa))));
end
plot(EggLength, NegFract.*100, 'Color', Colors(11).Color,'LineWidth',2.5,'LineStyle','-.');

[sorted_x, sorted_x_index] = sort(WholeNoise(5).CorrTotNoise(:));
smooth_y = smooth(WholeNoise(5).AllCoVarNoise(sorted_x_index),70);

%%


Perc2ProxIntra=(nanmean(WholeNoise(11).AllIntraNoise))-(nanmean(WholeNoise(5).AllIntraNoise))
Perc2ProxIntra=(Perc2ProxIntra./(nanmean(WholeNoise(5).AllIntraNoise))).*100;
Perc2ProxCoVar=(nanmean(WholeNoise(11).AllCoVarNoise))-(nanmean(WholeNoise(5).AllCoVarNoise));
Perc2ProxCoVar=(Perc2ProxCoVar./(nanmean(WholeNoise(5).AllCoVarNoise))).*100;
Perc2ProxTot=(nanmean(WholeNoise(11).TotalNoise))-(nanmean(WholeNoise(5).TotalNoise));
Perc2ProxTot=(Perc2ProxTot./(nanmean(WholeNoise(5).TotalNoise))).*100;
figure
plot(EggLength,Perc2ProxIntra,'LineWidth',2);
hold on 
plot(EggLength,Perc2ProxCoVar,'LineWidth',2);
plot(EggLength,Perc2ProxTot,'LineWidth',2);
CutOff=refline(0,0); CutOff.Color='k'; CutOff.LineWidth=1.5; CutOff.HandleVisibility='off'
Orig2ProxPeak=line([EggLength(NoiseContributions(5,4)) EggLength(NoiseContributions(5,4))],[-1 10]); Orig2ProxPeak.Color='g';
DoubProx32Peak=line([EggLength(NoiseContributions(11,4)) EggLength(NoiseContributions(11,4))],[-1 10]); DoubProx32Peak.Color='r';
legend('Inter-Allele noise', 'Co-Variance','Total noise','RT peak','32C Peak','Location','best');
xlabel('% Egg length');
ylabel('% noise increase');
title('Kr 2x Proximal');

figure
errorbar(EggLength,nanmean(WholeNoise(5).TotalNoise),WholeNoise(5).SETotalNoise,'Color',Colors(5).Color,'LineWidth',2);
hold on 
errorbar(EggLength,nanmean(WholeNoise(11).TotalNoise),WholeNoise(11).SETotalNoise,'Color',Colors(5).Color,'LineWidth',2,'LineStyle','-.');
Orig2ProxPeak=line([EggLength(NoiseContributions(5,4)) EggLength(NoiseContributions(5,4))],[1 10]); Orig2ProxPeak.Color='g';
DoubProx32Peak=line([EggLength(NoiseContributions(11,4)) EggLength(NoiseContributions(11,4))],[1 10]); DoubProx32Peak.Color='r';
RTPerc=prctile(WholeNoise(5).TotalNoise,75); HotPerc=prctile(WholeNoise(11).TotalNoise,75);
if (RTPerc(21)) > (HotPerc(21))
    SevFivPerc=RTPerc(21)
else
    SevFivPerc=HotPerc(21)
end
set(gca, 'FontSize',fontsize,'FontName',fontname);
xlim([0 100]);
%ylim([0 SevFivPerc]);
xlabel('% egg length');
ylabel('total noise');
print( [FigDirect filesep '2xProxTCompTotalNoiseAP'],'-dsvg');
%%
figure
bar(PercentNoiseContributions([1,7],[1:2]),'stacked')
title('Noise contribution at max mRNA production AP bin');
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
xticklabels({'Distal','Distal 32C'});
hold on
    errorbar(1,PercentNoiseContributions(1,1),(PercentNoiseContributions(1,5))*100,'Color','k','LineWidth',1.5);
    errorbar(1,PercentNoiseContributions(1,2),(PercentNoiseContributions(1,6))*100,'Color','k','LineWidth',1.5);
errorbar(2,PercentNoiseContributions(7,1),(PercentNoiseContributions(7,5))*100,'Color','k','LineWidth',1.5);
    errorbar(2,PercentNoiseContributions(7,2),(PercentNoiseContributions(7,6))*100,'Color','k','LineWidth',1.5);

figure
bar(PercentNoiseContributions([2,8],[1:2]),'stacked')
title('Noise contribution at max mRNA production AP bin');
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
xticklabels({'Proximal','Proximal 32C'});
hold on
    errorbar(1,PercentNoiseContributions(2,1),(PercentNoiseContributions(2,5))*100,'Color','k','LineWidth',1.5);
    errorbar(1,PercentNoiseContributions(2,2),(PercentNoiseContributions(2,6))*100,'Color','k','LineWidth',1.5);
errorbar(2,PercentNoiseContributions(8,1),(PercentNoiseContributions(8,5))*100,'Color','k','LineWidth',1.5);
    errorbar(2,PercentNoiseContributions(8,2),(PercentNoiseContributions(8,6))*100,'Color','k','LineWidth',1.5);
    
    
figure
bar(PercentNoiseContributions([3,9],[1:2]),'stacked')
title('Noise contribution at max mRNA production AP bin');
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
xticklabels({'Separated','Separated 32C'});
hold on
errorbar(1,PercentNoiseContributions(3,1),(PercentNoiseContributions(3,5))*100,'Color','k','LineWidth',1.5);
    errorbar(1,PercentNoiseContributions(3,2),(PercentNoiseContributions(3,6))*100,'Color','k','LineWidth',1.5);
errorbar(2,PercentNoiseContributions(9,1),(PercentNoiseContributions(9,5))*100,'Color','k','LineWidth',1.5);
    errorbar(2,PercentNoiseContributions(9,2),(PercentNoiseContributions(9,6))*100,'Color','k','LineWidth',1.5);

figure
bar(PercentNoiseContributions([6,10],[1:2]),'stacked')
title('Noise contribution at max mRNA production AP bin');
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
xticklabels({'Both','Both 32C'});
hold on
errorbar(1,PercentNoiseContributions(6,1),(PercentNoiseContributions(6,5))*100,'Color','k','LineWidth',1.5);
    errorbar(1,PercentNoiseContributions(6,2),(PercentNoiseContributions(6,6))*100,'Color','k','LineWidth',1.5);
errorbar(2,PercentNoiseContributions(10,1),(PercentNoiseContributions(10,5))*100,'Color','k','LineWidth',1.5);
    errorbar(2,PercentNoiseContributions(10,2),(PercentNoiseContributions(10,6))*100,'Color','k','LineWidth',1.5);

figure
bar(PercentNoiseContributions([5,11],[1:2]),'stacked')
title('Noise contribution at max mRNA production AP bin');
legend('inter-allele','co-variance','Location','best','AutoUpdate','off')
ylabel('Noise (A.U.)')
xticklabels({'2x Prox','2xProx 32C'});
hold on
errorbar(1,PercentNoiseContributions(5,1),(PercentNoiseContributions(5,5))*100,'Color','k','LineWidth',1.5);
    errorbar(1,PercentNoiseContributions(5,2),(PercentNoiseContributions(5,6))*100,'Color','k','LineWidth',1.5);
errorbar(2,PercentNoiseContributions(11,1),(PercentNoiseContributions(11,5))*100,'Color','k','LineWidth',1.5);
    errorbar(2,PercentNoiseContributions(11,2),(PercentNoiseContributions(11,6))*100,'Color','k','LineWidth',1.5);

PlotbyAP=input('Make graphs across AP?','s');
if PlotbyAP=='y'
%Plot points of each noise at given AP position
APbinUse=input('Which AP bin to use?')
for cc=1:length(ConstructList)
    plot(cc,nanmean(WholeNoise(cc).TotalNoise(:,APbinUse)),'o','Color','g','LineWidth',2.5);
    hold on
    plot(cc,nanmean(WholeNoise(cc).AllNoise(:,APbinUse)),'o','Color','b','LineWidth',2.5);
    plot(cc,nanmean(WholeNoise(cc).AllCoVarNoise(:,APbinUse)),'o','Color','y','LineWidth',2.5);
end
title([num2str(EggLength(APbinUse)),'% Egg length']);
xticks([1:length(ConstructList)]);
xticklabels({'Dist','Prox','BothSep','2xProx','2xDist','Both'});
%%
%Plot mean intra-allele noise across all nuclei for each construct
figure 
for cc=1:length(ConstructList)
    plot(EggLength, nanmean(WholeNoise(cc).TotalNoise),'Color',Colors(cc).Color,'LineWidth',1.5);
    hold on 
end
xlabel('% Egg length')
ylabel('Avg inter-allele noise');

%singles 
figure 
plot(EggLength, nanmean(WholeNoise(1).AllNoise),'Color',Colors(1).Color,'LineWidth',1.5);
hold on 
plot(EggLength, nanmean(WholeNoise(2).AllNoise),'Color',Colors(2).Color,'LineWidth',1.5);
plot(EggLength, nanmean(WholeNoise(6).AllNoise),'Color',Colors(6).Color,'LineWidth',1.5);
xlabel('% Egg length')
xlim([0 100]);
ylabel('Avg inter-allele noise');

%doubles
figure
plot(EggLength, nanmean(WholeNoise(5).AllNoise),'Color',Colors(5).Color,'LineWidth',1.5);
hold on
plot(EggLength, nanmean(WholeNoise(6).AllNoise),'Color',Colors(6).Color,'LineWidth',1.5);
xlabel('% Egg length')
xlim([0 100]);
ylabel('Avg inter-allele noise');

figure
plot(EggLength, nanmean(WholeNoise(3).AllNoise),'Color',Colors(3).Color,'LineWidth',1.5);
hold on
plot(EggLength, nanmean(WholeNoise(6).AllNoise),'Color',Colors(6).Color,'LineWidth',1.5);
xlabel('% Egg length')
xlim([0 100]);
ylabel('Avg inter-allele noise');

% plot Co-Variance noise
figure 
for cc=1:length(ConstructList)
    plot(EggLength, nanmean(WholeNoise(cc).AllCoVarNoise), 'Color', Colors(cc).Color, 'LineWidth', 1.5);
    hold on 
end
xlabel('% Egg length')
xlim([0 100]);
ylabel('Avg allele covariance');

% singles 
figure 
plot(EggLength, nanmean(WholeNoise(1).AllCoVarNoise), 'Color', Colors(1).Color, 'LineWidth', 1.5);
hold on 
plot(EggLength, nanmean(WholeNoise(2).AllCoVarNoise), 'Color', Colors(2).Color, 'LineWidth', 1.5);
xlabel('% Egg length')
xlim([0 100]);
ylabel('Avg allele covariance');

figure 
plot(EggLength, nanmean(WholeNoise(1).AllCoVarNoise), 'Color', Colors(1).Color, 'LineWidth', 1.5);
hold on 
plot(EggLength, nanmean(WholeNoise(2).AllCoVarNoise), 'Color', Colors(2).Color, 'LineWidth', 1.5);
plot(EggLength, nanmean(WholeNoise(6).AllCoVarNoise), 'Color', Colors(6).Color, 'LineWidth', 1.5);
xlabel('% Egg length')
xlim([0 100]);
ylabel('Avg allele covariance');



for cc=1:length(ConstructList)
    figure
    plot(EggLength,nanmean(WholeNoise(cc).AllIntraNoise),'LineWidth',1.5);
    hold on
    plot(EggLength,nanmean(WholeNoise(cc).AllCoVarNoise),'LineWidth',1.5);
    plot(EggLength,nanmean(WholeNoise(cc).TotalSystNoise),'LineWidth',1.5);
    legend('inter-allele noise','Co-variance','Total noise');
    title(ConstructList{cc});
    xlabel('% Egg length');
    xlim([0 100]);
end
end
save([DropboxFolder filesep 'Constructs' filesep 'AllInterNoise'],'WholeNoise');
%% Pairwise comparisions
ConsTotNoise=[];
ConsGroups=[];
for cc=[1,2,3,4,5,6];
    ConsTotNoise=[[ConsTotNoise]; [WholeNoise(cc).TotalNoise(:)]];
    for bb=1:length(WholeNoise(cc).TotalNoise(:))
    ConsGroups=[[ConsGroups]; cc];
    end
end
[p,tbl,stats]=kruskalwallis(ConsTotNoise,ConsGroups);
[c,m]=multcompare(stats,'CType','bonferroni');
xlswrite([DropboxFolder filesep 'Constructs' filesep 'TotNoiseMedianPairwise'],c);

ConsCoVar=[];
ConsCoVarGroups=[];
for cc=[1,2,3,4,5,6];
    ConsCoVar=[[ConsCoVar]; [WholeNoise(cc).AllCoVarNoise(:)]];
    for bb=1:length(WholeNoise(cc).AllCoVarNoise(:))
    ConsCoVarGroups=[[ConsCoVarGroups]; cc];
    end
end
[p,tbl,stats]=kruskalwallis(ConsCoVar,ConsCoVarGroups);
[c,m]=multcompare(stats,'CType','bonferroni');
xlswrite([DropboxFolder filesep 'Constructs' filesep 'CoVarMedianPairwise'],c);

ConsInterAllele=[];
ConsInterAlleleGroups=[];
for cc=[1,2,3,4,5,6];
    ConsInterAllele=[[ConsInterAllele]; [WholeNoise(cc).AllIntraNoise(:)]];
    for bb=1:length(WholeNoise(cc).AllIntraNoise(:))
    ConsInterAlleleGroups=[[ConsInterAlleleGroups]; cc];
    end
end
[p,tbl,stats]=kruskalwallis(ConsInterAllele,ConsInterAlleleGroups);
[c,m]=multcompare(stats,'CType','bonferroni');
xlswrite([DropboxFolder filesep 'Constructs' filesep 'InterAlleleMedianPairwise'],c);
%% Comparing Embryos 
EmbryoComp=[];
for ee=1:length(WholeNoise(4).Embryo)
    EmbryoComp=[EmbryoComp, WholeNoise(4).Embryo(ee).TotalNoiseCalc(:,22)];
end
    
%% ANOVA
DoAnova=input('Do ANOVA?','s');
if DoAnova=='y'
    
NoiseManovaVect=[];
APManovaVect=[];
ConManovaVect=[];
for cc=1:length(ConstructList)
    for bb=1:size(WholeNoise(cc).TotalNoise,1)
        for aa=1:length(APbinID)
            NoiseManovaVect=[NoiseManovaVect, WholeNoise(cc).TotalNoise(bb,aa)];
            APManovaVect=[APManovaVect, aa];
            ConManovaVect=[ConManovaVect,cc];
       
            
        end
    end
end

[p,tbl,stats]=anovan(NoiseManovaVect,{APManovaVect, ConManovaVect},'sstype',1,'varnames',{'AP bin','Construct'})
multcompare(stats,'Dimension',2)
end

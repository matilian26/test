%% Compare fraction of frames of each nucleus that have 0,1,2 spots for each construct
   %load constructs
ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Count for each construct
EvCounter=0;
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    Label = ConstructList(cc);
    AllNucFractZero=[];
    AllNucFractOne=[];
    AllNucFractTwo=[];
%Count for each embryo of construct 
Construct2Fract=[];
Construct0Fract=[];
Construct1Fract=[];
    for ee=1:NEmbryos
PrefixName=Data(ee).Prefix;
        filename=[DropboxFolder filesep PrefixName filesep 'APDetection.mat'];
        load(filename);
        filename2=[DropboxFolder filesep PrefixName filesep PrefixName '_lin.mat'];
        load(filename2);
        filename3=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        load(filename3);
        CompiledParticles=Data(ee).CompiledParticles;
        SpotLength(ee,cc)=length(SpotDiff);
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);
nuclei=unique([CompiledParticles_14.Nucleus]);
nc14Start=Data(ee).nc14;

for hh=1:length(schnitzcells)
    schnitzcells(hh).Nucleus=hh;
end
%nc14counter=0;                  %create schnitz structure only with spots that exist in nc14 (still has all of those nuclei's frames tho even if before nc14)
schnitzcells_nc14=schnitzcells;

for ff=1:length(schnitzcells_nc14)     %Now frames are only those in nc14
    schnitzcells_nc14(ff).frames=[schnitzcells_nc14(ff).frames(schnitzcells_nc14(ff).frames >= nc14Start)];
end

for ff=1:length(schnitzcells_nc14)  %used for histogram to see distribution of life in nc14- looks like 15 frames is good cutoff
    %frameslife(ff,ee,cc)=length(schnitzcells_nc14(ff).frames);
    if length(schnitzcells_nc14(ff).frames) >=15
        schnitzcells_nc14(ff)=schnitzcells_nc14(ff);   %limiting to nuclei with >=15 frames in nc14
    else
        %change nucs with too short life to nan's in schnitz and SpotDiff
        tempschnit=find([SpotDiff.Nucleus]==schnitzcells_nc14(ff).Nucleus);
        if ~isempty(tempschnit)
        SpotDiff(tempschnit(1)).Nucleus=nan;
        if length(tempschnit)==2
            SpotDiff(tempschnit(2)).Nucleus=nan;
        end
        end
        schnitzcells_nc14(ff).Nucleus=nan;
    end
end
%first translate x,y coordinates to AP position 
%Angle between the x-axis and the AP-axis
    if exist('coordPZoom', 'var')
        APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
    else
        error('coordPZoom not defined. Was AddParticlePosition.m run?')
    end
    APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
NucPosition=[];



for i=1:length(schnitzcells_nc14)     %find the AP position of each nucleus across time 
        for j=1:length(schnitzcells_nc14(i).frames)
            
            %Angle between the x-axis and the particle using the A position as a
            %zero
            Angles=atan((schnitzcells_nc14(i).ceny(j)-coordAZoom(2))./(schnitzcells_nc14(i).cenx(j)-coordAZoom(1)));
            
            %Distance between the points and the A point
            Distances=sqrt((coordAZoom(2)-schnitzcells_nc14(i).ceny(j)).^2+(coordAZoom(1)-schnitzcells_nc14(i).cenx(j)).^2);
            APPositions=Distances.*cos(Angles-APAngle);
            NucTimeFrame=schnitzcells_nc14(i).frames(j);  %making the columns the time frame, each row is a nucleus 
            NucPosition(i,NucTimeFrame)=APPositions/APLength;
            
        end
   %schnitzcells_nc14(i).APPos=[NucPosition(i,:)];     
end
for ii=1:length(NucPosition)
    schnitzcells_nc14(ii).APPos=[NucPosition(ii,:)];
end
for ii=1:length(schnitzcells_nc14)
    schnitzcells_nc14(ii).MeanAP=nanmean(schnitzcells_nc14(ii).APPos);
end
APbinID=[Data(ee).APbinID];                %estimate AP bin for later 
APEstm=[schnitzcells_nc14.MeanAP];
for j=1:length(schnitzcells_nc14)
    APEstm(j)=round(schnitzcells_nc14(j).MeanAP,2);
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            schnitzcells_nc14(j).APBin=APbinID(jj);
            break;
        end
    end
    if isempty(schnitzcells_nc14(j).frames)
        schnitzcells_nc14(j).APBin=nan;
    end
end

% find nuclei that never have spots 
NoSpotNucs=schnitzcells_nc14;          

for ii=1:length(nuclei)
    NucwSpot=find([NoSpotNucs.Nucleus]==nuclei(ii));
    if ~isempty(NucwSpot)
%         if length(NucwSpot) >= 2
%         NoSpotNucs(NucwSpot(1)).Nucleus=nan;    %change nuclei that exist in CompiledParticles(have spots in nc14) to nan
%         NoSpotNucs(NucwSpot(2)).Nucleus=nan;
%         else
            NoSpotNucs(NucwSpot).Nucleus=nan;
        
    end
end


for ii=1:length(NoSpotNucs)
    if ~isnan(NoSpotNucs(ii).Nucleus) & (~isempty(NoSpotNucs(ii).frames)) 
        SumNeverSpots(ii,1)=1;
        SumNeverSpots(ii,2)=NoSpotNucs(ii).Nucleus; %2nd column is nucleus   
        SumNeverSpots(ii,3)=NoSpotNucs(ii).APBin;
    end
end
NucFraction=[];     %rows are length of that compiledParticle life, columns are 0,1,or 2 spots, and pages are different nuclei
CalcNuc0Fract=[];
CalcNuc1Fract=[];
CalcNuc2Fract=[];
%See for each nucleus the fraction of frames with one vs two spots
for ii=1:length(SpotDiff)
    if (~isnan(SpotDiff(ii).Nucleus)) & (~isempty(SpotDiff(ii).Nucleus)) %want to only use nuc with life >=15 which above changed to nan if didn't meet
    for jj=1:length(SpotDiff(ii).SpotOne)   
        if (any(SpotDiff(ii).SpotOne(jj))) & (~isnan(SpotDiff(ii).SpotOne(jj))) & (any(SpotDiff(ii).SpotTwo(jj)))  & (~isnan(SpotDiff(ii).SpotTwo(jj)))
           NucFraction(jj, 3, ii)= 1;
         elseif (any(SpotDiff(ii).SpotOne(jj))) & (~isnan(SpotDiff(ii).SpotOne(jj))) & (SpotDiff(ii).SpotTwo(jj)==0) %| (isnan(SpotDiff(ii).SpotTwo(jj))))
             NucFraction(jj,2,ii)=1;
        elseif (SpotDiff(ii).SpotOne(jj) ==0) & (any(SpotDiff(ii).SpotTwo(jj))) & (~isnan(SpotDiff(ii).SpotTwo(jj)))
            NucFraction(jj,2,ii)=1;
         elseif (SpotDiff(ii).SpotOne(jj)==0) & (SpotDiff(ii).SpotTwo(jj)==0)
             NucFraction(jj,1,ii)=1;
             %If that frame is nan in SpotDiff means nucleus didn't exist
             %so making that nucleus' value nan at that frame 
%          elseif (isnan(SpotDiff(ii).SpotOne(jj))) & (isnan(SpotDiff(ii).SpotTwo(jj)))
%              NucFraction(jj,:,ii)=nan;
        end
    end
    elseif isempty(SpotDiff(ii).Nucleus) | isnan(SpotDiff(ii).Nucleus)
        NucFraction(:,:,ii)=nan;
    end
    
end
for pp=1:size(NucFraction,3)
    if (~isnan(SpotDiff(pp).Nucleus)) & (~isempty(SpotDiff(pp).Nucleus)) 
        Schnit14Nuc=find([schnitzcells_nc14.Nucleus]==SpotDiff(pp).Nucleus);
       if ~isempty(Schnit14Nuc)
        SpotDiff(pp).NoFract=(sum(NucFraction(:,1,pp)))/(length(schnitzcells_nc14(Schnit14Nuc(1)).frames));  % divide number of frames each nucleus has 0 spots by number of frames that nucleus exists for
        SpotDiff(pp).OneFract=(sum(NucFraction(:,2,pp)))/(length([schnitzcells_nc14(Schnit14Nuc(1)).frames]));%SpotDiff(pp).Nucleus).frames]);
        SpotDiff(pp).TwoFract=(sum(NucFraction(:,3,pp)))/(length([schnitzcells_nc14(Schnit14Nuc(1)).frames]));%([NoSpotNucs(SpotDiff(pp).Nucleus).frames]);
       end
    end 
end
    for pp=1:length(SpotDiff)       %3/7/18 still need to figure out why occassionally get fractions >1
        if SpotDiff(pp).NoFract >= 1
            SpotDiff(pp).NoFract=nan;
        end
        if SpotDiff(pp).OneFract >=1
            SpotDiff(pp).OneFract=nan;
        end
        if SpotDiff(pp).TwoFract >=1
            SpotDiff(pp).TwoFract=nan;
        end
%         if ((SpotDiff(pp).NoFract)+(SpotDiff(pp).OneFract)+(SpotDiff(pp).TwoFract)) ~= 1
%             ProblemSpot(ee,cc,pp)=1;
%              SpotDiff(pp).NoFract=nan;
%              SpotDiff(pp).OneFract=nan;
%             SpotDiff(pp).TwoFract=nan;
%         end
    end
    
%     CalcNuc0Fract(CalcNuc0Fract > 1)=nan;
%     CalcNuc1Fract(CalcNuc1Fract > 1)=nan;        %3/6/18 figuring out why sometimes more spot occurances than frames, so for now just remove when the fraction is greater than 1
%     CalcNuc2Fract(CalcNuc2Fract > 1)=nan;
% Construct2Fract=[Construct2Fract CalcNuc2Fract];
% Construct1Fract=[Construct1Fract CalcNuc1Fract];
% Construct0Fract=[Construct0Fract CalcNuc0Fract];
%for gg=1:length(SumNeverSpots)

for ii=1:length(SumNeverSpots)
    %spotcounter=spotcounter+ii;
    SpotDiff(end+1).NoFract=SumNeverSpots(ii,1);
    SpotDiff(end).OneFract=0;
    SpotDiff(end).TwoFract=0;
    SpotDiff(end).Nucleus=SumNeverSpots(ii,2);
    SpotDiff(end).APBin=SumNeverSpots(ii,3);
end

%CalcNuc0Fract=[CalcNuc0Fract SumNeverSpots2'];   %Add in the nuclei that never have a spot 
%end
%for ii=1:length(SpotDiff)
    Everything(cc).Embryo(ee).SpotDiff=SpotDiff;
%end

% for ii=1:length(SpotDiff)
%     EvCounter=EvCounter+1;
%     if (~isnan(SpotDiff(ii).Nucleus)) & (~isempty(SpotDiff(ii).Nucleus)) & (~isempty(SpotDiff(ii).APBin))
%             
%     Everything(EvCounter).Nucleus=SpotDiff(ii).Nucleus;  %need to say if Nucleus exists, etc
%      Everything(EvCounter).NoFract=SpotDiff(ii).NoFract;
%      Everything(EvCounter).OneFract=SpotDiff(ii).OneFract;
%      Everything(EvCounter).TwoFract=SpotDiff(ii).TwoFract;
%      Everything(EvCounter).APBin=SpotDiff(ii).APBin; 
%      Everything(EvCounter).Embryo=ee;
%      Everything(EvCounter).Construct=cc; %ConstructList{cc};
%         
%     else
%         Everything(EvCounter).Nucleus=nan;
%     end      
% end
%     Construct(ee).TwoSpot=[Construct2Fract];
%     Construct(ee).OneSpot=[Construct1Fract];
%     Construct(ee).NoSpot=[Construct0Fract];
          

%     AllFractSpots(cc).TwoSpot=[Construct(:).TwoSpot];%Construct2Fract;
%     AllFractSpots(cc).OneSpot=[Construct(:).OneSpot];
%     AllFractSpots(cc).NoSpot=[Construct(:).NoSpot];
%     for yy=1:length(AllFractSpots(cc).TwoSpot)
%         AllCons2Fract(yy,cc)=AllFractSpots(cc).TwoSpot(yy);
%     end
%     for yy=1:length(AllFractSpots(cc).OneSpot)
%         AllCons1Fract(yy,cc)=AllFractSpots(cc).OneSpot(yy);
%     end
%     for yy=1:length(AllFractSpots(cc).NoSpot)
%         AllCons0Fract(yy,cc)=AllFractSpots(cc).NoSpot(yy);
%     end
    

for ii=1:length(Everything(cc).Embryo(ee).SpotDiff)   %make empty APbin values nan's so indexing is right in later steps
    if isempty(Everything(cc).Embryo(ee).SpotDiff(ii).APBin)
        Everything(cc).Embryo(ee).SpotDiff(ii).APBin=nan;
    end
    if isempty(Everything(cc).Embryo(ee).SpotDiff(ii).TwoFract)
        Everything(cc).Embryo(ee).SpotDiff(ii).TwoFract=nan;
    end
    if isempty(Everything(cc).Embryo(ee).SpotDiff(ii).OneFract)
        Everything(cc).Embryo(ee).SpotDiff(ii).OneFract=nan;
    end
    if isempty(Everything(cc).Embryo(ee).SpotDiff(ii).NoFract)
        Everything(cc).Embryo(ee).SpotDiff(ii).NoFract=nan;
    end
%     if isempty(Everything(cc).Embryo(ee).SpotDiff(ii).Construct)   %same for construct for later on when use find on construct field
%         Everything(cc).Embryo(ee).SpotDiff(ii).Construct=nan;
%     end
end
end

%take the average of each spot fraction at each AP position for each
%construct
for zz=1:length([Everything(cc).Embryo])
    %for pp=1:length([Everything(cc).Embryo(zz).SpotDiff])
for ii=1:length(APbinID)  %maybe need to change empty 2/1/0 fraction fields to nans
    TempLocation=[];
    TempLocation=find([Everything(cc).Embryo(zz).SpotDiff.APBin]==APbinID(ii));%ismembertol([Everything.APBin],APbinID(ii))%find(abs([Everything.APBin]-APbinID(ii))<0.001);
    TwoSpotsVect=[]; OneSpotVect=[]; NoSpotVect=[];
    for jj=1:length(TempLocation)
        %if ~isempty(Everything(cc).Embryo(zz).SpotDiff(TempLocation(jj)).TwoFract)
        TwoSpotsVect(jj)=Everything(cc).Embryo(zz).SpotDiff(TempLocation(jj)).TwoFract;
        OneSpotVect(jj)=Everything(cc).Embryo(zz).SpotDiff(TempLocation(jj)).OneFract;
        NoSpotVect(jj)=Everything(cc).Embryo(zz).SpotDiff(TempLocation(jj)).NoFract;
        %end
    end
    AvgSpots(cc).Embryo(zz).MeanTwo(ii)=nanmean(TwoSpotsVect);  %column number is APbin
    AvgSpots(cc).Embryo(zz).SD(ii)=std(TwoSpotsVect);
    AvgSpots(cc).Embryo(zz).SE(ii)=AvgSpots(cc).Embryo(zz).SD(ii)/sqrt(length(TwoSpotsVect));
    AvgSpots(cc).Embryo(zz).MeanOne(ii)=nanmean(OneSpotVect);
    AvgSpots(cc).Embryo(zz).MeanNo(ii)=nanmean(NoSpotVect);
end
    %end
end

end

EggLength=APbinID .* 100;

%make bar graph for each construct of fraction of 1/2/0 spots at each AP
%bin

for cc=1%:length(ConstructList)
    figure %first take mean across all embryos 
    Construct2Avg=[];
    Construct1Avg=[];
    Construct0Avg=[];
    ConAvg=[];
    for ee=1:length([AvgSpots(cc).Embryo])
        Construct2Avg=[Construct2Avg;[AvgSpots(cc).Embryo(ee).MeanTwo]];
        Construct1Avg=[Construct1Avg;[AvgSpots(cc).Embryo(ee).MeanOne]];
        Construct0Avg=[Construct0Avg;[AvgSpots(cc).Embryo(ee).MeanNo]];
    end
    Construct2Avg=nanmean(Construct2Avg);
    Construct1Avg=nanmean(Construct1Avg);
    Construct0Avg=nanmean(Construct0Avg);
    ConAvg=[Construct2Avg', Construct1Avg', Construct0Avg'];
    bar(ConAvg,'stacked');
end
%end


%Now want to arrange by APbin (as the column and the fraction of whatever
%spot going down)
Every2Spots=nan(length(Everything),length(APbinID));
for ii=1:length(APbinID)     %right now skip AP bin 0 until decide how to deal with all 0's in structure
    APWant=find([Everything.APBin]==APbinID(ii));
    for jj=1:length(APWant)
        if ~isempty(Everything(APWant(jj)).TwoFract)
        Every2Spots(jj,ii)=Everything(APWant(jj)).TwoFract;
        end
    end
    
end
AllCons2Anova=[];
firsttime=1;
for uu=1:length(ConstructList)
    ConstructTwoSpots=nan(length([Everything.Construct]==uu),length(APbinID));
    ConstructOneSpot=nan(length([Everything.Construct]==uu),length(APbinID));
    ConstructNoSpots=nan(length([Everything.Construct]==uu),length(APbinID));


    Constructtemp=find([Everything.Construct]==uu);
    for kk=1:length(Constructtemp)
        ConstructSubset(kk)=Everything(Constructtemp(kk));
    end
    
    if ~isempty(AllCons2Anova)
        firsttime=0;
        for xx=1:length(AllCons2Anova)
            if ~isempty(ConstructSubset(xx).TwoFract)
            AllCons2Anova(xx,uu)=ConstructSubset(xx).TwoFract;
            end
        end
    end
    if firsttime==1
    for kk=1:length(ConstructSubset)                      %Trying to set up structure for anova by construct but need to limit so same amt data points for each structre- not sure this is the best way of doing so of data
        if ~isempty(ConstructSubset(kk).TwoFract)
        AllCons2Anova(kk,uu)=ConstructSubset(kk).TwoFract;
        end
    end
    end
for ii=1:length(APbinID)

    ConstructSpot=find([ConstructSubset.APBin]==APbinID(ii));
    for jj=1:length(ConstructSpot)
        if ~isempty(ConstructSubset(ConstructSpot(jj)).TwoFract)
            ConstructTwoSpots(jj,ii)=ConstructSubset(ConstructSpot(jj)).TwoFract;
        end
        if ~isempty(ConstructSubset(ConstructSpot(jj)).OneFract)
            ConstructOneSpot(jj,ii)=ConstructSubset(ConstructSpot(jj)).OneFract;
        end
        if ~isempty(ConstructSubset(ConstructSpot(jj)).NoFract)
            ConstructNoSpots(jj,ii)=ConstructSubset(ConstructSpot(jj)).NoFract;
        end
    end
end

%calculate SD and SE at each AP position for each construct
MeansStruct(uu).StdDev2=std(ConstructTwoSpots,'omitnan');
MeansStruct(uu).StdError2=(MeansStruct(uu).StdDev2)/(sqrt(length(ConstructTwoSpots)));   
MeansStruct(uu).StdDev1=std(ConstructOneSpot,'omitnan');
MeansStruct(uu).StdError1=(MeansStruct(uu).StdDev1)/(sqrt(length(ConstructOneSpot)));
MeansStruct(uu).StdDev0=std(ConstructNoSpots,'omitnan');
MeansStruct(uu).StdError0=(MeansStruct(uu).StdDev0)/(sqrt(length(ConstructNoSpots)));

ConstructMean2=nanmean(ConstructTwoSpots);      %take mean at each AP position so have row vector where each is mean of AP position
ConstructMean1=nanmean(ConstructOneSpot);
ConstructMeanNo=nanmean(ConstructNoSpots);
ConstructMean=[ConstructMeanNo',ConstructMean1',ConstructMean2'];     %Change so AP position going down 
figure 
bar(ConstructMean,'stacked');
xlabel('AP bin')
ylabel('Avg fraction of nucleus life');
legend('No Spot', '1 Spot', '2Spot');
title(ConstructList{uu});
ylim([0 1.1]);
xlim([0 41]);

MeansStruct(uu).Two=ConstructMean2;  %Mean fraction of life with each # of spots at each AP position 
MeansStruct(uu).One=ConstructMean1;
MeansStruct(uu).No=ConstructMeanNo;

%     [p,tbl,stats]=anova1(ConstructTwoSpots);
%     %title(ConstructList{uu});
%     %multcompare(stats);
%     xlim([-0.1 0.6])
%     ylim([1 28])
%     figure 
%     [p2,tbl2,stats2]=anova1(ConstructOneSpot);
%     xlim([-0.1 1]);
%     %multcompare(stats2)
%     figure
%     [p3,tbl3,stats3]=anova1(ConstructNoSpots)
%     %multcompare(stats3)
   clear ConstructMean ConstructMean1 ConstructMean2 ConstructMeanNo
end

%set up arrays for doing anovan analysis
TimewithTwo=[]; APPositionNuc=[]; ConstructValue=[];
for ii=1:length(Everything)
    if ~isempty(Everything(ii).TwoFract)
    TimewithTwo(ii)=Everything(ii).TwoFract;
    APPositionNuc(ii)=Everything(ii).APBin;
    ConstructValue(ii)=Everything(ii).Construct;
    end
end
[p,tbl,stats]=anovan(TimewithTwo,{APPositionNuc ConstructValue}, 'model','interaction','varnames',{'AP Position', 'Construct'},'sstype',2);
figure
for uu=1:length(ConstructList)
    TotalActivity(uu).activity=(MeansStruct(uu).Two) + (MeansStruct(uu).One);
    
    plot(EggLength,TotalActivity(uu).activity);
    hold on 
    %plot(EggLength,MeansStruct(uu).No);
    xlabel('% Egg Length')
    ylabel('Avg fraction nucleus life')
    title('Time with any transcriptional activity')
    legend('Both','Dist','Prox','BothSep','2xDist','2xProx','Location','best');
end
for uu=1:length(ConstructList)
    figure
    plot(EggLength,MeansStruct(uu).Two)
    hold on
    plot(EggLength,MeansStruct(uu).One)
    plot(EggLength,MeansStruct(uu).No)
    xlabel('% egg length')
    ylabel('Avg fraction of nucleus life');
    title(ConstructList{uu});
    legend('Two spots','One Spot','No Spots');
end
figure
for uu=1:length(ConstructList)
    plot(EggLength,MeansStruct(uu).Two)
    hold on 
end
ylim([0 1])
xlim([0 100]);
xlabel('% Egg Length')
ylabel('Avg fraction of nucleus life')
legend('KrBoth','Dist','Prox','BothSep','2xDist','2xProx');
title('Time with 2 spots');

figure
errorbar(EggLength,MeansStruct(4).Two,MeansStruct(4).StdError2,'LineWidth',1.5);
hold on 
errorbar(EggLength,MeansStruct(1).Two,MeansStruct(1).StdError2,'LineWidth',1.5);
ylim([0 1])
xlim([0 100]);
xlabel('% Egg Length')
ylabel('Mean fraction of nucleus life')
legend('KrBothSep','Both');
title('Time with both alleles active');

figure
errorbar(EggLength,MeansStruct(2).Two,MeansStruct(2).StdError2,'LineWidth',1.5)
hold on
errorbar(EggLength,MeansStruct(5).Two,MeansStruct(5).StdError2,'LineWidth',1.5)
ylim([0 1])
xlim([0 100]);
xlabel('% Egg Length')
ylabel('Mean fraction of nucleus life')
legend('Dist','2x Dist');
title('Time with both alleles active');

figure
errorbar(EggLength,MeansStruct(3).Two,MeansStruct(3).StdError2,'LineWidth',1.5)
hold on 
errorbar(EggLength,MeansStruct(6).Two,MeansStruct(6).StdError2,'LineWidth',1.5)
ylim([0 1])
xlim([0 100]);
xlabel('% Egg Length')
ylabel('Mean fraction of nucleus life')
legend('Prox','2x Prox');
title('Time with both alleles active');



figure
for uu=1:length(ConstructList)
    plot(EggLength,MeansStruct(uu).One)
    hold on 
end
ylim([0 1])
xlim([0 100]);
xlabel('% Egg Length')
ylabel('Avg fraction of nucleus life')
legend('KrBoth','Dist','Prox','BothSep','2xDist','2xProx');
title('Time with 1 spot');

figure
for uu=1:length(ConstructList)
    plot(EggLength,MeansStruct(uu).No)
    hold on 
end
ylim([0 1])
xlim([0 100]);
xlabel('% Egg Length')
ylabel('Avg fraction of nucleus life')
legend('KrBoth','Dist','Prox','BothSep','2xDist','2xProx');
title('Time with 0 spots');

% for uu=1:length(ConstructList)
%     figure
%     histogram(AllFractSpots(uu).TwoSpot,'FaceAlpha',0.3,'Normalization','probability');
%     hold on 
%     histogram(AllFractSpots(uu).OneSpot,'FaceAlpha',0.3,'Normalization','probability');
%     histogram(AllFractSpots(uu).NoSpot,'FaceAlpha',0.3,'Normalization','probability');
%     xlabel('Fraction of nucleus life with 2 spots');
%     ylabel('Number of Nuclei')
%     legend('TwoSpots','One','None','Location','best');
%     title(ConstructList{uu});
%     xlim([0 1.2])
%     ylim([0 1])
% end
% figure
% for uu=1:length(ConstructList)
%     %figure
%     histogram(AllFractSpots(uu).OneSpot);
%     hold on 
%     xlabel('Fraction of nucleus life with 1 spot');
%     ylabel('Number of Nuclei')
%     legend('Both', 'Dist', 'Prox', 'BothSep', '2xDist', '2xProx');
%     title(['1 spot']); %ConstructList{uu}]);
%     xlim([0 1.2]);
% end
% figure
% for uu=1:length(ConstructList)
%     histogram(AllFractSpots(uu).NoSpot);
%     hold on
%     xlabel('Fraction of nucleus life with 0 spots');
%     ylabel('Number of Nuclei')
%     legend('Both', 'Dist', 'Prox', 'BothSep', '2xDist', '2xProx','Location','best');
%     title(['0 spots']);
%     xlim([0 1.2]);
% end
% Embryos2Use=min(Datalength);
% SpotLength=SpotLength(1:Embryos2Use,:);
% Spots2Use=min(min(SpotLength));
% AllCons2Fract=AllCons2Fract(1:Spots2Use,:);
% [p,tbl,stats]=anova1(AllCons2Fract);
% multcompare(stats); 
% [p2,tbl2,stats2]=anova1(AllCons1Fract);
% [p3,tbl3,stats3]=anova1(AllCons0Fract);
% bar(ConstructFracts,'stacked');
% legend('Zero','One','Two');
% xlabel('Construct')
% ylabel('Mean Fraction of total frames');
% title('Mean fraction of frames with 0,1,or 2 spots per nucleus');

% [p,tbl,stats]=anova1(AnovaTwoSpots);
% [p2,tbl2,stats2]=anova1(EachNucAnovaTwoSpots);
    
% 
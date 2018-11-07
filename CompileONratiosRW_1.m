%% Look at fraction of 1, 2, or no spot nuclei vs AP position for all embryos of a construct
ConstructList= {'KrBoth';'KrBothSep';'KrDist';'KrProx';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
ConstructSpotsAllAP=[];
PrefixLabel={};
DispNuc=input('Fraction or number of nuclei? f/n','s');
j=input('10 or 20 frame threshold?', 's');
if j=='20'
for ii=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{ii});
    
    NEmbryos = length(Data);
    Label = ConstructList(ii);
   ConstructSpotsAP=[];
    for jj=1:NEmbryos
        PrefixName=Data(jj).Prefix;
        if DispNuc=='f'
        filename=[DropboxFolder filesep PrefixName filesep 'SpotsAPRatioFract.mat'];
        else
            filename=[DropboxFolder filesep PrefixName filesep 'SpotsAPRatio.mat'];
        end
        
        load(filename)
        ConstructSpotsAP(:,:,jj)=SpotsAP;
       
    end
    PrefixLabel{ii}=Label;
    ConstructSpotsAllAP(:,:,ii)=nanmean(ConstructSpotsAP,3);
    %ConstructSpotsAllAP(:,:,ii)=sum(ConstructSpotsAP,3);
end
end

 if j=='10'
        for ii=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{ii});
    
    NEmbryos = length(Data);
    Label = ConstructList(ii);
   ConstructSpotsAP=[];
    for jj=1:NEmbryos
        PrefixName=Data(jj).Prefix;
        filename=[DropboxFolder filesep PrefixName filesep 'SpotsAPRatio10Frame.mat'];
        
        load(filename)
        ConstructSpotsAP(:,:,jj)=SpotsAP;
       
    end
    PrefixLabel{ii}=Label;
    ConstructSpotsAllAP(:,:,ii)=nanmean(ConstructSpotsAP,3);
    %ConstructSpotsAllAP(:,:,ii)=sum(ConstructSpotsAP,3);
        end
    end


        
APbinID=Data(1).APbinID;

% BothSepSpotsAP=ConstructSpotsAllAP(:,:,1);
% DistSpotsAP=ConstructSpotsAllAP(:,:,2);
% ProxSpotsAP=ConstructSpotsAllAP(:,:,3);
for jj=1:length(ConstructList)
    figure
    bar(APbinID,ConstructSpotsAllAP(:,:,jj),'stacked');
    if DispNuc=='f'
        ylim([0 1.5]);
        ylabel('Fraction of nuclei');
    else
    ylim([0 40]);
    ylabel('Number of nuclei');
    end
    xlabel('Mean AP position'); 
    title(PrefixLabel{jj});
    legend('no spot','one spot','two spots','Location','best');
end
%% predict expected ON fraction for Kr BothSep 
ZerosConstruct=ConstructSpotsAllAP;
ZerosConstruct(ZerosConstruct==nan)=[];
Testtwo=[];
Test1ObsCorr=[];
Test1OtherCorr=[];
for ii=1:length(ConstructList)
    for jj=1:length(APbinID)                 %having nan's make all divisions go to nan if it is one of the values so need to change to 1 if exist
        allnuc(jj,ii)=nansum(ConstructSpotsAllAP(jj,:,ii));
%          
ProbON(ii).ON(jj)=(ConstructSpotsAllAP(jj,2,ii)+ConstructSpotsAllAP(jj,3,ii))./allnuc(jj,ii);
        ProbON(ii).NoSpot(jj)=ConstructSpotsAllAP(jj,1,ii)./allnuc(jj,ii);
         
        ProbON(ii).OneSpot(jj)=ConstructSpotsAllAP(jj,2,ii)./allnuc(jj,ii);
         %end
%       
Testtwo(jj,ii)=ProbON(ii).OneSpot(jj) * ProbON(ii).OneSpot(jj);
        ProbON(ii).TwoSpot(jj)=ConstructSpotsAllAP(jj,3,ii)./allnuc(jj,ii);
        CompProb(ii).OneSpot(jj)=1-((ProbON(ii).OneSpot(jj).^2) + (ProbON(ii).NoSpot(jj)).^2);
        CompProb(ii).TwoSpot(jj)=(ProbON(ii).OneSpot(jj)).^2;
        
        TestOne(jj,ii)=1-((ProbON(ii).OneSpot(jj).^2) + (ProbON(ii).NoSpot(jj)).^2);
        Test1ObsCorr=corrcoef(TestOne(:,ii),ProbON(ii).OneSpot(:),'Rows','complete');
        TestOneObsCorr(ii)=Test1ObsCorr(2);  % correlation between Probability of One spot calculated as fraction of 1 spot/total nuclei at AP to that calculated as 1-(P(one)^2 + P(OFF)^2)
        TestOneOther(jj,ii)=2.*(ProbON(ii).NoSpot(jj) .* ProbON(ii).OneSpot(jj));
        Test1OtherCorr=corrcoef(TestOneOther(:,ii),ProbON(ii).OneSpot(:),'Rows','complete');
        TestOneOtherCorr(ii)=Test1OtherCorr(2);
    FractProbNone(jj,ii)=ProbON(ii).NoSpot(jj);
    FractProbOne(jj,ii)=ProbON(ii).OneSpot(jj);
    FractProbTwo(jj,ii)=ProbON(ii).TwoSpot(jj);
    end
  
end
color2=['r','g','b'];
figure
for uu=1:length(ConstructList)
    color2=['r','g','b'];
  plot(APbinID,FractProbNone(:,uu),'o','Color',color2(uu));
  hold on 
end
legend('Kr both','Kr Dist','Kr Prox')
title('Probability No Spot')
figure
for uu=1:length(ConstructList)
    
    plot(APbinID,FractProbOne(:,uu),'Color',color2(uu));
     hold on 
end
for uu=1:length(ConstructList)
    plot(APbinID,CompProb(uu).OneSpot(:),'LineStyle','--','Color',color2(uu));
end

title('Probability One Spot by fraction')
legend('Kr both', 'Dist', 'Prox');

figure 
for uu=1:length(ConstructList)
    plot(APbinID,CompProb(uu).OneSpot(:),'Color',color2(uu))
    hold on 
end
title('Prob One Spot Comparision')



%Want to calculate predicted ON fractions of KRBothSep based on Dist/Prox
%probabilities 
for ii=1:length(ConstructList)
    for qq=1:length(APbinID)
        PredictedBoth2(qq)=ProbON(2).OneSpot(qq) .* ProbON(3).OneSpot(qq);
        PredictedBoth1(qq)=(ProbON(2).OneSpot(qq) + ProbON(3).OneSpot(qq)) - (ProbON(2).OneSpot(qq) .* ProbON(3).OneSpot(qq));
        PredictedBothZero(qq)=ProbON(2).NoSpot(qq) .* ProbON(3).NoSpot(qq);
    end
end
BothPredictions(:,1)=PredictedBothZero;
BothPredictions(:,2)=PredictedBoth1;
BothPredictions(:,3)=PredictedBoth2;
figure
bar(APbinID,BothPredictions,'stacked');
xlabel('Mean AP position')
title('Predicted Kr BothSep fraction ON nuclei')
legend('no spot','one spot','two spots','Location','best');

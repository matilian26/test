% %% plot mean correlation at each AP position for each embryo of a construct  (adj)
% ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
%     %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
%     
% [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
%  Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
%  ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
% 

%% Actually doing the comparison 
ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

ConCorrMean=[];
MeanAPCorr=[];
MeanAllCorr=[];
APCorr=[];


for ii=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{ii});
    
    NEmbryos = length(Data);
    Label = ConstructList(ii);
    figure
    APCorr=[];
    
    for jj=1:NEmbryos
        PrefixName=Data(jj).Prefix;
        filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        load(filename);
        APbinID=Data(1).APbinID;
        APCorr=nan(1,length(APBins));
        for aa=1:length(APbinID)
            APinfo=[SpotDiff.APBin];
for zz=1:length(Bins)
    temp3=find([APBins==Bins(zz)]);
    temp2=find([SpotDiff.APBin]==Bins(zz));%issue with find I think has to do with empty values 
    for qq=1:length(temp2)
        if ~isempty(SpotDiff(temp2(qq)).SpotCorr)
            if length(SpotDiff(temp2(qq)).SpotCorr)==1
                APCorr(qq,temp3,jj)=[SpotDiff(temp2(qq)).SpotCorr(1)];
            else
    APCorr(qq,temp3,jj)=[SpotDiff(temp2(qq)).SpotCorr(1,2)];
            end
        end
    end 
end
APCorr(APCorr==0)=nan;

MeanAPCorr(jj,:,ii)=nanmean(APCorr(:,:,jj));%mean of each embryo is a row and 3rd dimension is construct
MeanAllCorr(jj,ii)=nanmean(MeanAPCorr(jj,:,ii),2);
%setting up a structure to compare correlation and AP position (within each
%construct- working elsewhere on 2way ANOVA)

counter=0;
for tt=1:length(SpotDiff)
    APColumn=find([APBins==SpotDiff(tt).APBin]);
    if ~isempty([SpotDiff(tt).SpotCorr])
        counter=counter+1;
    if length(SpotDiff(tt).SpotCorr)==1
        SetupAP(counter,APColumn,jj)=SpotDiff(tt).SpotCorr(1)
    else
    SetupAP(counter,APColumn,jj)=SpotDiff(tt).SpotCorr(1,2);
    end
    %SetupAP(counter,APColumn,jj)=SpotDiff(tt).APBin;
    end
end
% ConstAPCorr=[];
%     for uu=1:length(Bins)
%         APColumn=find(APBins==Bins(uu));
% tempbins=find([SetupAP(:,2,:)]==Bins(uu));
% for zz=1:length(tempbins)
% ConstAPCor(zz,uu)=SetupAP(tempbins);
% end
%     end


%figure

    end
    
    ConstAPCor=[];
    for yy=1:NEmbryos
        ConstAPCor=vertcat(ConstAPCor, SetupAP(:,:,yy));
    end
ConstAPCor(ConstAPCor==0)=nan;

    %[p,tbl,stats]=anova1(ConstAPCor,'displayopt','off');
%     APStats(ii).p=p;
%     APStats(ii).tbl=tbl;
%     APStats(ii).stats=stats;
%     title(ConstructList{ii});
%     if ii >= 2
%         AllCorr(end:end+length(APCorr(:)),:)=nan;
%     end
    AllCorr(ii).AllValues=[APCorr(:)];
    %AllCorr(:,ii)=vertcat(AllCorr,[APCorr(:)]);

    ConstructAPCorr=[];
    for tt=1:NEmbryos
    ConstructAPCorr=vertcat(ConstructAPCorr,APCorr(:,:,tt));
    end
    %
    ConCorrMean(ii,:)=nanmean(ConstructAPCorr,1);
    
    h(ii)=plot(APBins, ConCorrMean(ii,:), 'LineWidth',2);
    hold on
    xlabel('Mean AP Position')
    ylabel('Average fluorescence correlation')
    Prefixname=ConstructList(ii)
    %legend([h(ii)],Prefixname);
    title(['2 Spot Correlation',Prefixname]);
    ylim([0 1]);
    xlim([0.4 0.75])
    
     
    
    clearvars -except APStats ConstructList DropboxFolder ConCorrMean AllCorr  APBins MeanAPCorr MeanAllCorr
end
for jj=1:length(ConstructList)
    NNucs(jj)=length(AllCorr(jj).AllValues);
    MaxNucs=max(NNucs);
end
AllCorrMat=zeros(MaxNucs,4);
for jj=1:length(ConstructList)
    for qq=1:length(AllCorr(jj).AllValues)
    AllCorrMat(qq,jj)=AllCorr(jj).AllValues(qq);
    end
end
AllCorrMat(AllCorrMat==0)=nan;  %ANOVA of spot correlation at all AP points across whole embryo for all embryos of a construct with columns being each different construct
[p,tbl,stats]=anova1(AllCorrMat);
title('Spot Correlation across whole embryo')
xlabel('construct')
figure
multcompare(stats);




%MeanAPCorr(MeanAPCorr==0)=nan;
% MeanAllCorr(MeanAllCorr==0)=nan;    %ANOVA of mean spot correlation across whole embryo with groups being the different constructs 
% Constructs=[ConstructList{1},'', ConstructList{2},'', ConstructList{3},'',ConstructList{4}];
% [p,tbl,stats]=anova1(MeanAllCorr)
% title('Mean Spot Correlation across whole embryo')
% xlabel('construct')
% figure
% multcompare(stats);

BothData=LoadMS2SetsCS(ConstructList{1});
DistData=LoadMS2SetsCS(ConstructList{2});
ProxData=LoadMS2SetsCS(ConstructList{3});
BothSepData=LoadMS2SetsCS(ConstructList{4});
DoubDistData=LoadMS2SetsCS(ConstructList{5});
DoubProxData=LoadMS2SetsCS(ConstructList{6});
maxNEmbryos=max(size(MeanAPCorr,1));
MeanAPCorrAll=[];
for ii=1:length(ConstructList)         %setting up array for 2way ANOVA
MeanAPCorrAll=vertcat(MeanAPCorrAll,MeanAPCorr(:,:,ii));  %right now just look at 3 embryos (since have that for all constructs, until figure out unbalanced anova)
end



figure
C={'r','g','b','y','c','m'};
for jj=1:length(ConstructList)  %want to plot construct means all on one graph for comparison 
    h(jj)=plot(APBins,ConCorrMean(jj,:),'LineWidth',2,'Color',C{jj});%plot(APBins,h(jj),'LineWidth',2)
    hold on
end
legend('Kr Both', 'Kr Dist', 'Kr Prox','Kr BothSep','Kr 2xDist','Kr 2xProx');
title('Mean spot correlation all constructs')
xlabel('Mean AP position')
ylabel('Mean spot fluorescence correlation')
ylim([-0.3 1]);
xlim([0.3 0.75]);
            
figure 
%NEmbryos=[];
Data=[];

% for qq=1:length(ConstructList)
%     Data(qq)= LoadMS2SetsCS(ConstructList{qq});
% end
BothEmbryos=length(BothData);
for ii=1:BothEmbryos
    q(ii)=plot(APBins,MeanAPCorr(ii,:,1),'o','LineWidth',2,'Color','r');
    hold on
end
DistEmbryos=length(DistData);
for qq=1:DistEmbryos
    w(qq)=plot(APBins,MeanAPCorr(qq,:,2),'o','LineWidth',2,'Color','g');
end
ProxEmbryos=length(ProxData);
for zz=1:ProxEmbryos
    z(zz)=plot(APBins,MeanAPCorr(zz,:,3),'o','LineWidth',2,'Color','b');
end
BothSepEmbryos=length(BothSepData);
for hh=1:BothSepEmbryos
    h(hh)=plot(APBins,MeanAPCorr(hh,:,4), 'o', 'LineWidth', 2, 'Color', 'y');
end
DoubDistEmbryos=length(DoubDistData);
for rr=1:DoubDistEmbryos
    r(rr)=plot(APBins,MeanAPCorr(rr,:,5),'o','LineWidth',2,'Color','c');
end
DoubProxEmbryos=length(DoubProxData);
for bb=1:DoubProxEmbryos
    b(bb)=plot(APBins,MeanAPCorr(bb,:,5),'o','LineWidth',2,'Color','m');
end
legend([q(1), w(1), z(1), h(1),r(1),b(1)],{'Both','Dist','Prox','BothSep','2xDist','2xProx'});
 xlabel('Mean AP position'); ylabel('Avg correlation of spot fluorescence');
 title('Mean spot correlation all embryos');
 ylim([0 1]);

 figure %plot the mean of correlation across the whole embryo for each embryo color-coded by construct 
 for uu=1:BothEmbryos 
     u(uu)=plot(1,MeanAllCorr(uu,1),'o','Color','r');
     hold on 
 end
 for jj=1:DistEmbryos
     t(jj)=plot(2,MeanAllCorr(jj,2),'o','Color','g');
 end
 for hh=1:ProxEmbryos
     b(hh)=plot(3,MeanAllCorr(hh,3),'o','Color','b');
 end
 for h=1:BothSepEmbryos
     as(h)=plot(4,MeanAllCorr(h,4),'o', 'Color', 'y');
 end
 for r=1:DoubDistEmbryos
     ra(r)=plot(5,MeanAllCorr(r,5),'o','Color','c');
 end
 for v=1:DoubProxEmbryos
     va(v)=plot(6,MeanAllCorr(v,6),'o','color','m');
 end
 legend([u(1),t(1),b(1),as(1),ra(1),va(1)],{'Both','Dist','Prox','BothSep','2xDist','2xProx'},'Location','best');
 ylim([0 1]); xlim([0 7]);
 ylabel('Mean spot correlation');
 
 title('Mean spot correlation across whole embryo');

            
            
        
        